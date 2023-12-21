function [G,theta,u,l,numIter]=pcm_NR_diag(Y,Z,varargin)
% function [G,theta,u,l,numIter]=pcm_NR_diag(y,Z,varargin)
% Estimate random-effects variance component coefficients using
% Newton-Raphson gradient descent.
% For this form of the algorithm, the G-matrix components need to be diagnoal 
% But otherwise can be arbitrarily constrained 
% with diag(G)=sum(Gc*hc))
% The diagonal structure can be specified by giving either Gd or Gc
% 
% For models with non-diagnoal Gc matrices, use pcm_NR.  
% In many applications pcm_NR_diag should be faster. 
% 
% Estimates the variance coefficients of the model described in:
% Diedrichsen, Ridgway, Friston & Wiestler (2011).
%
% y_n = X b_n + Z u_n + e,
%         u_n ~ N(0, G)
%           G = sum (exp(theta_c) * Gc)
%         e_n ~ N(0,S exp(theta_e)) 
%
% Y: N x P observations
% Z: N x Q random effects matrix
% X: N x L fixed effects matrix
%
% N: numbers of observations (trials)
% P: numbers of experimental units (voxels)
% Q: number of random effects
%
% VARARGIN:
%   'numIter'     : Maximal number of iterations
%   'Gc'           : (QxQxH) tensor of variance components
%                    matrix G = sum(h_m Gc(:,:,i)), with m=1...H
%   'Gd'           : Q x H matrix with each column referring the the diagonal of the
%                    variance component
%   'h0'           : Starting values for the parameters (Hx1 vector)
%                    (otherwise uses starting guess based on Laird & Lange)
%   'TolL'         : Tolerance of the likelihood (l-l'), where l' is
%                    projected likelihood
%   'X'            : Fixed effects matrix that will be removed from y
%                    In this case, ReML will be used for the estimation of
%                    G.
%   'S':            Explicit noise covariance matrix structure matrix. The For speed, 
%                   this is a cell array that contains
%                   S.S:     Structure of noise 
%                   S.invS:  inverse of the noise covariance matrix 
%   'HessReg':      Regulariser on the Hessian to increase the stability of
%                   the fit (set to 1/256)
%
% OUTPUT:
%   G       : variance-covariance matrix
%   theta   : (log-)parameters
%   u       : hidden patterns components
%   l       : log-likelihood
%   numIter : Numer of iterations 
%
% Examples:
% v.3.0:
%
% Copyright 2017 Joern Diedrichsen, joern.diedrichsen@googlemail.com 

% Defaults
%--------------------------------------------------------------------------
ac      = [];                      % Which terms do I want to include?
numIter = 32;                 % Maximal number of iterations
Gc      = [];
Gd      = [];
h0      = [];
low     = -16;
thres   = 1e-2;
X       = [];
S       = [];
HessReg = 1/256;                  % Regulariser on Hessian
i  = 0; 
% Variable argument otions
%--------------------------------------------------------------------------
vararginoptions(varargin, ...
    {'Gc','Gd','meanS','h0','ac','thres','X','numIter','S','HessReg'});

% check input size
%--------------------------------------------------------------------------
Y=double(Y);
Z=double(Z);
[N,P]  =  size(Y);
[N2,Q] =  size(Z);
if N2  ~= N
    error('Mismatched numbers of rows in data (%d) and design (%d)', N, N2)
end

% Intialize the Model structure
%--------------------------------------------------------------------------
if (isempty(Gd) && ~isempty(Gc))
    H     =  size(Gc,3)+1;       % Number of Hyperparameters (+ 1 for noise term)
    for i = 1:H-1
        Gd(:,i)=diag(Gc(:,:,i));
    end;
elseif (isempty(Gc) && ~isempty(Gd))
    H    = size(Gd,2)+1;
    for h = 1:H-1
        Gc(:,:,h)=diag(Gd(:,h));
    end;
elseif (~isempty(Gc) && ~isempty(Gd))
    error('Only Gd or Gc option should be given');
else
    error('Either Gd or Gc option should be given');
end;

YY   =  Y*Y';                        % This is Suffient stats 1 (S1)
trYY =  sum(diag(YY));

% Figure out which terms to include into the model (set others to -32)
%--------------------------------------------------------------------------
if (isempty(ac))
    as    = 1:H;
    nas   = [];
else
    as    = find(ac);
    nas   = find(ac==0);
end;

% initialise h
%--------------------------------------------------------------------------
if (isempty(h0))
    rs     = Y-X*pinv(X)*Y; 
    h      =  ones(H,1)*low;
    u      =  pinv(Z)*rs;
    h(H,1) =  (trYY-traceAB(u'*Z',rs))/(P*N);
    D      =  u*u'/P-h(H)*pinv(Z'*Z)/P;
    hD     =  diag(D);
    % Starting values for constrained estimates
    h(1:H-1,1) = real(log((Gd'*Gd)\(Gd'*hD(:))));
    h(h<low) = low+1;
else
    h      = h0;
end;
h(nas) = low;

% Initialize Interations and hyperpriors
%--------------------------------------------------------------------------
dF    = Inf;
dFdh  = zeros(H,1);
dFdhh = zeros(H,H);
HessReg = HessReg*speye(H,H);          % Regulariser on Hessian 


for k = 1:numIter
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    G     = zeros(Q,Q);
    for i = as(1:end-1);
        G = G + Gc(:,:,i)*exp(h(i));
    end
    
    % Find the inverse of V - while dropping the zero dimensions in G
    [u,s] = eig(full(G));
    dS    = diag(s);
    idx   = dS>eps;
    Zu     = Z*u(:,idx);
    if (isempty(S))
        iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(h(H))+Zu'*Zu)*Zu')./exp(h(H)); % Matrix inversion lemma
    else
        iV    = (S.invS-S.invS*Zu/(diag(1./dS(idx))*exp(h(H))+Zu'*S.invS*Zu)*Zu'*S.invS)./exp(h(H)); % Matrix inversion lemma
    end; 

    if (~isempty(X))
        iVX   = iV * X;
        iVr   = iV - iVX*((X'*iVX)\iVX');  % Correction for the fixed effects
    else
        iVr   = iV;
    end;
    
    % ReML estimate of hyperparameters
    
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    A     = iVr*Z;
    B     = YY*A/P;
    dFdh(1:H-1) = -P/2*(sum(A.*Z)-sum(A.*B)) * Gd;
    if (isempty(S)) 
        dFdh(H) = -P/2*traceABtrans(iVr,(speye(N)-YY*iVr/P));
    else 
        dFdh(H) = -P/2*traceABtrans(iVr*S.S,(speye(N)-YY*iVr/P));
    end; 
    
    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    dFdhh(1:H-1,1:H-1)  = -P/2*Gd'*((Z'*A).*(Z'*A))*Gd; % Short form
    if (isempty(S)) 
        dFdhh(H,H)          = -P/2*traceABtrans(iVr,iVr);
        dFdhh(1:H-1,H)      = -P/2*Gd'*sum((Z'*iVr).*A',2); %  -P/2*Gd'*diag(Z'*iVr*A);
    else
        iVrS = iVrS; 
        dFdhh(H,H)          = -P/2*traceABtrans(iVrS,iVrS);
        dFdhh(1:H-1,H)      = -P/2*Gd'*sum((Z'*iVrS).*A',2);
    end; 
    dFdhh(H,1:H-1)      = dFdhh(1:H-1,H)';
    
    
    % modulate
    %----------------------------------------------------------------------
    dFdh  = dFdh.*exp(h);
    dFdhh = dFdhh.*(exp(h)*exp(h)');
        
    % Regularise second moment by a little bit 
    dFdhh = dFdhh - HessReg;

    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    % dh    = spm_dx(dFdhh(as,as),dFdh(as),{t}); % This is actually much
    % less accurate it seems than
    dh   = -dFdhh(as,as)\dFdh(as);
    h(as) = h(as) + dh;
    
    
    % predicted change in F - increase regularisation if increasing
    %----------------------------------------------------------------------
    pF    = dFdh(as)'*dh;
    dF    = pF;
    
    % convergence
    %----------------------------------------------------------------------
    % fprintf('%s %-23d: %10s%e [%+3.2f]\n','  ReML Iteration',k,'...',full(dF),t);
    if dF < thres
        break;
    else
        % eliminate redundant components (automatic selection)
        %------------------------------------------------------------------
        as  = find(h > low);
        h(h<low)=low;
        as  = as(:)';
    end
    hh(:,k)=h;
end

% return exp(h) and rescale
%--------------------------------------------------------------------------
theta  = h;
G     = zeros(Q);
for i = as(1:end-1);
    G = G + full(Gc(:,:,i))*exp(theta(i));
end

% Find the inverse of V - while dropping the zero dimensions in G
% V = Z*G*Z' + I sigma 
[u,s] = eig(G);
dS    = diag(s);
idx   = dS>eps;
Zu     = Z*u(:,idx);
if (isempty(S))
    iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(theta(H))+Zu'*Zu)*Zu')./exp(theta(H)); % Matrix inversion lemma
else 
    iV    = (S.invS-S.invS*Zu/(diag(1./dS(idx))*exp(theta(H))+Zu'*S.invS*Zu)*Zu'*S.invS)./exp(theta(H)); % Matrix inversion lemma
end; 

if (~isempty(X))
    iVX   = iV * X;
    iVr   = iV - iVX*((X'*iVX)\iVX');  % Correction for the fixed effects
else
    iVr   = iV;
end;

if nargout > 2
    u=G*Z'*(iV)*Y;
end;
if ~all(isreal(u(:)));
    keyboard;
end;

if nargout > 3
    ldet  = -2*sum(log(diag(chol(iV))));        % Safe computation of the log determinant (V) Thanks to code from D. lu
    l     = -P/2*(ldet)-0.5*traceABtrans(iVr,YY);
    if (~isempty(X)) % Correct for ReML estimates
        l = l - P*sum(log(diag(chol(X'*iV*X))));  % - 0.5 * log(det());
    end;
end;
