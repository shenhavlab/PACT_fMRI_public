function pcm_vararginoptions(options,allowed_vars,allowed_flags);
% function pcm_vararginoptions(options,allowed_vars,allowed_flags);
% Deals with variable argument in   
% INPUTS
%   options: cell array of a argument list passed to a function
%   allowed_vars: Variables that can be set 
%   allowed_flags: Flags that can be set 
%  vararginoptions assigns the value of the option to a variable with the
%  name option (in called workspace).
%  Flags are set to one:
% EXAMPLE:
%   the option-string 'var1',4,'var2',10,'flag'
%   causes the var1 and var2 to be set to 4 and 10 and flag to 1
%   if allowedvars are not given, all variables are allowed
% Joern Diedrichsen 
% v1.0 9/13/05
% v2.0 7/06/07 Extentension to understand struct-list
checkflags=1;
checkvars=1;
if nargin<2
    checkvars=0;
end;
if nargin<3
    checkflags=0;
end;

c=1;
while c<=length(options)
    a=[];
    if ~ischar(options{c})
        error(sprintf('Options must be strings on argument %d',c));
    end;
    if checkflags
        a=find(strcmp(options{c},allowed_flags));
    end;
    if ~isempty(a)
        assignin('caller',options{c},1);
        c=c+1;
    else
        if checkvars
            a=find(strcmp(options{c},allowed_vars));
            if (isempty(a))
                error(['unknown option:' options{c}]);
            end;
        end;
        if (c==length(options))
            error(sprintf('Option %s must be followed by a argument',options{c}));
        end;
        p=strfind(options{c},'.');
        if (~isempty(p))
            structname=options{c}(1:p-1); 
            fieldname=options{c}(p+1:end); 
            if (evalin('caller',['exist(''' structname ''')']));
                S=evalin('caller',structname);
            else 
                S=[]; 
            end; 
            S.(fieldname)=options{c+1};
            assignin('caller',structname,S);
        else         
            assignin('caller',options{c},options{c+1});
        end; 
        c=c+2;
    end;
end;