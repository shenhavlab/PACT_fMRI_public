function flip_nii(path, isFilt)

doFilt = 0;
if nargin >1
    doFilt = isFilt;
end

if doFilt % loop over filt

    dr = dir([path, '*.nii']);

    for dd = 1:length(dr)

        path = fullfile(dr(dd).folder, dr(dd).name);

        vo = spm_vol(path);
        newvo      = struct(...
            'fname',    vo.fname,...
            'dim',      vo(1).dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      vo(1).mat,...
            'n',        [1 1],...
            'descrip',  [vo.descrip, ' - LR flipped']);

        spm_write_vol(newvo, flip(spm_read_vols(vo)));

    end

else % just do single image

    vo = spm_vol(path);
    newvo      = struct(...
        'fname',    [vo.fname(1:end-4), '_flip.nii'],...
        'dim',      vo(1).dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      vo(1).mat,...
        'n',        [1 1],...
        'descrip',  [vo.descrip, ' - LR flipped']);

    spm_write_vol(newvo, flip(spm_read_vols(vo)));

end

end