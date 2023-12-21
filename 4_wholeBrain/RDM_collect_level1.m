function RDM_collect_level1
% Harrison Ritz 2021


pts = 8004:8031;


root_dir = "/users/hritz/data/mri-data/RDM2"
modelName = ''
save_dir = fullfile(root_dir, 'spm-data', 'collect', modelName)



fileNames = {'spmT_0002.nii', 'spmT_0004.nii', 'spmT_0006.nii', 'spmT_0008.nii', 'spmT_0010.nii'}


% clear folder
if exist(save_dir, 'dir')
    rmdir(save_dir, 's'); % remove existing results folder
end
mkdir(save_dir);


for pp = 1:length(pts)
    
    mkdir(fullfile(save_dir, sprintf('sub-%d', pts(pp))));
    
    % get anat
    endPath = fullfile(save_dir, sprintf('sub-%d', pts(pp)));
    startPath = fullfile(root_dir, 'spm-data', sprintf('sub-%d', pts(pp)), 'anat', '*desc-preproc_T1w*');
    
    copyfile(startPath, endPath);
    
    % get unprocessed anaot
    endPath = fullfile(save_dir, sprintf('sub-%d', pts(pp)));
    startPath = fullfile(root_dir, 'spm-data', sprintf('sub-%d', pts(pp)), 'anat', '*desc-preproc_T1w*');
    
    copyfile(startPath, endPath);
    
    
    % get func
    for ff = 1:length(fileNames)
        
        endPath = fullfile(save_dir, sprintf('sub-%d', pts(pp)), fileNames{ff});
        startPath = fullfile(root_dir, 'spm-data', sprintf('sub-%d', pts(pp)),'level 1', modelName, fileNames{ff});
        
        % copy
        copyfile(startPath, endPath);
        
    end
    
    
end




end
