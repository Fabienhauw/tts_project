%binarize maps for each subject of his activation > others controls for
%normal > reverse speech; then it intersects it with language mask,
%individual reading mask, and count number of voxel for each + ratio over
%total num of voxels.
clear;clc;
D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/one_vs_all/';
cd(D)
S=dir('*vs_con_age_covariate');
mask = ismember({S.name}, {'.', '..'});
S(mask)='';
scans = {};

%% count the number of voxels for each subject difference map
% first, need to binarize maps

for k = 1 : numel(S)
    cd(fullfile(D,S(k).name,'phonology'));
    diff_activation = dir('*others_10-3_fwe510-2*');
    matlabbatch{1}.spm.util.imcalc.input = {diff_activation.name};
    output = sprintf('binarized_%s.nii', diff_activation.name(1:end-4));
    matlabbatch{1}.spm.util.imcalc.output = output;
    matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(D,S(k).name,'phonology')};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1>0';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run', matlabbatch);
end

%%
for k = 1 : numel(S)
    cd(fullfile(D,S(k).name,'phonology'));
    diff_activation = dir('binarized*others_10-3_fwe510-2*');
    image = diff_activation.name;
    [status vox] = num_voxel_image(image);
    num_vox(k) = vox;
end

for num = 1:length(num_vox)
    if num_vox(num) == 463623
        num_vox(num)=0;
    end
end
