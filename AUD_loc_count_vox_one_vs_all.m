%binarize maps for each subject of his activation > others controls for
%normal > reverse speech; then it intersects it with language mask,
%individual reading mask, and count number of voxel for each + ratio over
%total num of voxels.
clear;clc;

D       = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/one_vs_all_s8';
cd(D)
S       = dir('*vs_con_age_covariate');
Syn     = S(~cellfun(@isempty,(regexp({S.name},'Sujet'))));
Con     = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
mask    = ismember({S.name}, {'.', '..'});
S(mask) = '';
scans   = {};

%% threshold the maps
pvalue = 0.01;
for k = 1 : numel(S)
    cd(fullfile(D,S(k).name,'words + pseudowords + numbers + normal_speech'));
    con_diff = 'spmT_0003.nii';
    con_header      = spm_vol(con_diff);
    con_vol         = spm_read_vols(con_header);
    
    load('SPM.mat');
    criticalt = tinv(1-pvalue,SPM.xX.erdf);
    con_vol(con_vol<criticalt) = 0;
    con_vol(con_vol>=criticalt) = 1;
    nvox_thr(k) = length(find(con_vol>0));
    
    con_header_copy = con_header;
    con_header_copy.fname = 'spmT_0003_threshold_10-2.nii';
    spm_write_vol(con_header_copy,con_vol);
    
end

%% add maps in both groups
scans1 = {};
scans2 = {};

% for syn
for k = 1 : numel(Syn)
    cd(fullfile(D,Syn(k).name,'words + pseudowords + numbers + normal_speech'));
    thr_con_diff = 'spmT_0003_threshold_10-2.nii';
    scans1 = [scans1; fullfile(D,Syn(k).name,'words + pseudowords + numbers + normal_speech', thr_con_diff)];
end
for k = 1 : numel(Con)
    cd(fullfile(D,Con(k).name,'words + pseudowords + numbers + normal_speech'));
    thr_con_diff = 'spmT_0003_threshold_10-2.nii';
    scans2 = [scans2; fullfile(D,Con(k).name,'words + pseudowords + numbers + normal_speech', thr_con_diff)];
end

total = '';
for k = 1:length(scans1)
    map = sprintf('%s%d','i',k);
    if k == 1
        total = map;
    else
        total = [total '+' map];
    end
end
total = ['(' total ')'];

outdir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/one_vs_all_s8';
if ~isdir(outdir)
    mkdir(outdir)
end

matlabbatch{1}.spm.util.imcalc.input = scans1;
matlabbatch{1}.spm.util.imcalc.output = sprintf('mean_speech>baseline_syn_10-2');
matlabbatch{1}.spm.util.imcalc.outdir = {outdir};
matlabbatch{1}.spm.util.imcalc.expression = [total '/' num2str(length(scans1))];
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = -7;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

matlabbatch{2}.spm.util.imcalc.input = scans2;
matlabbatch{2}.spm.util.imcalc.output = sprintf('mean_speech>baseline_con_10-2');
matlabbatch{2}.spm.util.imcalc.outdir = {outdir};
matlabbatch{2}.spm.util.imcalc.expression = [total '/' num2str(length(scans2))];
matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = -7;
matlabbatch{2}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run', matlabbatch);
    



% 
% %% count the number of voxels for each subject difference map
% % first, need to binarize maps
% 
% for k = 1 : numel(S)
%     cd(fullfile(D,S(k).name,'phonology'));
%     diff_activation = dir('*others_10-3_fwe510-2*');
%     matlabbatch{1}.spm.util.imcalc.input = {diff_activation.name};
%     output = sprintf('binarized_%s.nii', diff_activation.name(1:end-4));
%     matlabbatch{1}.spm.util.imcalc.output = output;
%     matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(D,S(k).name,'phonology')};
%     matlabbatch{1}.spm.util.imcalc.expression = 'i1>0';
%     matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
%     matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
%     matlabbatch{1}.spm.util.imcalc.options.mask = 0;
%     matlabbatch{1}.spm.util.imcalc.options.interp = 0;
%     matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
%     spm_jobman('run', matlabbatch);
% end
% 
% %%
% for k = 1 : numel(S)
%     cd(fullfile(D,S(k).name,'phonology'));
%     diff_activation = dir('binarized*others_10-3_fwe510-2*');
%     image = diff_activation.name;
%     [status vox] = num_voxel_image(image);
%     num_vox(k) = vox;
% end
% 
% for num = 1:length(num_vox)
%     if num_vox(num) == 463623
%         num_vox(num)=0;
%     end
% end
