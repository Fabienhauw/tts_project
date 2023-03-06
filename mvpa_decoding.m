%% Decoding
clc; clear;
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/decoding_toolbox_v3.997'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12'))
addpath(genpath('/network/lustre/iss02/apps/software/scit/AFNI'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/abin'))

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
S = dir(D);
mask = ismember({S.name}, {'.', '..'});
S(mask) = [];
wd = pwd;
radius = 10;

clear matlabbatch
i=1;

label1 = 'Words';
label2 =  'Numbers';

% label1 = 'ScrSpeech';
% label2 =  'NormalSpeech';

% label1 = 'PW';
% label2 =  'Words';

a = 1; b = numel(S);
% a = 45; b = 45;

for k = a:b
    % Conditions: Motor, NormalSpeech, Numbers, Odds, PW, ScrSpeech, Words    
    matlabbatch{i}.spm.tools.tdt.decod.subj.dir = {sprintf('%s/%s/Aud/loc/mvpa_without_resting_%dmm', D, S(k).name, radius)};
%     matlabbatch{i}.spm.tools.tdt.decod.subj.dir = {fullfile(D, S(k).name, 'Aud/loc/mvpa')};
    matlabbatch{i}.spm.tools.tdt.decod.subj.conds.cond1 = label1;
    matlabbatch{i}.spm.tools.tdt.decod.subj.conds.cond2 = label2;
    matlabbatch{i}.spm.tools.tdt.decod.options.nrun = 8;
    matlabbatch{i}.spm.tools.tdt.decod.options.anal.searchlight.unit = 'mm';
    matlabbatch{i}.spm.tools.tdt.decod.options.anal.searchlight.rad = radius;
    matlabbatch{i}.spm.tools.tdt.decod.options.anal.searchlight.mask = {sprintf('%s/%s/Aud/loc/mvpa_without_resting_%dmm/mask.nii', D, S(k).name, radius)};    
    matlabbatch{i}.spm.tools.tdt.decod.options.meth = 'kernel';
    matlabbatch{i}.spm.tools.tdt.decod.options.analysis = 'accuracy';
    matlabbatch{i}.spm.tools.tdt.decod.options.display = 'no';
    
    i = i + 1;
    
end

% spm_jobman('run', matlabbatch)


%% .....................................................
cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'mvpa_decod';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);
















%% Permutations

clear matlabbatch
i=1;

results_folder = sprintf('results_%s_vs_%s', label1, label2);
results_fold_inv = sprintf('results_%s_vs_%s', label2, label1);


% Conditions: Motor, NormalSpeech, Numbers, Odds, PW, ScrSpeech, Words    

for k = a:b
%     mvpa_res_dir = fullfile(D, S(k).name, 'Aud/loc/mvpa');
    mvpa_res_dir = sprintf('%s/%s/Aud/loc/mvpa%dmm', D, S(k).name, radius);
    cd(mvpa_res_dir)
    files = dir;
    dirFlags = [files.isdir];
    subDirs = files(dirFlags);
    for dirname = 1 : length(subDirs)
        if ~isempty(regexp(subDirs(dirname).name,results_fold_inv))
            tmp_results_folder=results_fold_inv;
            fprintf('This analysis already exists, with %s labelled Condition1 and conversely (same decoding). Overwriting in the corresponding folder. \n', label2)
        elseif ~isempty(regexp(subDirs(dirname).name,results_folder))
            tmp_results_folder = results_folder;
        end
    end
    
    matlabbatch{i}.spm.tools.tdt.perm.res_cfg(1) = {fullfile(mvpa_res_dir, tmp_results_folder, 'res_cfg.mat')};
    matlabbatch{i}.spm.tools.tdt.perm.res_dir = {fullfile(mvpa_res_dir, tmp_results_folder)};
    matlabbatch{i}.spm.tools.tdt.perm.suffix_res_dir = 'perm';
    matlabbatch{i}.spm.tools.tdt.perm.options.all_perm = 'all';
    
    i = i + 1;
end

% spm_jobman('interactive', matlabbatch)

%% ....................................
cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '12:00:00'; % usually takes 7-8h on cluster for 124 permutations x 8-fold classifier
par.jobname  = 'mvpa_perm';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

%% T-tests
clear matlabbatch
i=1;

results_folder = sprintf('results_%s_vs_%s', label1, label2);
results_fold_inv = sprintf('results_%s_vs_%s', label2, label1);

% Conditions: Motor, NormalSpeech, Numbers, Odds, PW, ScrSpeech, Words    

for k = a:b
%     mvpa_res_dir = fullfile(D, S(k).name, 'Aud/loc/mvpa');
    mvpa_res_dir = sprintf('%s/%s/Aud/loc/mvpa%dmm', D, S(k).name, radius);
    cd(mvpa_res_dir)
    files = dir;
    dirFlags = [files.isdir];
    subDirs = files(dirFlags);
    for dirname = 1 : length(subDirs)
        if ~isempty(regexp(subDirs(dirname).name,results_fold_inv))
            tmp_results_folder = results_fold_inv;
            fprintf('This analysis already exists, with %s labelled Condition1 and conversely (same decoding). Overwriting in the corresponding folder. \n', label2)
        elseif ~isempty(regexp(subDirs(dirname).name,results_folder))
            tmp_results_folder = results_folder;
        end
    end
    
    matlabbatch{i}.spm.tools.tdt.ttest.dcdg_dir(1) = {fullfile(mvpa_res_dir, tmp_results_folder)};
    matlabbatch{i}.spm.tools.tdt.ttest.perm_dir(1) = {fullfile(mvpa_res_dir, tmp_results_folder, 'perm')};
    matlabbatch{i}.spm.tools.tdt.ttest.options = 'right';
    
    i = i + 1;
end

% spm_jobman('interactive', matlabbatch)

%% ....................................
cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')
spm_jobman('run', matlabbatch)

% par.run = 0;
% par.sge = 1;
% par.sge_queu = 'normal,bigmem';
% par.pct = 1;
% par.walltime = '24:00:00'; % probably 24h is enough
% par.jobname  = 'mvpa_ttest';
% %%%%%%%% this line below to comment to avoid re estimating
% %%%%%%%% models
% 
% job_ending_rountines(matlabbatch, [], par);

%% NIFTI Writing
clear matlabbatch
i=1;

results_folder = sprintf('results_%s_vs_%s', label1, label2);
results_fold_inv = sprintf('results_%s_vs_%s', label2, label1);

% Conditions: Motor, NormalSpeech, Numbers, Odds, PW, ScrSpeech, Words    

for k = a:b
%     mvpa_res_dir = fullfile(D, S(k).name, 'Aud/loc/mvpa');
    mvpa_res_dir = sprintf('%s/%s/Aud/loc/mvpa%dmm', D, S(k).name, radius);
    cd(mvpa_res_dir)
    files = dir;
    dirFlags = [files.isdir];
    subDirs = files(dirFlags);
    for dirname = 1 : length(subDirs)
        if ~isempty(regexp(subDirs(dirname).name,results_fold_inv))
            tmp_results_folder = results_fold_inv;
            fprintf('This analysis already exists, with %s labelled Condition1 and conversely (same decoding). Overwriting in the corresponding folder. \n', label2)
        elseif ~isempty(regexp(subDirs(dirname).name,results_folder))
            tmp_results_folder = results_folder;
        end
    end
    
    matlabbatch{i}.spm.tools.tdt.nft.dir(1) = {fullfile(mvpa_res_dir, tmp_results_folder)};
    matlabbatch{i}.spm.tools.tdt.nft.pmat(1:2,1) = {
        fullfile(mvpa_res_dir, tmp_results_folder, 'perm/p_value.mat')
        fullfile(mvpa_res_dir, tmp_results_folder, 'perm/1-p_value.mat')
        };
    
    i = i + 1;
    
end

% spm_jobman('interactive', matlabbatch)

%% ....................................
cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

spm_jobman('run', matlabbatch)

% par.run = 0;
% par.sge = 1;
% par.sge_queu = 'normal,bigmem';
% par.pct = 1;
% par.walltime = '24:00:00'; % probably 24h is enough
% par.jobname  = 'mvpa_perm';
% %%%%%%%% this line below to comment to avoid re estimating
% %%%%%%%% models
% 
% job_ending_rountines(matlabbatch, [], par);