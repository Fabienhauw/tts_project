% to run decoding between 2 conditions, within regions of interest, defined
% with spheres around peak coordinates of activation in normal > scr speech
% in auditive run
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

label1 = 'ScrSpeech';
label2 =  'NormalSpeech';

% label1 = 'PW';
% label2 =  'Words';

a = 1; b = numel(S);
% a = 45; b = 45;

for k = a : b
    cd(fullfile(D, S(k).name, 'Aud/loc/stats_s5'))
    ROIs = dir('spmT_0011_10mm_sphere*');

    for tmp_roi = 1 : length(ROIs)
        % Conditions: Motor, NormalSpeech, Numbers, Odds, PW, ScrSpeech, Words
        beta_loc = fullfile(D, S(k).name, sprintf('Aud/loc/mvpa%dmm_rois', radius));
        dir_roi = ROIs(tmp_roi).name(1:end-4);
        res_dir = fullfile(beta_loc, dir_roi);
        matlabbatch{i}.spm.tools.tdt.decod.subj.dir = {beta_loc};
        matlabbatch{i}.spm.tools.tdt.decod.subj.res_dir = {res_dir};
        matlabbatch{i}.spm.tools.tdt.decod.subj.conds.cond1 = label1;
        matlabbatch{i}.spm.tools.tdt.decod.subj.conds.cond2 = label2;
        matlabbatch{i}.spm.tools.tdt.decod.options.nrun = 8;
        matlabbatch{i}.spm.tools.tdt.decod.options.anal.ROI.mask_roi = {fullfile(D, S(k).name, 'Aud/loc/stats_s5', ROIs(tmp_roi).name)};
        matlabbatch{i}.spm.tools.tdt.decod.options.meth = 'kernel';
        matlabbatch{i}.spm.tools.tdt.decod.options.analysis = 'accuracy';
        matlabbatch{i}.spm.tools.tdt.decod.options.display = 'no';
        matlabbatch{i}.spm.tools.tdt.decod.options.overwrite = 1;
        i = i + 1;
    end
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

%%
clear matlabbatch
i=1;
for k = a : b
    cd(fullfile(D, S(k).name, 'Aud/loc/stats_s5'))
%     ROIs = dir('spmT_0011_best_vox*');
    ROIs = dir('spmT_0011_best_vox*');

    for tmp_roi = 1 : length(ROIs)
        % Conditions: Motor, NormalSpeech, Numbers, Odds, PW, ScrSpeech, Words
        beta_loc = fullfile(D, S(k).name, sprintf('Aud/loc/mvpa%dmm_rois', radius));
        dir_roi = ROIs(tmp_roi).name(1:end-4);
        res_dir = fullfile(beta_loc, dir_roi);
        matlabbatch{i}.spm.tools.tdt.decod.subj.dir = {beta_loc};
        matlabbatch{i}.spm.tools.tdt.decod.subj.res_dir = {res_dir};
        matlabbatch{i}.spm.tools.tdt.decod.subj.conds.cond1 = label1;
        matlabbatch{i}.spm.tools.tdt.decod.subj.conds.cond2 = label2;
        matlabbatch{i}.spm.tools.tdt.decod.options.nrun = 8;
        matlabbatch{i}.spm.tools.tdt.decod.options.anal.ROI.mask_roi = {fullfile(D, S(k).name, 'Aud/loc/stats_s5', ROIs(tmp_roi).name)};
        matlabbatch{i}.spm.tools.tdt.decod.options.meth = 'kernel';
        matlabbatch{i}.spm.tools.tdt.decod.options.analysis = 'accuracy';
        matlabbatch{i}.spm.tools.tdt.decod.options.display = 'no';
        matlabbatch{i}.spm.tools.tdt.decod.options.overwrite =1;
        i = i + 1;
    end
end

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

