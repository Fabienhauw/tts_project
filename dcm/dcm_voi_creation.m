%%

clear;
addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))

wd = pwd;
all_con = [16];
roi_sph_rad = 4;
F_con = 29;


res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/MVPA/Aud/loc/syn_vs_con_rh_s5_10mm/roi_decoding_comparisons';
% res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/MVPA/Aud/loc/syn_vs_con_rh_s5';

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet'))));
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17


gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% gaucher_appar = {'Control02|Control04|Control07|Control17|Sujet'};
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));

mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con = S;
S_con(mask_gauch_con) = [];

S_syn = S;
S_syn(mask_gauch) = [];

S_effect = [S_syn ; S_con];
% S_effect = S_effect(18);

for j = 1 : size(S,1)
    if ~isempty(find(~cellfun(@isempty,(regexp({S_effect.name},S(j).name)))))
        mask_cov(j,1) = 1;
    else
        mask_cov(j,1) = 0;
    end
end

% coord from TTS network (speech > baseline, syn > con) except: common_vwfa
% from con sent>scrambled in syn + con and smg common from speech >
% baseline for syn + con
% ROIs_names = {'SMG',    'VWFA',      'lIPS',     'lprecent',   'MFG',  'lpSTG',   'laSTG',  'rpSTG',   'raSTG',   'lSTGm',    'AngG',    'SMG_common'};
% ROIs_coord = [-48 -44 23; -45 -51 -10; -40 -41 46; -50 -16 50; -50 6 53; -70 -28 3; -60 12 -7; 50 -24 16; 65 4 0; -45 -24 8; -30 -71 33; -48 -41 16]; 

ROIs_names = { 'VWFA_common', 'lSTG_common'};
ROIs_coord = [-42 -38 -20; -58 -11 0];

i=1;
for tmp_con = 1 : length(all_con)
    select_con = all_con(tmp_con);
    for k = 1 : numel(S_effect)
        for tmp_ROI = 1 : numel(ROIs_names)
            %% batch for VOI, resulting from all subj, speech>baseline
            matlabbatch{i}.spm.util.voi.spmmat = {fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5_without_resting/SPM.mat')}; % use this .mat to extract time series.
            matlabbatch{i}.spm.util.voi.adjust = F_con;
            matlabbatch{i}.spm.util.voi.session = 1;
            matlabbatch{i}.spm.util.voi.name = sprintf('%s_%d_%d_%d_%dmm_sph_DCM_ROI_adapted_to_aud_con%d_adj_eoi5', ROIs_names{tmp_ROI}, ROIs_coord(tmp_ROI, :), roi_sph_rad, select_con);
            matlabbatch{i}.spm.util.voi.roi{1}.spm.spmmat = {fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5_without_resting/SPM.mat')}; % use this .mat to deduct the highest activated voxel.
            matlabbatch{i}.spm.util.voi.roi{1}.spm.contrast = select_con;
            matlabbatch{i}.spm.util.voi.roi{1}.spm.conjunction = 1;
            matlabbatch{i}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
            matlabbatch{i}.spm.util.voi.roi{1}.spm.thresh = 1;
            matlabbatch{i}.spm.util.voi.roi{1}.spm.extent = 0;
            matlabbatch{i}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
            matlabbatch{i}.spm.util.voi.roi{2}.sphere.centre = ROIs_coord(tmp_ROI,:);
            matlabbatch{i}.spm.util.voi.roi{2}.sphere.radius = 6;
            matlabbatch{i}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
            matlabbatch{i}.spm.util.voi.roi{3}.sphere.centre = ROIs_coord(tmp_ROI,:);
            matlabbatch{i}.spm.util.voi.roi{3}.sphere.radius = roi_sph_rad;
            matlabbatch{i}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
            matlabbatch{i}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';
            matlabbatch{i}.spm.util.voi.expression = 'i1&i3';
            i=i+1;
        end
    end
end

%%
% spm_jobman('run', matlabbatch)

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'dcm_voi';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

%% for the STS resulting from the mvpa analysis, spherical ROI:
i=1;
clear matlabbatch
for tmp_con = 1 : length(all_con)
    select_con = all_con(tmp_con);
    for k = 1 : numel(S_effect)
        %% batch for VOI, resulting from all subj, norm>scr speech
        matlabbatch{i}.spm.util.voi.spmmat = {fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5_without_resting/SPM.mat')}; % use this .mat to extract time series.
        matlabbatch{i}.spm.util.voi.adjust = F_con;
        matlabbatch{i}.spm.util.voi.session = 1;
        matlabbatch{i}.spm.util.voi.name = sprintf('lSTS_-65_-41_6_%dmm_sph_DCM_ROI_adapted_to_aud_con%d_adj_eoi5', roi_sph_rad, select_con);
        matlabbatch{i}.spm.util.voi.roi{1}.sphere.centre = [-65 -41 6];
        matlabbatch{i}.spm.util.voi.roi{1}.sphere.radius = roi_sph_rad;
        matlabbatch{i}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
        matlabbatch{i}.spm.util.voi.roi{2}.mask.image = {fullfile(D, S_effect(k).name,'Aud/loc/stats_s5_without_resting/mask.nii')};
        matlabbatch{i}.spm.util.voi.roi{2}.mask.threshold = 0.5;
        matlabbatch{i}.spm.util.voi.expression = 'i1 & i2';
        i=i+1;
    end
end

%%
% spm_jobman('run', matlabbatch)

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'dcm_voi';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

%% for the STS resulting from the mvpa analysis:
i=1;
clear matlabbatch
for tmp_con = 1 : length(all_con)
    select_con = all_con(tmp_con);
    for k = 1 : numel(S_effect)
        %% batch for VOI, resulting from all subj, norm>scr speech
        matlabbatch{i}.spm.util.voi.spmmat = {fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5_without_resting/SPM.mat')}; % use this .mat to extract time series.
        matlabbatch{i}.spm.util.voi.adjust = F_con;
        matlabbatch{i}.spm.util.voi.session = 1;
        matlabbatch{i}.spm.util.voi.name = sprintf('lSTS_ROI_from_MVPA_sph_DCM_ROI_adapted_to_aud_con%d_adj_eoi5', select_con);
        matlabbatch{i}.spm.util.voi.roi{1}.mask.image = {'/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/mvpa_without_resting_ttest_s8/results_PW_vs_Words/all_participants_pw_vs_words_left_cluster_10-3_510-2_mvpa.nii'};
        matlabbatch{i}.spm.util.voi.roi{1}.mask.threshold = 0.5;        
        matlabbatch{i}.spm.util.voi.roi{2}.mask.image = {fullfile(D, S_effect(k).name,'Aud/loc/stats_s5_without_resting/mask.nii')};
        matlabbatch{i}.spm.util.voi.roi{2}.mask.threshold = 0.5;
        matlabbatch{i}.spm.util.voi.expression = 'i1 & i2';
        i=i+1;
    end
end

%%
% spm_jobman('run', matlabbatch)

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'dcm_voi';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

%% for the left angular gyrus from univariate comparison words > PW:
i=1;
clear matlabbatch
for tmp_con = 1 : length(all_con)
    select_con = all_con(tmp_con);
    for k = 1 : numel(S_effect)
        %% batch for VOI, resulting from all subj, norm>scr speech
        matlabbatch{i}.spm.util.voi.spmmat = {fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5_without_resting/SPM.mat')}; % use this .mat to extract time series.
        matlabbatch{i}.spm.util.voi.adjust = F_con;
        matlabbatch{i}.spm.util.voi.session = 1;
        matlabbatch{i}.spm.util.voi.name = sprintf('lAG_ROI_DCM_ROI_adapted_to_aud_con%d_adj_eoi5', select_con);
        matlabbatch{i}.spm.util.voi.roi{1}.mask.image = {'/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/mvpa_without_resting_ttest_s8/results_PW_vs_Words/all_participants_pw_vs_words_left_cluster_10-3_510-2_mvpa.nii'};
        matlabbatch{i}.spm.util.voi.roi{1}.mask.threshold = 0.5;        
        matlabbatch{i}.spm.util.voi.roi{2}.mask.image = {fullfile(D, S_effect(k).name,'Aud/loc/stats_s5_without_resting/mask.nii')};
        matlabbatch{i}.spm.util.voi.roi{2}.mask.threshold = 0.5;
        matlabbatch{i}.spm.util.voi.expression = 'i1 & i2';
        i=i+1;
    end
end

%%
% spm_jobman('run', matlabbatch)

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'dcm_voi';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);


return
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% visual rois
% coord from TTS network (speech > baseline, syn > con), 
% VFWA from w > htf in syn + con
% SMG is from w > fh in syn + con
% occip is from all_stim in syn + con;

% ROIs_names = {'common_SMG', 'common_VWFA',  'common_occip',  'lIPS',     'MFG',    'lpSTG',   'lSTGm'};
% ROIs_coord = [-52 -44 20;   -50 -54 -14;    -28 -84 -12;    -32 -48 40; -40 6 30; -62 -31 3;  -45 -24 8]; 

ROIs_names = {'common_SMG', 'common_VWFA',  'common_occip'};
ROIs_coord = [-52 -44 20;   -50 -54 -14;    -28 -84 -12]; 

F_con = 21;
select_con = 12;
clear matlabbatch
i=1;
for k = 1 : numel(S_effect)
    for tmp_ROI = 1 : numel(ROIs_names)
        %% batch for VOI, resulting from all subj, norm>scr speech
        matlabbatch{i}.spm.util.voi.spmmat = {fullfile(D, S_effect(k).name, 'Vis/loc/stats_s5_without_resting/SPM.mat')}; % use this .mat to extract time series.
        matlabbatch{i}.spm.util.voi.adjust = F_con;
        matlabbatch{i}.spm.util.voi.session = 1;
        matlabbatch{i}.spm.util.voi.name = sprintf('%s_%d_%d_%d_%dmm_sph_DCM_ROI_adapted_to_vis_con%d_adj_eoi5', ROIs_names{tmp_ROI}, ROIs_coord(tmp_ROI, :), roi_sph_rad, select_con);
        matlabbatch{i}.spm.util.voi.roi{1}.spm.spmmat = {fullfile(D, S_effect(k).name, 'Vis/loc/stats_s5_without_resting/SPM.mat')}; % use this .mat to deduct the highest activated voxel.
        matlabbatch{i}.spm.util.voi.roi{1}.spm.contrast = select_con;
        matlabbatch{i}.spm.util.voi.roi{1}.spm.conjunction = 1;
        matlabbatch{i}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
        matlabbatch{i}.spm.util.voi.roi{1}.spm.thresh = 1;
        matlabbatch{i}.spm.util.voi.roi{1}.spm.extent = 0;
        matlabbatch{i}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
        matlabbatch{i}.spm.util.voi.roi{2}.sphere.centre = ROIs_coord(tmp_ROI,:);
        matlabbatch{i}.spm.util.voi.roi{2}.sphere.radius = 6;
        matlabbatch{i}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
        matlabbatch{i}.spm.util.voi.roi{3}.sphere.centre = ROIs_coord(tmp_ROI,:);
        matlabbatch{i}.spm.util.voi.roi{3}.sphere.radius = roi_sph_rad;
        matlabbatch{i}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
        matlabbatch{i}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';
        matlabbatch{i}.spm.util.voi.expression = 'i1&i3';
        i=i+1;
    end
end

% spm_jobman('run', matlabbatch)

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'dcm_voi';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);
