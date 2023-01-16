%% creates VOI for DCM for controls, based on the locations of controls normal>reverse activations

clear;clc;
i=1;
Batches.path_to_subject;
mask = ismember({S.name}, {'.', '..', 'FLC_func'});
S(mask) = [];
subjs={};
nb_subj=numel(S);
path_to_scan = {};
repet_time = {};
for k=1:numel(S)
    subjs=[subjs;S(k).name];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
    repet_time = [repet_time;1.3];
end

Batches.rosso.path_to_subject;
nb_subj=nb_subj + numel(S);
for k=1:numel(S)
    subjs=[subjs;S(k).name];
    repet_time = [repet_time;3];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
end

for k = 1:nb_subj
    cd(fullfile(path_to_scan{k},'anat'));
    mask = dir('brain*.nii');
    mask = fullfile(path_to_scan{k},'anat',mask(1).name);
    
    if exist(fullfile(path_to_scan{k}, 'Aud/stats/chaperon/stc'))
        path_to_stats1=fullfile(path_to_scan{k}, 'Aud/stats/chaperon/stc');
        path_to_stats2=fullfile(path_to_scan{k}, 'Aud/stats/chaperon/stc/dcm_model_param_modul_endroit');
    else
        path_to_stats1=fullfile(path_to_scan{k}, 'Aud/stats/chaperon');
        path_to_stats2=fullfile(path_to_scan{k}, 'Aud/stats/chaperon/dcm_model_param_modul_endroit');
    end
        
    cd (path_to_stats2)
    stat_mask = dir('mask.nii');
    stat_mask = fullfile(path_to_stats2,stat_mask(1).name);
    
    path_to_spm = fullfile(path_to_stats1, 'SPM.mat');
    path_to_spm2 = fullfile(path_to_stats2, 'SPM.mat');
    load(path_to_spm2);
    
    cd(path_to_stats2)
    for name = 1 : length(SPM.xY.VY)
        SPM.xY.VY(name).fname = strrep(SPM.xY.VY(name).fname,'/media/biggest_drive','/media/fabien.hauw/biggest_drive');
    end
    SPM.swd = strrep(SPM.swd,'/media/biggest_drive','/media/fabien.hauw/biggest_drive');
    save('SPM.mat', 'SPM');
    
    for con=1:length(SPM.xCon)
        if SPM.xCon(con).STAT=='F'
            f_con=con;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%   VOI creations    %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %% batch for VOI stgm: localisation of the VOI from the reverse>rest in controls only
    matlabbatch{i}.spm.util.voi.spmmat = {path_to_spm2}; % the first spm, to extract the time series from;
    matlabbatch{i}.spm.util.voi.adjust = f_con;
    matlabbatch{i}.spm.util.voi.session = 1;
    matlabbatch{i}.spm.util.voi.name = 'dcm_STGm_perso';
    matlabbatch{i}.spm.util.voi.roi{1}.spm.spmmat = {path_to_spm}; % the second spm, used to define the peak of activations in a certain contrast;
    matlabbatch{i}.spm.util.voi.roi{1}.spm.contrast = 2; %this is the contrast reverse-rest in the "standard" SPM in the chaperon run
    matlabbatch{i}.spm.util.voi.roi{1}.spm.conjunction = 1;
    matlabbatch{i}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
    matlabbatch{i}.spm.util.voi.roi{1}.spm.thresh = 1;
    matlabbatch{i}.spm.util.voi.roi{1}.spm.extent = 0;
    matlabbatch{i}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
    matlabbatch{i}.spm.util.voi.roi{2}.sphere.centre = [-48 -16 5];
    matlabbatch{i}.spm.util.voi.roi{2}.sphere.radius = 6;
    matlabbatch{i}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
%     matlabbatch{i}.spm.util.voi.roi{3}.sphere.centre = [0 0 0];
%     matlabbatch{i}.spm.util.voi.roi{3}.sphere.radius = 6;
%     matlabbatch{i}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
%     matlabbatch{i}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';
    matlabbatch{i}.spm.util.voi.roi{3}.mask.image = {mask};
    matlabbatch{i}.spm.util.voi.roi{3}.mask.threshold = 0.5;
    matlabbatch{i}.spm.util.voi.expression = 'i1&i3&i2';
    i=i+1;
    
    %% batch for VOI smg: localisation of the VOI from the normal>reverse in flc>controls
    matlabbatch{i}.spm.util.voi.spmmat = {path_to_spm2};
    matlabbatch{i}.spm.util.voi.adjust = f_con;
    matlabbatch{i}.spm.util.voi.session = 1;
    matlabbatch{i}.spm.util.voi.name = 'dcm_SMG'; %comes from the visual localizer because the vwfa from FLC-controls is too small;
%     matlabbatch{i}.spm.util.voi.roi{1}.mask.image = {'/media/biggest_drive/Fabien/second_level_rosso/Aud/chaperon_stc_minus_bad_temoins/FLC_vs_all_minus_bad_temoins_age_covariate/normal-reverse/smg_flc-others_fwd-bwd_10-3_510-2.nii'};
    matlabbatch{i}.spm.util.voi.roi{1}.mask.image = {'/media/fabien.hauw/biggest_drive/Fabien/second_level_rosso/Aud/chaperon_stc_minus_sujet08/FLC_vs_all_age_covariate/normal-reverse/smg_flc-others_fwd-bwd_10-3_510-2.nii'};
    matlabbatch{i}.spm.util.voi.roi{1}.mask.threshold = 0.5;
    matlabbatch{i}.spm.util.voi.expression = 'i1';
    i=i+1;
    
    %% batch for VOI VWFA
    matlabbatch{i}.spm.util.voi.spmmat = {path_to_spm2};
    matlabbatch{i}.spm.util.voi.adjust = f_con;
    matlabbatch{i}.spm.util.voi.session = 1;
    matlabbatch{i}.spm.util.voi.name = 'dcm_vwfa_sphere';
    matlabbatch{i}.spm.util.voi.roi{1}.mask.image = {'/media/fabien.hauw/biggest_drive/manip_kasia/first_level/5001/Vis/stats/stc/vwfa_vis_loc_sphere_6.nii'};
    matlabbatch{i}.spm.util.voi.roi{1}.mask.threshold = 0.5;
    matlabbatch{i}.spm.util.voi.expression = 'i1';
    i=i+1;
    
    %% batch for VOI mtg: intersection of the MTG-vwfa cluster decoding words vs pw in flc only with a box with wide ranges to exclude vwfa (x<-51, y>-74, z>-13);
    matlabbatch{i}.spm.util.voi.spmmat = {path_to_spm2}; % the first spm, to extract the time series from;
    matlabbatch{i}.spm.util.voi.adjust = f_con;
    matlabbatch{i}.spm.util.voi.session = 1;
    matlabbatch{i}.spm.util.voi.name = 'dcm_mtg_mvpa';
    matlabbatch{i}.spm.util.voi.roi{1}.mask.image = {'/media/fabien.hauw/biggest_drive/Fabien/FinalMRI/FLC/Aud/mvpa/run/balanced_data/lexicality/results_Words_vs_PW_mask_words_pw_flc/perm/cluster/cluster_MTG_intersection_box.nii,1'};
    matlabbatch{i}.spm.util.voi.roi{1}.mask.threshold = 0.5;
    matlabbatch{i}.spm.util.voi.roi{2}.mask.image = {stat_mask};
    matlabbatch{i}.spm.util.voi.roi{2}.mask.threshold = 0.5;
    matlabbatch{i}.spm.util.voi.expression = 'i1&i2';
    i=i+1;
end

spm_jobman('run', matlabbatch)