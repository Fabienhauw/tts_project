% This script is to calculate the % of overlapping between each subject speech activation and reading activation.

clear;clc;

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

%% parameters
tts_net = fullfile(D, 'Aud/loc/ANOVA_s5_without_resting_Control01_to_Sujet22_s8/bin_syn>con_speech>rest_10-3_510-2_mask_syn_10-3.nii');
mask_nifti = fullfile(D, 'Aud/loc/ANOVA_s5_without_resting_Control01_to_Sujet22_s8/mask.nii');
dys_roi_dir = fullfile(D, 'dyslexia_feng_fh');
cd(dys_roi_dir)
dys_rois = dir ('8mm*.nii');
percent_intersection = '';

%% second we loop over subjects, with multiple p-value threshold.
% this part is just if you want to intersect with another mask, eg. language mask
select_tts_file    = tts_net;
select_tts_header  = spm_vol(select_tts_file);
select_tts_vol     = spm_read_vols(select_tts_header);
select_tts_vol     = round(select_tts_vol);

% read mask nifti
select_mask_file        = mask_nifti;
select_mask_header      = spm_vol(select_mask_file);
select_mask_vol         = spm_read_vols(select_mask_header);
    
for i_roi = 1 : length(dys_rois)
    % read aud nifti
    select_roi_file         = fullfile(dys_roi_dir, dys_rois(i_roi).name);
    select_roi_header       = spm_vol(select_roi_file);
    select_roi_vol          = spm_read_vols(select_roi_header);
    
    %% speech network
    % select voxels from reading network which are in the language mask
    for x=1:size(select_roi_vol,1)
        for y=1:size(select_roi_vol,2)
            for z=1:size(select_roi_vol,3)
                if select_mask_vol(x,y,z)==1 & select_roi_vol(x,y,z)==1 % to exclude voxels which are outside the statistical map mask.nii
                    roi_mask_vol(x,y,z)=1;
                else
                    roi_mask_vol(x,y,z)=0;
                end
            end
        end
    end
    %% tts network
    % select intersection from roi and tts network
    for x=1:size(select_tts_vol,1)
        for y=1:size(select_tts_vol,2)
            for z=1:size(select_tts_vol,3)
                if select_tts_vol(x,y,z)==1 & select_roi_vol(x,y,z)==1 % to exclude voxels which are outside the statistical map mask.nii
                    intersection(x,y,z)=1;
                else
                    intersection(x,y,z)=0;
                end
            end
        end
    end
    
    for x=1:size(select_tts_vol,1)
        for y=1:size(select_tts_vol,2)
            for z=1:size(select_tts_vol,3)
                if select_tts_vol(x,y,z)==1 & roi_mask_vol(x,y,z)==1 % to exclude voxels which are outside the statistical map mask.nii
                    intersection_mask(x,y,z)=1;
                else
                    intersection_mask(x,y,z)=0;
                end
            end
        end
    end
    
    cd(dys_roi_dir)
    select_tts_header_copy = select_tts_header;
    select_tts_header_copy.fname = sprintf('intersection_tts_networkmap_roi_%s.nii', dys_rois(i_roi).name(1:end-4));
    spm_write_vol(select_tts_header_copy,intersection);

    num_vox_intersec = length(find(intersection==1));
    num_vox_tts = length(find(select_tts_vol ==1));
    num_vox_roi = length(find(select_roi_vol ==1));
    
    num_vox_intersec_mask = length(find(intersection_mask==1));
    num_vox_tts_mask = length(find(select_tts_vol ==1));
    num_vox_roi_mask = length(find(roi_mask_vol ==1));

%     dys_rois(i_roi).name(1:end-4) = strrep(dys_rois(i_roi).name(1:end-4), '.', '_');

    percent_intersection(i_roi).name  = dys_rois(i_roi).name(1:end-4);
    percent_intersection(i_roi).(sprintf('num_vox_tts_network'))     = num_vox_tts;
    percent_intersection(i_roi).(sprintf('num_vox_commun_roi_tts'))  = num_vox_intersec;
    percent_intersection(i_roi).(sprintf('num_vox_roi_%s'))  = num_vox_roi;
    percent_intersection(i_roi).(sprintf('percentage_intersection_tts_roi_%s')) = num_vox_intersec/num_vox_roi*100;
    percent_intersection(i_roi).(sprintf('num_vox_commun_roi_mask_tts'))  = num_vox_intersec_mask;
    percent_intersection(i_roi).(sprintf('num_vox_roi_mask_%s'))  = num_vox_roi_mask;
    percent_intersection(i_roi).(sprintf('percentage_intersection_tts_roi_mask_%s')) = num_vox_intersec_mask/num_vox_roi_mask*100;
    
end

%% stats

if use_best_vox
    fn = fieldnames(best_vox_percent_inter);
    for k = 1 : numel(fn)
        mean_val (1, k) = mean([best_vox_percent_inter(1:end/2).(fn{k})]);
        mean_val (2, k) = mean([best_vox_percent_inter(end/2+1:end).(fn{k})]);
        syn_val = [best_vox_percent_inter(1:end/2).(fn{k})];
        con_val = [best_vox_percent_inter(end/2+1:end).(fn{k})];
        [h(k,1).fn(k).name, p(k,1).fn(k).name, t(k,1).fn(k).name, df(k,1).fn(k).name] = ttest2(syn_val, con_val);
    end
elseif ~use_best_vox
    fn = fieldnames(percent_intersection);
    for k = 1 : numel(fn)
        mean_val(1).(fn{k}) = mean([percent_intersection(1:end/2).(fn{k})]);
        mean_val(2).(fn{k}) = mean([percent_intersection(end/2+1:end).(fn{k})]);
        syn_val = [percent_intersection(1:end/2).(fn{k})];
        con_val = [percent_intersection(end/2+1:end).(fn{k})];
        [h.(fn{k}), p.(fn{k}), t.(fn{k}), df.(fn{k})] = ttest2(syn_val, con_val);
    end
end

res_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/overlap_comp';
if ~isdir(res_dir)
    mkdir(res_dir)
end
save(fullfile(res_dir, sprintf('result_comparison_overlap_aud_con%d_reading.mat',aud_con)),'percent_intersection','h','p','t', 'df');
