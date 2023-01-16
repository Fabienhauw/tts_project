% This script is to calculate the % of overlapping between each subject speech activation and reading activation.

clear;clc;

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

%% parameters
p_values        = [0.05 0.01 0.005 0.001 0.0005 0.0001];
p_values_text   = {'510-2' '10-2' '510-3' '10-3' '510-4' '10-4'};
p_values_field  = {'05' '01' '005' '001' '0005' '0001'};
per_vox = 10; % percentage of voxels you want to keep
use_best_vox = 0; % decide if you want to keep only X% best voxels or just use a general threshold;

%% subjects selection
Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17 

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
% gaucher_appar = {'Control02|Control04|Control05|Control07|Control17|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con_app = S;
S_con_app(mask_gauch_con) = [];
con_group = repmat({'control'}, 1, length(S_con_app));

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_droit = S;
S_droit(mask_gauch) = [];
syn_group = repmat({'syn'}, 1, length(S_droit));

categ = [syn_group, con_group];
nsub(1) = length(S_droit); % synesthetes
nsub(2) = length(S_con_app); % controls
ngroups = 2;

S_effect = [S_droit ; S_con_app];

percent_intersection = '';

%% first is to catch all speech activation and reading activation maps.

for i_subj = 1 : length(S_effect)
    speech_spmfiles_select{i_subj}  = fullfile(D, S_effect(i_subj).name, 'Aud/loc/stats_s5_without_resting/SPM.mat');
    speech_nifti{i_subj}            = fullfile(D, S_effect(i_subj).name, 'Aud/loc/stats_s5_without_resting/spmT_0011.nii'); % con 16 = all speech > rest ; con 11 = phonology
    read_spmfiles_select{i_subj}    = fullfile(D, S_effect(i_subj).name, 'Vis/loc/stats_s5_without_resting/SPM.mat');
    read_nifti{i_subj}              = fullfile(D, S_effect(i_subj).name, 'Vis/loc/stats_s5_without_resting/spmT_0012.nii'); % con 12
    mask_nifti{i_subj}              = fullfile(D, S_effect(i_subj).name, 'Aud/loc/stats_s5_without_resting/mask.nii');
end

%% second we loop over subjects, with multiple p-value threshold.
% this part is just if you want to intersect with another mask, eg. language mask
select_mask_file    = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/masks/left_language_mask_from_aal3.nii';
select_mask_header  = spm_vol(select_mask_file);
select_mask_vol     = spm_read_vols(select_mask_header);
        

for i_subj = 1 : length(S_effect)
    for tmp_p = 1 : length(p_values)
        % p_val values:
        p_val       = p_values(tmp_p); % good numeric format
        p_val_text  = p_values_text(tmp_p); % good format for name (txt)
        p_val_field = p_values_field(tmp_p); % good format for struct field...
        
        % read aud nifti
        select_speech_file    = speech_nifti{i_subj};
        select_speech_header  = spm_vol(select_speech_file);
        select_speech_vol     = spm_read_vols(select_speech_header);
        
        % read vis nifti
        select_read_file    = read_nifti{i_subj};
        select_read_header  = spm_vol(select_read_file);
        select_read_vol     = spm_read_vols(select_read_header);
        
        % read mask nifti
        select_mask_file    = mask_nifti{i_subj};
        select_mask_header  = spm_vol(select_mask_file);
        select_mask_vol     = spm_read_vols(select_mask_header);
        
        
        %% speech network
        load(speech_spmfiles_select{i_subj}); %%% load  SPM.mat
        criticalt = tinv(1-p_val,SPM.xX.erdf);  % t value threshold for the corresponding p value
        % select speech network voxels below p-val threshold
        if ~use_best_vox
            select_speech_vol(select_speech_vol<criticalt) = 0;
            select_speech_vol(select_speech_vol>=criticalt) = 1;
        elseif use_best_vox
            % OR select speech network X% voxels
            % first intersection speech X brain mask for each subject
            for x=1:size(select_speech_vol,1)
                for y=1:size(select_speech_vol,2)
                    for z=1:size(select_speech_vol,3)
                        if select_mask_vol(x,y,z)==1 % to exclude voxels which are outside the statistical map mask.nii
                            searchvol(x,y,z)=1;
                        else
                            searchvol(x,y,z)=0;
                        end
                    end
                end
            end
            
            select_speech_vol(searchvol(:)==0) = Inf; % every voxel outside of the mask is set to 'Inf'
            tvol = select_speech_vol(:) .* (searchvol(:)>0); % keep only the voxels of the speech network (select_speech_vol) which are within the brain mask (searchvol)
            
            [tvalues,xyz] = sort(tvol,'descend');
            xyz(isnan(tvalues)) = [];
            tvalues(isnan(tvalues)) = [];
            nvox = round(length(find(tvol > 0))*(per_vox/100)); % this keeps X% of the spmT voxels, within the mask, not the 10% voxel of the mask;
            xyz = xyz(1:nvox);
            select_speech_vol_copy = select_speech_vol;
            select_speech_vol_copy(:) = 0;
            select_speech_vol_copy(xyz) = 1;
            select_speech_header_copy = select_speech_header;
            select_speech_header_copy.fname = sprintf('%s_%d_percent_best_vox.nii',select_speech_header_copy.fname(1:end-4), per_vox);
            
            spm_write_vol(select_speech_header_copy,select_speech_vol_copy);
            select_speech_vol = select_speech_vol_copy;
        end

        %% reading network 
        load(read_spmfiles_select{i_subj}); %%% load  SPM.mat
        criticalt = tinv(1-p_val,SPM.xX.erdf);  % t value threshold for the corresponding p value
        % select reading network above p-val threshold
        select_read_vol(select_read_vol<criticalt) = 0;
        select_read_vol(select_read_vol>=criticalt) = 1;
        
        
        % select intersection from speech and reading networks
        for x=1:size(select_read_vol,1)
            for y=1:size(select_read_vol,2)
                for z=1:size(select_read_vol,3)
                    if select_read_vol(x,y,z)==1 & select_speech_vol(x,y,z)==1 % to exclude voxels which are outside the statistical map mask.nii
                        intersection(x,y,z)=1;
                    else
                        intersection(x,y,z)=0;
                    end
                end
            end
        end
        
%         % select voxels from reading network which are in the language mask
%         for x=1:size(select_read_vol,1)
%             for y=1:size(select_read_vol,2)
%                 for z=1:size(select_read_vol,3)
%                     if select_mask_vol(x,y,z)==1 & select_read_vol(x,y,z)==1 % to exclude voxels which are outside the statistical map mask.nii
%                         read_mask_vol(x,y,z)=1;
%                     else
%                         read_mask_vol(x,y,z)=0;
%                     end
%                 end
%             end
%         end
        
%         % select voxels from speech network which are in the language mask
%         for x=1:size(select_speech_vol,1)
%             for y=1:size(select_speech_vol,2)
%                 for z=1:size(select_speech_vol,3)
%                     if select_mask_vol(x,y,z)==1 & select_speech_vol(x,y,z)==1 % to exclude voxels which are outside the statistical map mask.nii
%                         speech_mask_vol(x,y,z)=1;
%                     else
%                         speech_mask_vol(x,y,z)=0;
%                     end
%                 end
%             end
%         end
        
%         % select intersection from speech and reading networks within the language mask
%         for x=1:size(read_mask_vol,1)
%             for y=1:size(read_mask_vol,2)
%                 for z=1:size(read_mask_vol,3)
%                     if read_mask_vol(x,y,z)==1 & speech_mask_vol(x,y,z)==1 % to exclude voxels which are outside the statistical map mask.nii
%                         intersection_mask(x,y,z)=1;
%                     else
%                         intersection_mask(x,y,z)=0;
%                     end
%                 end
%             end
%         end
        
        cd(fullfile(D, S_effect(i_subj).name, 'Vis/loc/stats_s5_without_resting'))
        select_read_header_copy = select_read_header;
        select_read_header_copy.fname = sprintf('intersection_map_reading_speech_networks_%s_%s.nii',S_effect(i_subj).name, p_val_text{1});
        spm_write_vol(select_read_header_copy,intersection);
        
%         select_read_header_copy = select_read_header;
%         select_read_header_copy.fname = sprintf('intersection_map_reading_speech_networks_language_mask_%s_%s.nii',S_effect(i_subj).name, p_val_text{1});
%         spm_write_vol(select_read_header_copy,intersection_mask);

        num_vox_intersec = length(find(intersection==1));
        num_vox_reading = length(find(select_read_vol ==1));
        num_vox_speech = length(find(select_speech_vol ==1));
%         num_vox_intersec_mask = length(find(intersection_mask==1));
%         num_vox_reading_mask = length(find(read_mask_vol ==1));
%         num_vox_speech_mask = length(find(speech_mask_vol ==1));

        percent_intersection(i_subj).(sprintf('num_vox_commun_%s', p_val_field{1}))             = num_vox_intersec;
        percent_intersection(i_subj).(sprintf('num_vox_reading_network_%s', p_val_field{1}))    = num_vox_reading;
        percent_intersection(i_subj).(sprintf('num_vox_speech_network_%s', p_val_field{1}))     = num_vox_speech;
        
%         percent_intersection(i_subj).(sprintf('num_vox_commun_mask_%s', p_val_field{1}))             = num_vox_intersec_mask;
%         percent_intersection(i_subj).(sprintf('num_vox_reading_network_mask_%s', p_val_field{1}))    = num_vox_reading_mask;
%         percent_intersection(i_subj).(sprintf('num_vox_speech_network_mask_%s', p_val_field{1}))     = num_vox_speech_mask;

        percent_intersection(i_subj).(sprintf('percentage_intersection_reading_net_%s', p_val_field{1}))    = num_vox_intersec/num_vox_reading*100;
        percent_intersection(i_subj).(sprintf('percentage_intersection_speech_net_%s', p_val_field{1}))    = num_vox_intersec/num_vox_speech*100;
        percent_intersection(i_subj).(sprintf('jaccard_index_%s', p_val_field{1}))              = num_vox_intersec/(num_vox_reading ...
            + num_vox_speech - num_vox_intersec) *100;
        
%         percent_intersection(i_subj).(sprintf('percentage_intersection_reading_net_mask_%s', p_val_field{1}))    = num_vox_intersec_mask/num_vox_reading_mask*100;
%         percent_intersection(i_subj).(sprintf('percentage_intersection_speech_net_mask_%s', p_val_field{1}))    = num_vox_intersec_mask/num_vox_speech_mask*100;
%         percent_intersection(i_subj).(sprintf('jaccard_index_mask_%s', p_val_field{1}))              = num_vox_intersec_mask/(num_vox_reading_mask ...
%             + num_vox_speech_mask - num_vox_intersec_mask) *100;
        
    end
end
if use_best_vox
    best_vox_percent_inter = percent_intersection;
    clear percent_intersection;
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

