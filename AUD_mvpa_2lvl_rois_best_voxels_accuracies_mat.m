% to compare matrices of accuracy between controls and synesthetes in ROIs defined as spheres around peak coordinates of activation by norm>scr
% speech in auditive run.
clear;
addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))

wd = pwd;

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

con_dir = sprintf('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/mvpa10mm_rois/', S_effect(1).name);
cd(con_dir)
ROIs = dir('*best_vox*');
cd(ROIs(1).name)
results = dir('results*ROI');
for tmp_comp = 1 : length(results)
    str_path = strsplit(results(tmp_comp).name, '_');
    if ~isempty(regexp(results(tmp_comp).name, 'ROI'))
        all_comp_inv{tmp_comp} = strjoin({str_path{1}, str_path{4}, str_path{3}, str_path{2}, str_path{5}}, '_');
    else
        all_comp_inv{tmp_comp} = strjoin({str_path{1}, str_path{4}, str_path{3}, str_path{2}}, '_');
    end
    all_comp{tmp_comp} = results(tmp_comp).name;
end

%% stats:

for tmp_roi = 1 : length(ROIs)
    for tmp_comp = 1 : length(all_comp)
        %% group 1:
        for k = 1:length(S_syn)
            if ~isdir(fullfile(D, S_syn(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp_inv{tmp_comp}))
                tmp_file = fullfile(D, S_syn(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp{tmp_comp}, '/res_accuracy_minus_chance.mat');
            else
                tmp_file = fullfile(D, S_syn(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp_inv{tmp_comp}, '/res_accuracy_minus_chance.mat');
            end
            load(tmp_file);
            summary(tmp_comp).comparison_name = all_comp{tmp_comp};
            summary(tmp_comp).analytics(tmp_roi).name = ROIs(tmp_roi).name;
            summary(tmp_comp).analytics(tmp_roi).accuracies1{k,1} = results.accuracy_minus_chance.output;
            summary(tmp_comp).analytics(tmp_roi).paths1{k,1} = tmp_file;
        end
        summary(tmp_comp).analytics(tmp_roi).mean_accuracies1 = mean([summary(tmp_comp).analytics(tmp_roi).accuracies1{:}]);
        
        %% for group 2:
        for k = 1:length(S_con) %S_con if all controls
            if ~isdir(fullfile(D, S_con(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp_inv{tmp_comp}))
                tmp_file = fullfile(D, S_con(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp{tmp_comp}, '/res_accuracy_minus_chance.mat');
            else
                tmp_file = fullfile(D, S_con(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp_inv{tmp_comp}, '/res_accuracy_minus_chance.nii');
            end
            load(tmp_file);
            summary(tmp_comp).analytics(tmp_roi).accuracies2{k,1} = results.accuracy_minus_chance.output;
            summary(tmp_comp).analytics(tmp_roi).paths2{k,1} = tmp_file;
        end
        summary(tmp_comp).analytics(tmp_roi).mean_accuracies2 = mean([summary(tmp_comp).analytics(tmp_roi).accuracies2{:}]);
        
        [h{tmp_roi,tmp_comp},p{tmp_roi,tmp_comp},ci{tmp_roi,tmp_comp},stats{tmp_roi,tmp_comp}] = ...
            ttest2(...
            [summary(tmp_comp).analytics(tmp_roi).accuracies1{:}], ...
            [summary(tmp_comp).analytics(tmp_roi).accuracies2{:}]);
        tval{tmp_roi,tmp_comp} = stats{tmp_roi,tmp_comp}.tstat;
        
        summary(tmp_comp).analytics(tmp_roi).p_val = p{tmp_roi,tmp_comp};
        summary(tmp_comp).analytics(tmp_roi).t_val = tval{tmp_roi,tmp_comp};
    end
end
cd(res_dir_base)
save(fullfile(res_dir_base, 'results_decoding_10_voxels_roi.mat'), 'summary');

clear summary
%%
cd(con_dir)
ROIs = dir('*10mm_sphere*');
cd(ROIs(1).name)
results = dir('results*ROI');
for tmp_comp = 1 : length(results)
    str_path = strsplit(results(tmp_comp).name, '_');
    if ~isempty(regexp(results(tmp_comp).name, 'ROI'))
        all_comp_inv{tmp_comp} = strjoin({str_path{1}, str_path{4}, str_path{3}, str_path{2}, str_path{5}}, '_');
    else
        all_comp_inv{tmp_comp} = strjoin({str_path{1}, str_path{4}, str_path{3}, str_path{2}}, '_');
    end
    all_comp{tmp_comp} = results(tmp_comp).name;
end

%% stats:


for tmp_comp = 1 : length(all_comp)
    for tmp_roi = 1 : length(ROIs)
        %% group 1:
        for k = 1:length(S_syn)
            if ~isdir(fullfile(D, S_syn(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp_inv{tmp_comp}))
                tmp_file = fullfile(D, S_syn(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp{tmp_comp}, '/res_accuracy_minus_chance.mat');
            else
                tmp_file = fullfile(D, S_syn(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp_inv{tmp_comp}, '/res_accuracy_minus_chance.mat');
            end
            load(tmp_file);
            summary(tmp_comp).comparison_name = all_comp{tmp_comp};
            summary(tmp_comp).analytics(tmp_roi).name = ROIs(tmp_roi).name;
            summary(tmp_comp).analytics(tmp_roi).accuracies1{k,1} = results.accuracy_minus_chance.output;
            summary(tmp_comp).analytics(tmp_roi).paths1{k,1} = tmp_file;
        end
        summary(tmp_comp).analytics(tmp_roi).mean_accuracies1 = mean([summary(tmp_comp).analytics(tmp_roi).accuracies1{:}]);
        
        %% for group 2:
        for k = 1:length(S_con) %S_con if all controls
            if ~isdir(fullfile(D, S_con(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp_inv{tmp_comp}))
                tmp_file = fullfile(D, S_con(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp{tmp_comp}, '/res_accuracy_minus_chance.mat');
            else
                tmp_file = fullfile(D, S_con(k).name,'Aud/loc/mvpa10mm_rois/', ROIs(tmp_roi).name, all_comp_inv{tmp_comp}, '/res_accuracy_minus_chance.mat');
            end
            load(tmp_file)
            summary(tmp_comp).analytics(tmp_roi).accuracies2{k,1} = results.accuracy_minus_chance.output;
            summary(tmp_comp).analytics(tmp_roi).paths2{k,1} = tmp_file;
        end
        summary(tmp_comp).analytics(tmp_roi).mean_accuracies2 = mean([summary(tmp_comp).analytics(tmp_roi).accuracies2{:}]);
        
        [h{tmp_roi,tmp_comp},p{tmp_roi,tmp_comp},ci{tmp_roi,tmp_comp},stats{tmp_roi,tmp_comp}] = ...
            ttest2(...
            [summary(tmp_comp).analytics(tmp_roi).accuracies1{:}], ...
            [summary(tmp_comp).analytics(tmp_roi).accuracies2{:}]);
        tval{tmp_roi,tmp_comp} = stats{tmp_roi,tmp_comp}.tstat;
        summary(tmp_comp).analytics(tmp_roi).p_val = p{tmp_roi,tmp_comp};
        summary(tmp_comp).analytics(tmp_roi).t_val = tval{tmp_roi,tmp_comp};
        
    end
end

cd(res_dir_base)
save(fullfile(res_dir_base, 'results_decoding_whole_roi.mat'), 'summary');
