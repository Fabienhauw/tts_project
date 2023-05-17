% This script is to calculate the % of overlapping between each subject speech activation and reading activation.

clear;clc;

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

%% parameters
p_values        = [0.01 0.005 0.001 0.0001];
p_values_text   = {'10-2' '510-3' '10-3' '10-4'};
p_values_field  = {'01' '005' '001' '0001'};
aud_con = 16;

%% subjects selection
Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17 

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
% gaucher_appar = {'Control02|Control04|Control05|Control07|Control17|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con = S;
S_con(mask_gauch_con) = [];
con_group = repmat({'control'}, 1, length(S_con));

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_syn = S;
S_syn(mask_gauch) = [];
syn_group = repmat({'syn'}, 1, length(S_syn));

categ = [syn_group, con_group];
nsub(1) = length(S_syn); % synesthetes
nsub(2) = length(S_con); % controls
ngroups = 2;

S_effect = [S_syn ; S_con];

percent_intersection = '';

%% first is to catch all speech activation and reading activation maps.

for i_subj = 1 : length(S_effect)
    speech_spmfiles_select{i_subj}  = fullfile(D, S_effect(i_subj).name, 'Aud/loc/stats_s5_without_resting/SPM.mat');
    speech_nifti{i_subj}            = fullfile(D, S_effect(i_subj).name, sprintf('Aud/loc/stats_s5_without_resting/spmT_00%d.nii', aud_con)); % con 16 = all speech > rest ; con 11 = phonology
    mask_nifti{i_subj}              = fullfile(D, S_effect(i_subj).name, 'Aud/loc/stats_s5_without_resting/mask.nii');
end

%% second we loop over subjects, with multiple p-value threshold.

for tmp_p = length(p_values)
    for i_subj = 1 : length(S_effect)
        % p_val values:
        p_val       = p_values(tmp_p); % good numeric format
        p_val_text  = p_values_text(tmp_p); % good format for name (txt)
        p_val_field = p_values_field(tmp_p); % good format for struct field...
        
        % read aud nifti
        select_speech_file    = speech_nifti{i_subj};
        select_speech_header  = spm_vol(select_speech_file);
        select_speech_vol     = spm_read_vols(select_speech_header);
        
        % read mask nifti
        select_mask_file    = mask_nifti{i_subj};
        select_mask_header  = spm_vol(select_mask_file);
        select_mask_vol     = spm_read_vols(select_mask_header);
        
        
        %% speech network
        load(speech_spmfiles_select{i_subj}); %%% load  SPM.mat
        criticalt = tinv(1-p_val,SPM.xX.erdf);  % t value threshold for the corresponding p value
        % select speech network voxels below p-val threshold
        select_speech_vol(select_speech_vol<criticalt) = 0;
        select_speech_vol(select_speech_vol>=criticalt) = 1;
        
        select_speech_header_copy = select_speech_header;
        select_speech_header_copy.fname = sprintf('%s_thresh_%s.nii', select_speech_header.fname(1:end-4), p_val_text{1});
        spm_write_vol(select_speech_header_copy,select_speech_vol);
        
        final_images(i_subj).(sprintf('thresh_%s', p_val_field{1})) = select_speech_header_copy.fname;
    end
    
    outdir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/sum_thresholded_maps_con';
    if ~isdir(outdir)
        mkdir(outdir)
    end
    scans1 = {};
    scans2 = {};
    scans3 = {};
    for tmp_img = 1 : numel(final_images)/2
        scans1 = [scans1; final_images(tmp_img).(sprintf('thresh_%s', p_val_field{1}))];
    end
    for tmp_img = numel(final_images)/2 + 1 : numel(final_images)
        scans2 = [scans2; final_images(tmp_img).(sprintf('thresh_%s', p_val_field{1}))];
    end
    
    total1 = '';
    for k = 1:length(scans1)
        map = sprintf('i%d',k);
        if k == 1
            total1 = map;
        else
            total1 = [total1 '+' map];
        end
    end
    total1 = ['(' total1 ')'];
    
    total2 = '';
    for k = 1:length(scans2)
        map = sprintf('i%d',k);
        if k == 1
            total2 = map;
        else
            total2 = [total2 '+' map];
        end
    end
    total2 = ['(' total2 ')'];
    
    clear matlabbatch
    
    matlabbatch{1}.spm.util.imcalc.input = scans1;
    matlabbatch{1}.spm.util.imcalc.output = sprintf('sum_speech>baseline_syn_%s.nii', p_val_field{1});
    matlabbatch{1}.spm.util.imcalc.outdir = {outdir};
    matlabbatch{1}.spm.util.imcalc.expression = [total1];
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = -7;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    matlabbatch{2}.spm.util.imcalc.input = scans2;
    matlabbatch{2}.spm.util.imcalc.output = sprintf('sum_speech>baseline_con_%s.nii', p_val_field{1});
    matlabbatch{2}.spm.util.imcalc.outdir = {outdir};
    matlabbatch{2}.spm.util.imcalc.expression = [total2];
    matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{2}.spm.util.imcalc.options.mask = 0;
    matlabbatch{2}.spm.util.imcalc.options.interp = -7;
    matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
    
    spm_jobman('run', matlabbatch)
    
    clear matlabbatch
    scans3 = {
        fullfile(outdir, sprintf('sum_speech>baseline_syn_%s.nii', p_val_field{1}))
        fullfile(outdir, sprintf('sum_speech>baseline_con_%s.nii', p_val_field{1}))
        };
    matlabbatch{1}.spm.util.imcalc.input = scans3;
    matlabbatch{1}.spm.util.imcalc.output = sprintf('diff_speech>baseline_syn>con_%s', p_val_field{1});
    matlabbatch{1}.spm.util.imcalc.outdir = {outdir};
    matlabbatch{1}.spm.util.imcalc.expression = ['(i1 - i2)'];
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = -7;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    matlabbatch{2}.spm.util.imcalc.input = scans3;
    matlabbatch{2}.spm.util.imcalc.output = sprintf('diff_speech>baseline_con>syn_%s', p_val_field{1});
    matlabbatch{2}.spm.util.imcalc.outdir = {outdir};
    matlabbatch{2}.spm.util.imcalc.expression = ['(i2 - i1)'];
    matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{2}.spm.util.imcalc.options.mask = 0;
    matlabbatch{2}.spm.util.imcalc.options.interp = -7;
    matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
    
    spm_jobman('run', matlabbatch)
    
end
