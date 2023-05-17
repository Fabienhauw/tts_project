clear;clc;
addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))

wd = pwd;

% res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/ttest_syn_multi_vs_one_s5_without_resting';
% res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/ttest_syn_proj_vs_asso_s5_without_resting';
res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/multi_reg_s5_pca';
% res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/syn_vs_con_s5';

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
S = [Syn];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17 
% syn_pcs = 'Sujet03|Sujet04|Sujet06|Sujet08|Sujet09|Sujet10|Sujet15|Sujet17|Sujet19|Sujet20|Sujet22';

% syn_one = {'Sujet02|Sujet03|Sujet04|Sujet09|Sujet20'};
% mask_one = cellfun(@isempty,(regexp({S.name},syn_one)));
mask_syn =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% S_one(mask_one) = [];
S(mask_syn) = [];

vector_pca = [
    1.46695396813551;-13.5431937157681;1.41741598407124;26.5667386823073;-3.58433618681199;1.42514228090863;6.54959587962896;-13.4982323294988;-8.51452331666684;21.4660769814270;1.36396126465718;-8.59389647259469;21.4747879876487;-13.5367001585281;-13.4354695343489;-8.53931091167326;1.51498959710619;...
    ];

S_effect = [S];

for j = 1 : size(S,1)
    if ~isempty(find(~cellfun(@isempty,(regexp({S_effect.name},S(j).name)))))
        mask_cov(j,1) = 1;
    else
        mask_cov(j,1) = 0;
    end
end

% vector_cov1 = vector_age(mask_cov==1);
vector_cov1 = vector_pca;

clear matlabbatch
i=1;

%% to erase previous models
erase_model = input('Do you want to erase previous 2nd lvl auditive models? [yes/no] ', 's');
if isequal(erase_model,'yes')
    redo_model = 1;
elseif isequal(erase_model,'no')
    redo_model = 0;
end

%% second level:
clear matlabbatch
i=1;

cd(fullfile(D, S(1).name,'Aud/loc/stats_s5_without_resting'))
all_contrast = dir(sprintf('s8con*.nii'));
all_contrast = all_contrast(~cellfun(@isempty,(regexp({all_contrast.name}, '^s8con_\d+\.nii'))));

names = {...
    'words', 'pseudowords', 'numbers', 'normal_speech', 'scramble_speech',...
    'odds', 'motor',...
    'lexicality', '-lexicality','(words + normal_speech + numbers) - pseudowords',...
    'phonology', '-phonology', 'words + pseudowords', 'numbers - (words + pseudowords)', ...
    'normal_speech - words', 'words + pseudowords + numbers + normal_speech', ...
    'all', 'words - normal_speech', 'words - scramble_speech', 'pseudowords - scramble_speech', '(words+pseudowords) -  scramble_speech', ...
    '(words + normal_speech) -  scramble_speech', '(normal_speech + pseudowords) -  scramble_speech', '(normal_speech + pseudowords + words) -  scramble_speech',...
    '(normal_speech + pseudowords + words + numbers) - scramble_speech',...
    'numbers - words', 'words - numbers', ...
    'EOI', 'EOI2', ...
    };

load ('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/comportemental_exp/res_syn_exp123.mat')

for c = 16
%     for tmp_cov = 1 : size(sum_res_ter,2)
        %% subject of interest:
        contrast = all_contrast(c).name;
        scans1 = {};
        for k = 1:length(S_effect)
            vol_name = fullfile(D, S_effect(k).name,'Aud/loc/stats_s5_without_resting/');
            vol_name = sprintf('%s%s,1',vol_name,contrast);
            scans1 = [scans1;vol_name];
        end
        res_dir = fullfile(res_dir_base,names{c});
%         res_dir = fullfile(res_dir_base,names{c}, sprintf('cov%d', tmp_cov));
        if ~isdir(res_dir)
            mkdir(res_dir)
        end
        dinfo = dir(res_dir);
        dinfo([dinfo.isdir]) = [];   %skip directories
        filenames = fullfile(res_dir, {dinfo.name});
        if ~exist('redo_model','var')
            redo_model=0;
        end
        if redo_model & ~isempty (filenames)
            delete (filenames{:});
        end
        
%         vector_cov1 = sum_res_ter(:, tmp_cov);
    
        matlabbatch{i}.spm.stats.factorial_design.dir = {res_dir};
        matlabbatch{i}.spm.stats.factorial_design.des.mreg.scans = scans1;
        
        matlabbatch{i}.spm.stats.factorial_design.des.mreg.mcov.c = vector_cov1;
        matlabbatch{i}.spm.stats.factorial_design.des.mreg.mcov.cname = 'cov';
        matlabbatch{i}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;
        matlabbatch{i}.spm.stats.factorial_design.des.mreg.incint = 1;
        
        matlabbatch{i}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{i}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{i}.spm.stats.factorial_design.masking.im = 1;
        exp_mask = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks/Sujet01_Control21_aud_loc_mask_thr_s5.nii');
        matlabbatch{i}.spm.stats.factorial_design.masking.em = {exp_mask};
        matlabbatch{i}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{i}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{i}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        i = i+1;
        
        matlabbatch{i}.spm.stats.fmri_est.spmmat = {fullfile(res_dir,'SPM.mat')};
        matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
        
        i=i+1;
        
        matlabbatch{i}.spm.stats.con.spmmat = {fullfile(res_dir, 'SPM.mat')};
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'cov';
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = [0 1];
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        
        matlabbatch{i}.spm.stats.con.consess{2}.tcon.name = '-cov';
        matlabbatch{i}.spm.stats.con.consess{2}.tcon.weights = [0 -1];
        matlabbatch{i}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        
        matlabbatch{i}.spm.stats.con.delete = 1;
        
        i = i+1;
%     end
end

spm_jobman('run', matlabbatch)

cd(res_dir)