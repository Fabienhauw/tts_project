% to compare FLC chaperon rouge run to CRosso controls and M2 in the same
% time.
clear;
addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))

wd = pwd;

res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/one_proj_vs_all_assoc_s8';

%% to erase previous smoothed contrasts;
erase_model = input('Do you want to erase previous auditive models? [yes/no] ', 's');
if isequal(erase_model,'yes')
    redo_model = 1;
elseif isequal(erase_model,'no')
    redo_model = 0;
end

%%

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

mask_proj =  cellfun(@isempty,(regexp({S.name},'Sujet01|Sujet03|Sujet06|Sujet15')));
mask_assoc = ~cellfun(@isempty,(regexp({S.name},'Sujet01|Sujet03|Sujet06|Sujet15|Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));

S_proj = S;
S_proj(mask_proj) = [];

S_assoc = S;
S_assoc(mask_assoc) = [];

S_effect = [S_proj ; S_assoc];
subname = {S_effect.name};

totsub=length(subname);

path_to_scan_proj = {};
path_to_scan_assoc = {};
for k = 1 : numel(S_proj)
    path_to_subj = fullfile(D, S_proj(k).name, 'Aud/loc/stats_s5_without_resting');
    path_to_scan_proj = [path_to_scan_proj;path_to_subj];
end

for k = 1 : numel(S_assoc)
    path_to_subj = fullfile(D, S_assoc(k).name, 'Aud/loc/stats_s5_without_resting');
    path_to_scan_assoc = [path_to_scan_assoc;path_to_subj];
end

vector_age = [
    25.1013699; 70.8219178; 23.6821918; 24.3342466; 21.109589; 31.7260274; 18.8438356; ...
    36.3424658; 49.7753425; 27.2767123; 26.3616438; 40.8876712; 22.8246575; 43.1013699; 44.309589; ...
    51.2246575; 40.7561644; 18.5835616; 39.939726; 58.9534247; 43.1753425; 39.9260274; ... % end of synesthetes
    35.8219178; 23.2054795; 21.9726027; 31.7945205; 30.7589041; 70.3808219; 26.3589041; ... %start of controls
    22.660274; 42.0027397; 45.5589041; 19.6027397; 55.9041096; 19.4958904; 38.3643836; ...
    49.9178082; 46.0739726; 51.8630137; 25.1945205; 23.0547945; 41.4383562; 41.5972603; ...
    29.865753; 26.57534247; 24.67945205; 29.72328767; 58.97534247; ...
    ];

vector_hand = [
    0; 0; 0; 0; 1; 0; 1; 0; 0; 0; 1; 0; 0; 1; 0; 1; 0; 0; 0; 0; 0; 0;... % end of synesthetes
    zeros(21,1); 1; 1; 1; 1; 1; ... % end of controls
    ]; %0 = right, 1 = left;

pairs = {
    'S*01', 'C*19';
    'S*02', 'C*06';
    'S*03', 'C*18';
    'S*04', 'C*08';
    'S*05', 'C*24';
    'S*06', 'C*01';
    'S*07', 'C*23';
    'S*08', 'C*09';
    'S*09', 'C*15';
    'S*10', 'C*11';
    'S*11', 'C*25';
    'S*12', 'C*21';
    'S*13', 'C*03';
    'S*14', 'C*22';
    'S*15', 'C*16';
    'S*16', 'C*26';
    'S*17', 'C*20';
    'S*18', 'C*13';
    'S*19', 'C*14';
    'S*20', 'C*12';
    'S*21', 'C*10';
    'S*22', 'C*05';
    };

for j = 1 : size(S,1)
    if ~isempty(find(~cellfun(@isempty,(regexp({S_proj.name},S(j).name)))))
        mask_cov_proj(j,1) = 1;
    else
        mask_cov_proj(j,1) = 0;
    end
end
for j = 1 : size(S,1)
    if ~isempty(find(~cellfun(@isempty,(regexp({S_assoc.name},S(j).name)))))
        mask_cov_assoc(j,1) = 1;
    else
        mask_cov_assoc(j,1) = 0;
    end
end

vector_cov_proj1 = vector_age(mask_cov_proj==1);
vector_cov_proj2 = vector_hand(mask_cov_proj==1);
vector_cov_assoc1 = vector_age(mask_cov_assoc==1);
vector_cov_assoc2 = vector_hand(mask_cov_assoc==1);

if ~exist ('i', 'var')
    i=1;
end

% get image files names
cd(path_to_scan_assoc{1})
cont = dir('con_00*.nii');
all_cont = [1 : length(cont)];
all_cont_names = {cont.name};
ncon = length(cont);
P={};
for con=1:ncon
    for s=1:totsub
        sub=subname{s};
        P{(con-1)*totsub+s} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/stats_s5_without_resting/s8%s'], sub, all_cont_names{con});
    end
end

j=0;
for i=1:length(P)
    if ~exist(P{i})
        j=j+1;
        filestosmooth{j}=strrep(P{i},'s8con','con');
    end
end

if j>0
    for u=1:j
        spm_smooth(filestosmooth{u},strrep(filestosmooth{u},'con_','s8con_'),[8 8 8],0); 
    end
end

%% second level:
clear matlabbatch
i=1;

cd(fullfile(D, S_effect(1).name,'Aud/loc/stats_s5_without_resting'))
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

for k = 1 : length(S_proj)
    subj_of_int         = path_to_scan_proj{k};
    name_subj_of_int    = S_proj(k).name;
    assoc_subj       = path_to_scan_assoc;
    cov_g1         = vector_cov_proj1(k);
    cov_g2         = vector_cov_assoc1;
    age_vector          = [cov_g1;cov_g2];
    
    for c = 1:length(all_contrast)
        %% subject of interest:
        
        contrast = all_contrast(c).name;
        scans1 = {};
        vol_name = sprintf('%s/%s,1',subj_of_int,contrast);
        scans1 = [scans1;vol_name];
        
        res_dir = fullfile(res_dir_base,name_subj_of_int);
        res_dir = sprintf('%s%s/%s', res_dir, '_vs_assoc_age_covariate', names{c});
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
        
        %% for control group:
        scans2 = {};
        for k = 1:numel(assoc_subj)
            vol_name = fullfile(assoc_subj{k});
            vol_name = sprintf('%s/%s,1', assoc_subj{k}, contrast);
            scans2 = [scans2;vol_name];
        end
        
%         matlabbatch{i}.spm.stats.factorial_design.dir = {res_dir};
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.scans1 = scans1;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.scans2 = scans2;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.dept = 0;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.variance = 1;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.gmsca = 0;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.ancova = 0;
% %         matlabbatch{i}.spm.stats.factorial_design.cov(1).c = age_vector;
% %         matlabbatch{i}.spm.stats.factorial_design.cov(1).cname = 'age';
% %         matlabbatch{i}.spm.stats.factorial_design.cov(1).iCFI = 1;
% %         matlabbatch{i}.spm.stats.factorial_design.cov(1).iCC = 1;
%         
%         matlabbatch{i}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
%         matlabbatch{i}.spm.stats.factorial_design.masking.tm.tm_none = 1;
%         matlabbatch{i}.spm.stats.factorial_design.masking.im = 1;
% %         exp_mask = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks',S_effect(1).name);
% %         exp_mask = sprintf('%s_%s_aud_loc_mask_thr_s5.nii',exp_mask, S_effect(end).name);
%         exp_mask = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/anat/Sujet01_Control21_mean_anat.nii';
%         matlabbatch{i}.spm.stats.factorial_design.masking.em = {exp_mask};
%         matlabbatch{i}.spm.stats.factorial_design.globalc.g_omit = 1;
%         matlabbatch{i}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
%         matlabbatch{i}.spm.stats.factorial_design.globalm.glonorm = 1;
%         
%         i = i+1;
        
%         matlabbatch{i}.spm.stats.fmri_est.spmmat = {fullfile(res_dir,'SPM.mat')};
%         matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
%         matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
%         
%         i=i+1;
        
        matlabbatch{i}.spm.stats.con.spmmat = {fullfile(res_dir, 'SPM.mat')};
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'Positive';
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{i}.spm.stats.con.consess{2}.tcon.name = 'Negative';
        matlabbatch{i}.spm.stats.con.consess{2}.tcon.weights = -1;
        matlabbatch{i}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{i}.spm.stats.con.consess{3}.tcon.name = 'subject - all others';
        matlabbatch{i}.spm.stats.con.consess{3}.tcon.weights = [1 -1];
        matlabbatch{i}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{i}.spm.stats.con.delete = 1;
        
        i = i+1;
    end
end

%%
cd '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts'

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'aud_glm';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

% job_ending_rountines(matlabbatch, [], par);

%% part to compare one control vs other controls

% clear matlabbatch
% i=1;

for k = 1 : length(S_assoc)
    subj_of_int         = path_to_scan_assoc{k};
    name_subj_of_int    = S_assoc(k).name;
    assoc_subj       = path_to_scan_assoc;
    assoc_subj(k) = '';    
    cov_g1         = vector_cov_assoc1(k);
    cov_g2         = vector_cov_assoc1;
    cov_g2(k)      = '';
    age_vector          = [cov_g1;cov_g2];
    
    for c = 1:length(all_contrast)
        %% subject of interest:
        
        contrast = all_contrast(c).name;
        scans1 = {};
        vol_name = sprintf('%s/%s,1',subj_of_int,contrast);
        scans1 = [scans1;vol_name];
        
        res_dir = fullfile(res_dir_base,name_subj_of_int);
        res_dir = sprintf('%s%s/%s', res_dir, '_vs_assoc_age_covariate', names{c});
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
        
        %% for control group:
        scans2 = {};
        for k = 1:numel(assoc_subj)
            vol_name = fullfile(assoc_subj{k});
            vol_name = sprintf('%s/%s,1', assoc_subj{k}, contrast);
            scans2 = [scans2;vol_name];
        end
        
%         matlabbatch{i}.spm.stats.factorial_design.dir = {res_dir};
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.scans1 = scans1;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.scans2 = scans2;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.dept = 0;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.variance = 1;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.gmsca = 0;
%         matlabbatch{i}.spm.stats.factorial_design.des.t2.ancova = 0;
% %         matlabbatch{i}.spm.stats.factorial_design.cov(1).c = age_vector;
% %         matlabbatch{i}.spm.stats.factorial_design.cov(1).cname = 'age';
% %         matlabbatch{i}.spm.stats.factorial_design.cov(1).iCFI = 1;
% %         matlabbatch{i}.spm.stats.factorial_design.cov(1).iCC = 1;
%         
%         matlabbatch{i}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
%         matlabbatch{i}.spm.stats.factorial_design.masking.tm.tm_none = 1;
%         matlabbatch{i}.spm.stats.factorial_design.masking.im = 1;
% %         exp_mask = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks',S_effect(1).name);
% %         exp_mask = sprintf('%s_%s_aud_loc_mask_thr_s5.nii',exp_mask, S_effect(end).name);
%         exp_mask = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/anat/Sujet01_Control21_mean_anat.nii';
%         matlabbatch{i}.spm.stats.factorial_design.masking.em = {exp_mask};
%         matlabbatch{i}.spm.stats.factorial_design.globalc.g_omit = 1;
%         matlabbatch{i}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
%         matlabbatch{i}.spm.stats.factorial_design.globalm.glonorm = 1;
%         
%         i = i+1;
        
%         matlabbatch{i}.spm.stats.fmri_est.spmmat = {fullfile(res_dir,'SPM.mat')};
%         matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
%         matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
%         
%         i=i+1;
% 
        matlabbatch{i}.spm.stats.con.spmmat = {fullfile(res_dir, 'SPM.mat')};
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'Positive';
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{i}.spm.stats.con.consess{2}.tcon.name = 'Negative';
        matlabbatch{i}.spm.stats.con.consess{2}.tcon.weights = -1;
        matlabbatch{i}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{i}.spm.stats.con.consess{3}.tcon.name = 'subject - all others';
        matlabbatch{i}.spm.stats.con.consess{3}.tcon.weights = [1 -1];
        matlabbatch{i}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{i}.spm.stats.con.delete = 1;
        
        i = i+1;
    end
end

%%
cd '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts'

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'aud_glm';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

