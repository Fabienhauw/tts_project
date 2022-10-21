clear
addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))

wd = pwd;

res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/syn_vs_con_rh_s5_with_pca_score';
% res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/syn_vs_con_s5';

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
% gaucher_appar = {'Control02|Control04|Control07|Control17|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con_app = S;
S_con_app(mask_gauch_con) = [];

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_droit = S;
S_droit(mask_gauch) = [];

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

S_effect = [S_droit ; S_con_app];

for j = 1 : size(S,1)
    if ~isempty(find(~cellfun(@isempty,(regexp({S_effect.name},S(j).name)))))
        mask_cov(j,1) = 1;
    else
        mask_cov(j,1) = 0;
    end
end

load ('/home/fabien.hauw/Desktop/Fabien/NeoTopLex/+tts_group/cpt_data/pca_score.mat')
vector_cov1 = score_pca;
vector_cov2 = vector_age(mask_cov==1);
vector_cov3 = vector_hand(mask_cov==1);

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

cd(fullfile(D, S(1).name,'Aud/loc/stats_s5'))
all_contrast = dir(sprintf('con*.nii'));
all_contrast = all_contrast(~cellfun(@isempty,(regexp({all_contrast.name}, '^con_\d+\.nii'))));

names = {...
            'words', 'pseudowords', 'numbers', 'normal_speech', 'scramble_speech',...
            'odds', 'motor',...
            'lexicality', '-lexicality','(words + normal_speech + numbers) - pseudowords',...
            'phonology', '-phonology', 'numbers - (words + pseudowords)', ...
            'normal_speech - words', 'words + pseudowords + numbers + normal_speech', ...
            'words - normal_speech', 'words - scramble_speech', 'pseudowords - scramble_speech', '(words+pseudowords) -  scramble_speech', ...
            '(words + normal_speech) -  scramble_speech', '(normal_speech + pseudowords) -  scramble_speech', '(normal_speech + pseudowords + words) -  scramble_speech',...
            '(normal_speech + pseudowords + words + numbers) - scramble_speech',...
            'numbers - words', 'words - numbers', ...
            'EOI', ...
            };

for c = 1:length(all_contrast)
    %% subject of interest:
    contrast = all_contrast(c).name;
    scans = {};
    for k = 1:length(S_effect)
        vol_name = fullfile(D, S_effect(k).name,'Aud/loc/stats_s5/');
        vol_name = sprintf('%s%s,1',vol_name,contrast);
        scans = [scans;vol_name];
    end
    res_dir = fullfile(res_dir_base,names{c});
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
    
    matlabbatch{i}.spm.stats.factorial_design.dir = {res_dir};
    matlabbatch{i}.spm.stats.factorial_design.des.t1.scans = scans;
        
    matlabbatch{i}.spm.stats.factorial_design.cov(1).c = vector_cov1;
    matlabbatch{i}.spm.stats.factorial_design.cov(1).cname = 'pca_score';
    matlabbatch{i}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{i}.spm.stats.factorial_design.cov(1).iCC = 1;
    
    matlabbatch{i}.spm.stats.factorial_design.cov(2).c = vector_cov2;
    matlabbatch{i}.spm.stats.factorial_design.cov(2).cname = 'age';
    matlabbatch{i}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{i}.spm.stats.factorial_design.cov(2).iCC = 1;
    
    matlabbatch{i}.spm.stats.factorial_design.cov(3).c = vector_cov3;
    matlabbatch{i}.spm.stats.factorial_design.cov(3).cname = 'handedness';
    matlabbatch{i}.spm.stats.factorial_design.cov(3).iCFI = 1;
    matlabbatch{i}.spm.stats.factorial_design.cov(3).iCC = 1;

    
%     matlabbatch{i}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{i}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{i}.spm.stats.factorial_design.masking.im = 1;
    exp_mask = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks',S_effect(1).name);
    exp_mask = sprintf('%s_%s_aud_loc_mask_thr_s5.nii',exp_mask, S_effect(end).name);
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
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'all';
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.name = '- all';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.name = 'effect of pca';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.weights = [0 1];
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.name = 'effect of - pca';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.delete = 1;
    
    i = i+1;
end

spm_jobman('run', matlabbatch)

cd(res_dir)