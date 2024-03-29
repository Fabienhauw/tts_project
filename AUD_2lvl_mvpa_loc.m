clear; clc;
addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))

wd = pwd;

res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/mvpa_without_resting_ttest_s8';
% res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/MVPA/Aud/loc/syn_vs_con_rh_s5';

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

radius = 10;

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17 

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
% gaucher_appar = {'Control02|Control04|Control07|Control17|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con = S;
S_con(mask_gauch_con) = [];

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_syn = S;
S_syn(mask_gauch) = [];

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

S_effect = [S_syn ; S_con];
subname = {S_effect.name};

totsub=length(subname);

for j = 1 : size(S,1)
    if ~isempty(find(~cellfun(@isempty,(regexp({S_effect.name},S(j).name)))))
        mask_cov(j,1) = 1;
    else
        mask_cov(j,1) = 0;
    end
end

vector_cov1 = vector_age(mask_cov==1);
vector_cov2 = vector_hand(mask_cov==1);

con_dir = sprintf('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/mvpa_without_resting_%dmm/',S_effect(1).name, radius);
cd(con_dir)
results = dir('results_Scr*NormalSpeech');
for tmp_comp = 1 : length(results)
    str_path = strsplit(results(tmp_comp).name, '_');
    all_comp_inv{tmp_comp} = strjoin({str_path{1}, str_path{4}, str_path{3}, str_path{2}}, '_');
    all_comp{tmp_comp} = results(tmp_comp).name;
end
ncon = 1;
all_comp_names = all_comp;
cont = 1; %[1 : length(all_comp)];

% get image files names
P={};
for con=1:ncon
    for s=1:totsub
        sub=subname{s};
%         P{(con-1)*totsub+s} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/mvpa_without_resting_%dmm/%s/res_accuracy_minus_chance.nii'],sub, radius, all_comp_names{con});
%         P{(con-1)*totsub+s} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/mvpa_without_resting_%dmm/%s/s6res_accuracy_minus_chance.nii'],sub, radius, all_comp_names{con});
        P{(con-1)*totsub+s} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/mvpa_without_resting_%dmm/%s/s8res_accuracy_minus_chance.nii'],sub, radius, all_comp_names{con});
    end
end

j=0;
for i=1:length(P)
    if ~exist(P{i})
        j=j+1;
%         filestosmooth{j}=strrep(P{i},'s6res','res');
        filestosmooth{j}=strrep(P{i},'s8res','res');
    end
end

if j>0
    for u=1:j
%         spm_smooth(filestosmooth{u},strrep(filestosmooth{u},'res_','s6res_'),[6 6 6],0); 
        spm_smooth(filestosmooth{u},strrep(filestosmooth{u},'res_','s8res_'),[8 8 8],0); 
    end
end

clear matlabbatch
i=1;

%% to erase previous smoothed contrasts;
erase_scon = sprintf(input('Do you want to erase previous auditive scon? [yes/no] ', 's'));
if isequal(erase_scon,'yes')
    redo_scon = 1;
elseif isequal(erase_scon,'no')
    redo_scon = 0;
end

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

for tmp_comp = 1 : length(all_comp)
    %% group 1:
    scans1 = {};
    for k = 1:length(S_syn)
        if ~isdir(fullfile(D, S_syn(k).name,sprintf('Aud/loc/mvpa_without_resting_%dmm', radius),all_comp_inv{tmp_comp}))
            vol_name = fullfile(D, S_syn(k).name,sprintf('Aud/loc/mvpa_without_resting_%dmm', radius),all_comp{tmp_comp}, 's8res_accuracy_minus_chance.nii');
        else
            vol_name = fullfile(D, S_syn(k).name,sprintf('Aud/loc/mvpa_without_resting_%dmm', radius),all_comp_inv{tmp_comp}, 's8res_accuracy_minus_chance.nii');
        end
        scans1 = [scans1;vol_name];
    end
    res_dir = fullfile(res_dir_base,all_comp{tmp_comp});
%     res_dir = sprintf('%s/smooth_%s', res_dir_base,all_comp{tmp_comp});

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
    
    %% for group 2:
    scans2 = {};
    for k = 1:length(S_con) %S_con if all controls
        if ~isdir(fullfile(D, S_con(k).name,sprintf('Aud/loc/mvpa_without_resting_%dmm', radius),all_comp_inv{tmp_comp}))
            vol_name = fullfile(D, S_con(k).name,sprintf('Aud/loc/mvpa_without_resting_%dmm', radius),all_comp{tmp_comp}, 's8res_accuracy_minus_chance.nii');
        else
            vol_name = fullfile(D, S_con(k).name,sprintf('Aud/loc/mvpa_without_resting_%dmm', radius),all_comp_inv{tmp_comp}, 's8res_accuracy_minus_chance.nii');
        end
        scans2 = [scans2;vol_name];
    end
    
    matlabbatch{i}.spm.stats.factorial_design.dir = {res_dir};
    matlabbatch{i}.spm.stats.factorial_design.des.t2.scans1 = scans1;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.scans2 = scans2;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.dept = 0;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.variance = 1;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.gmsca = 0;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.ancova = 0;
        
%     matlabbatch{i}.spm.stats.factorial_design.cov(1).c = vector_cov1;
%     matlabbatch{i}.spm.stats.factorial_design.cov(1).cname = 'age';
%     matlabbatch{i}.spm.stats.factorial_design.cov(1).iCFI = 1;
%     matlabbatch{i}.spm.stats.factorial_design.cov(1).iCC = 1;
%     
%     matlabbatch{i}.spm.stats.factorial_design.cov(2).c = vector_cov2;
%     matlabbatch{i}.spm.stats.factorial_design.cov(2).cname = 'handedness';
%     matlabbatch{i}.spm.stats.factorial_design.cov(2).iCFI = 1;
%     matlabbatch{i}.spm.stats.factorial_design.cov(2).iCC = 1;
    
    matlabbatch{i}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{i}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{i}.spm.stats.factorial_design.masking.im = 1;
    exp_mask = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/MVPA/mask',S_effect(1).name);
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
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'Synesthetes';
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.name = 'Controls';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.weights = [0 1];
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.name = '-Synesthetes';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.weights = -1;
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.name = '-Controls';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{5}.tcon.name = 'syn-con';
    matlabbatch{i}.spm.stats.con.consess{5}.tcon.weights = [1 -1];
    matlabbatch{i}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{6}.tcon.name = 'syn+con';
    matlabbatch{i}.spm.stats.con.consess{6}.tcon.weights = [1 1];
    matlabbatch{i}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{7}.tcon.name = 'con-syn';
    matlabbatch{i}.spm.stats.con.consess{7}.tcon.weights = [-1 1];
    matlabbatch{i}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    
    matlabbatch{i}.spm.stats.con.delete = 1;
    
    i = i+1;
end

spm_jobman('run', matlabbatch)

cd(res_dir)