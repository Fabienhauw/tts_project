clc;clear;

addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))

wd = pwd;

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];


res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/MVPA/Aud/loc';

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet'))));
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control07|Control17

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
    0; 0; 0; 0; 1; 0; 1; 0; 0; 0; 1; 0; 0; 1; 0; 1; 0; 0; 0; 0; 0; 0; ... % end of synesthetes
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

vector_cov1 = vector_age(mask_cov==1);
vector_cov2 = vector_hand(mask_cov==1);

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
erase_model = input('Do you want to erase previous auditive models? [yes/no] ', 's');
if isequal(erase_model,'yes')
    redo_model = 1;
elseif isequal(erase_model,'no')
    redo_model = 0;
end

%% first, get M2 subj scans
scans1={};
smoothcontrast = 0;

if ~exist ('i', 'var')
    i=1;
end

scans_to_smooth = {};

con_dir = sprintf('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/mvpa/',S(1).name);
cd(con_dir)
comparisons = dir('results*');

a = 1; b = numel(S);

for tmp_comp = 1 : length(comparisons)
    for k = a:b %for each subject;
        con_dir = sprintf('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/mvpa/',S(k).name);
        cd(fullfile(con_dir,comparisons(tmp_comp).name));
        all_smooth_contrast = dir('s*minus_chance.nii*');
        if ~exist('redo_scon','var')
            redo_scon=0;
        end
        if redo_scon & ~isempty (all_smooth_contrast)
            try
                delete (all_smooth_contrast.name)
                all_smooth_contrast={};
            end
        end
        
        all_contrast = dir('*minus_chance.nii');
        allsmooth = {};
        for asc = 1:length(all_smooth_contrast)
            smoothcon = all_smooth_contrast(asc).name;
            allsmooth = [allsmooth;smoothcon];
        end
        
        allcon = {};
        for ac = 1:length(all_contrast)
            con = all_contrast(ac).name;
            allcon = [allcon;con];
        end
        
        desmoothcon={};
        if ~isempty(allsmooth)
            for y=1:length(allsmooth)
                desmoothcon(y,:) = strrep(allsmooth(y,:),'sres_AUC*.nii','res_AUC*.nii');
            end
        end
        if ~isempty(desmoothcon)
            contrast = setdiff(allcon,desmoothcon);
        elseif isempty(desmoothcon) | ~exist(desmoothcon)
            contrast = allcon;
        end
        
        for c = 1:length(contrast)
            con = contrast(c);
            vol_name=fullfile(con_dir, comparisons(tmp_comp).name, con);
            scans_to_smooth=[scans_to_smooth;vol_name];
        end
    end
end

matlabbatch{i}.spm.spatial.smooth.data = scans_to_smooth;
matlabbatch{i}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{i}.spm.spatial.smooth.dtype = 0;
matlabbatch{i}.spm.spatial.smooth.im = 0;
matlabbatch{i}.spm.spatial.smooth.prefix = 's';
i=i+1;


cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/mvpa_scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '01:00:00'; % usually takes 7-8h on cluster for 124 permutations x 8-fold classifier
par.jobname  = 'smoothing_of_con';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

return


%% second level:
clear matlabbatch
i=1;

for tmp_comp = 1:length(comparisons)
    %% subject of interest:
    scans1 = {};
    for k = 1 : numel(S_droit)
        mvpa_res_dir = fullfile(D, S_droit(k).name, 'Aud/loc/mvpa', comparisons(tmp_comp).name);
        cd(mvpa_res_dir);
        contrast = dir('s*_minus_chance.nii');
        vol_name = fullfile(mvpa_res_dir, contrast(1).name);
        scans1 = [scans1 ; vol_name];
    end
    
    res_dir = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/MVPA/Aud/loc/syn_vs_con', comparisons(tmp_comp).name);
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
    
    %% for group, M2 controls:
    scans2 = {};
    for k = 1 : numel(S_con_app)
        mvpa_res_dir = fullfile(D, S_con_app(k).name, 'Aud/loc/mvpa', comparisons(tmp_comp).name);
        cd(mvpa_res_dir);
        contrast = dir('res*_minus_chance.nii');
        vol_name = fullfile(mvpa_res_dir, contrast(1).name);
        scans2 = [scans2 ; vol_name];
    end
    
    matlabbatch{i}.spm.stats.factorial_design.dir = {res_dir};
    matlabbatch{i}.spm.stats.factorial_design.des.t2.scans1 = scans1;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.scans2 = scans2;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.dept = 0;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.variance = 1;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.gmsca = 0;
    matlabbatch{i}.spm.stats.factorial_design.des.t2.ancova = 0;
    
    matlabbatch{i}.spm.stats.factorial_design.cov(1).c = vector_cov1;
    matlabbatch{i}.spm.stats.factorial_design.cov(1).cname = 'age';
    matlabbatch{i}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{i}.spm.stats.factorial_design.cov(1).iCC = 1;
    
    matlabbatch{i}.spm.stats.factorial_design.cov(2).c = vector_cov2;
    matlabbatch{i}.spm.stats.factorial_design.cov(2).cname = 'handedness';
    matlabbatch{i}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{i}.spm.stats.factorial_design.cov(2).iCC = 1;
    
    %     matlabbatch{i}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{i}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{i}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{i}.spm.stats.factorial_design.masking.im = 1;
    exp_mask = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/MVPA/mask',S_effect(1).name);
    exp_mask = sprintf('%s_%s_aud_mask_thresholded.nii',exp_mask, S_effect(end).name);
    matlabbatch{i}.spm.stats.factorial_design.masking.em = {exp_mask};
    matlabbatch{i}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{i}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{i}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    i = i+1;
end

%%

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/mvpa_scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00'; % usually takes 7-8h on cluster for 124 permutations x 8-fold classifier
par.jobname  = 'mvpa_second_lvl_aud_spec';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

    
    %%
clear matlabbatch
i=1;

for tmp_comp = 1:length(comparisons)
    res_dir = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/MVPA/Aud/loc/syn_vs_con', comparisons(tmp_comp).name);
    matlabbatch{i}.spm.stats.fmri_est.spmmat = {fullfile(res_dir,'SPM.mat')};
    matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
    
    i=i+1;

end

%%

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/mvpa_scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00'; % usually takes 7-8h on cluster for 124 permutations x 8-fold classifier
par.jobname  = 'mvpa_second_lvl_aud_est';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

%%

clear matlabbatch
i=1;

for tmp_comp = 1:length(comparisons)

    res_dir = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/MVPA/Aud/loc/syn_vs_con', comparisons(tmp_comp).name);
    matlabbatch{i}.spm.stats.con.spmmat = {fullfile(res_dir, 'SPM.mat')};
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'syn';
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.name = 'con';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.weights = [0 1];
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.name = 'syn-con';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.weights = [1 -1];
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.name = 'con-syn';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.weights = [-1 1];
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.delete = 1;
    
    i = i+1;
end

%%

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/mvpa_scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00'; % usually takes 7-8h on cluster for 124 permutations x 8-fold classifier
par.jobname  = 'mvpa_second_lvl_aud_con';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);


% spm_jobman('interactive', matlabbatch)
