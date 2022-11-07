%%batch for VOI mSTG, SMG, VWFA, IFG for auditive run.

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

ROIs_names = {'SMG'; 'VWFA'; 'mIFG'; 'mSTG'}; % VOI mSTG, resulting from all subj, scr speech > silence
ROIs_coord = [-58 -44 23; -52 -51 -20; -45 22 23; -52 -14 6];

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% gaucher_appar = {'Control02|Control04|Control07|Control17|Sujet'};
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));

mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con = S;
S_con(mask_gauch_con) = [];

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

for j = 1 : size(S,1)
    if ~isempty(find(~cellfun(@isempty,(regexp({S_effect.name},S(j).name)))))
        mask_cov(j,1) = 1;
    else
        mask_cov(j,1) = 0;
    end
end

vector_cov1 = vector_age(mask_cov==1);
vector_cov2 = vector_hand(mask_cov==1);

i=1;

for k = 1 : numel(S_effect)
    for tmp_ROI = 1 : numel(ROIs_names)
    %% batch for VOI SMG, resulting from all subj, norm>scr speech
    matlabbatch{i}.spm.util.voi.spmmat = {fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5/SPM.mat')}; % use this .mat to extract time series.
    matlabbatch{i}.spm.util.voi.adjust = NaN;
    matlabbatch{i}.spm.util.voi.session = 1;
    matlabbatch{i}.spm.util.voi.name = sprintf('%s_PPI_adapted', ROIs_names{tmp_ROI});
    matlabbatch{i}.spm.util.voi.roi{1}.spm.spmmat = {fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5/SPM.mat')}; % use this .mat to deduct the highest activated voxel.
    matlabbatch{i}.spm.util.voi.roi{1}.spm.contrast = 11;
    matlabbatch{i}.spm.util.voi.roi{1}.spm.conjunction = 1;
    matlabbatch{i}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
    matlabbatch{i}.spm.util.voi.roi{1}.spm.thresh = 1;
    matlabbatch{i}.spm.util.voi.roi{1}.spm.extent = 0;
    matlabbatch{i}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
    matlabbatch{i}.spm.util.voi.roi{2}.sphere.centre = ROIs_coord(tmp_ROI,:);
    matlabbatch{i}.spm.util.voi.roi{2}.sphere.radius = 6;
    matlabbatch{i}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
    matlabbatch{i}.spm.util.voi.roi{3}.sphere.centre = [0 0 0];
    matlabbatch{i}.spm.util.voi.roi{3}.sphere.radius = 4;
    matlabbatch{i}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
    matlabbatch{i}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';
    matlabbatch{i}.spm.util.voi.expression = 'i1&i3';
    i=i+1;
    end
end

%%

% spm_jobman('run', matlabbatch)

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'ppi_voi';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

%%
clear matlabbatch
i=1;
ROI_files = {
    'spmT_0011_best_vox_sph_-58_-44_23_based_on_auditive_con_11.nii';
    'spmT_0011_best_vox_sph_-52_-51_-20_based_on_auditive_con_11.nii';
    'spmT_0011_best_vox_sph_-45_22_23_based_on_auditive_con_11.nii';
    'spmT_0011_best_vox_sph_-52_-14_6_based_on_auditive_con_11.nii';
    };

for k = 1 : numel(S_effect)
    for tmp_ROI = 1 : numel(ROI_files)
    %% batch for VOI SMG, resulting from all subj, norm>scr speech
    matlabbatch{i}.spm.util.voi.spmmat = {fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5/SPM.mat')}; % use this .mat to extract time series.
    matlabbatch{i}.spm.util.voi.adjust = NaN;
    matlabbatch{i}.spm.util.voi.session = 1;
    matlabbatch{i}.spm.util.voi.name = sprintf('%s_best_vox_PPI_adapted', ROIs_names{tmp_ROI});
    matlabbatch{i}.spm.util.voi.roi{1}.mask.image = {fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5', ROI_files{tmp_ROI})};
    matlabbatch{i}.spm.util.voi.roi{1}.mask.threshold = 0.5;
    matlabbatch{i}.spm.util.voi.expression = 'i1';
    i=i+1;
    end
end

%%

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'ppi_voi';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

%%
clear matlabbatch
i=1;
for k = 1 : numel(S_effect)
    for tmp_ROI = 1 : numel(ROIs_names)
        %% batch for PPI ROI X (norm-scr)
        
        SPM_dir = fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5/SPM.mat');
        
        matlabbatch{i}.spm.stats.ppi.spmmat = {SPM_dir};
        matlabbatch{i}.spm.stats.ppi.type.ppi.voi = {fullfile(D, S_effect(k).name, sprintf('/Aud/loc/stats_s5/VOI_%s_PPI_adapted_1.mat', ROIs_names{tmp_ROI}))};
        matlabbatch{i}.spm.stats.ppi.type.ppi.u = [4 1 1
            5 1 -1];
        matlabbatch{i}.spm.stats.ppi.name = sprintf('%s_adaptx(norm-scr)', ROIs_names{tmp_ROI});
        matlabbatch{i}.spm.stats.ppi.disp = 0; %1 to display
        i=i+1;
    end
end

%%
for k = 1 : numel(S_effect)
    for tmp_ROI = 1 : numel(ROI_files)
        %% batch for PPI ROI X (norm-scr)
        
        SPM_dir = fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5/SPM.mat');
        
        matlabbatch{i}.spm.stats.ppi.spmmat = {SPM_dir};
        matlabbatch{i}.spm.stats.ppi.type.ppi.voi = {fullfile(D, S_effect(k).name, sprintf('/Aud/loc/stats_s5/VOI_%s_best_vox_PPI_adapted_1.mat', ROIs_names{tmp_ROI}))};
        matlabbatch{i}.spm.stats.ppi.type.ppi.u = [4 1 1
            5 1 -1];
        matlabbatch{i}.spm.stats.ppi.name = sprintf('%s_best_vox_adaptx(norm-scr)', ROIs_names{tmp_ROI});
        matlabbatch{i}.spm.stats.ppi.disp = 0; %1 to display
        i=i+1;
    end
end
%%
cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 1;
par.sge = 0;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'ppi_voi_x_conds';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);

%%
clear matlabbatch
i=1;
for k = 1 : numel(S_effect)
    for tmp_ROI = 1 : numel(ROIs_names)
        %% PPI model with ROI VOI
        PPI_SPM_dir = sprintf('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/PPI/PPI_%s/PPI_%s_norm_scr', S_effect(k).name, ROIs_names{tmp_ROI});
        if ~isdir(PPI_SPM_dir)
            mkdir(PPI_SPM_dir)
        end
        
        filename = fullfile(D,S_effect(k).name,'Aud/loc/param');
        cd (filename);
        json=dir('*.json');
        json=json.name;
        
        res = get_string_from_json(json, {'RepetitionTime'}, {'num'});
        if res{1}>100
            TR          = res{1}/1000; % millisecond -> second
        else
            TR          = res{1}; % second
        end
        
        scans={};
        cd(fullfile(D, S_effect(k).name,'Aud/loc/swf'))
        vol_name = dir('s5*wts_OC.nii');
        vol_name=fullfile(D, S_effect(k).name,'Vis/loc/swf', vol_name(1).name);
        nb_vol = size(spm_vol(vol_name),1);
        
        
        matlabbatch{i}.spm.stats.fmri_spec.dir = {PPI_SPM_dir};
        matlabbatch{i}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{i}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        for v=1:nb_vol
            volume=sprintf('%s,%d',vol_name,v);
            scans=[scans;volume];
        end
        
        cd (fullfile(D,S_effect(k).name,'Aud/loc/stats_s5'))
        load(sprintf('PPI_%s_adaptx(norm-rev)', ROIs_names{tmp_ROI}));
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.scans = scans;
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{i}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(1).name = 'PPI-interaction';
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(1).val = PPI.ppi;
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(2).name = sprintf('%s-BOLD', ROIs_names{tmp_ROI});
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(2).val = PPI.Y;
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(3).name = 'Psych_norm-scr';
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(3).val = PPI.P;
        
        multi_reg=fullfile(D,S_effect(k).name,'Aud/loc/param');
        cd (multi_reg);
        mr = dir('multiple*.txt');
        mr = mr.name;
        multi_reg = fullfile(multi_reg,mr);
        matlabbatch{i}.spm.stats.fmri_spec.sess.multi_reg = {multi_reg};
        matlabbatch{i}.spm.stats.fmri_spec.sess.hpf = 192;
        matlabbatch{i}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{i}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{i}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{i}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{i}.spm.stats.fmri_spec.mthresh = 0;
        mask = fullfile(D,S_effect(k).name,'anat/brain_extraction_mask.nii');
        matlabbatch{i}.spm.stats.fmri_spec.mask = {mask};
        matlabbatch{i}.spm.stats.fmri_spec.cvi = 'AR(1)';
        i=i+1;
        
        matlabbatch{i}.spm.stats.fmri_est.spmmat(1) = fullfile(PPI_SPM_dir, 'SPM.mat');
        matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
        i=i+1;
        
        matlabbatch{i}.spm.stats.con.spmmat(1) = fullfile(PPI_SPM_dir, 'SPM.mat');
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'PPI-interaction';
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{i}.spm.stats.con.delete = 0;
        
        clear(sprintf('PPI_%s_adaptx(norm-scr)', ROIs_names{tmp_ROI}));
    end
end

%%
cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 1;
par.sge = 0;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '01:00:00';
par.jobname  = 'ppi_glm_models';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);


%% batch for each VOI x each condition separately
clear matlabbatch
i=1;
for k = 1 : numel(S_effect)
    for tmp_ROI = 1 : numel(ROIs_names)
        SPM_dir = fullfile(D, S_effect(k).name, 'Aud/loc/stats_s5/SPM.mat');
        
        matlabbatch{i}.spm.stats.ppi.spmmat = {SPM_dir};
        matlabbatch{i}.spm.stats.ppi.type.ppi.voi = {fullfile(D, S_effect(k).name, sprintf('/Aud/loc/stats_s5/VOI_%s_PPI_adapted_1.mat', ROIs_names{tmp_ROI}))};
        matlabbatch{i}.spm.stats.ppi.type.ppi.u = [4 1 1];
        matlabbatch{i}.spm.stats.ppi.name = sprintf('%s_adapt_x_norm', ROIs_names{tmp_ROI});
        matlabbatch{i}.spm.stats.ppi.disp = 0;
        i=i+1;
        
        matlabbatch{i}.spm.stats.ppi.spmmat = {SPM_dir};
        matlabbatch{i}.spm.stats.ppi.type.ppi.voi = {fullfile(D, S_effect(k).name, sprintf('/Aud/loc/stats_s5/VOI_%s_PPI_adapted_1.mat', ROIs_names{tmp_ROI}))};
        matlabbatch{i}.spm.stats.ppi.type.ppi.u = [5 1 1];
        matlabbatch{i}.spm.stats.ppi.name = sprintf('%s_adapt_x_scr', ROIs_names{tmp_ROI});
        matlabbatch{i}.spm.stats.ppi.disp = 0;
        i=i+1;
    end
end

%%
cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'ppi_voi_x_conds';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);


%% PPI model with SMG VOI
clear matlabbatch
i=1;
for k = 1 : numel(S_effect)
    for tmp_ROI = 1 : numel(ROIs_names)
        PPI_SPM_dir = sprintf('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/PPI/PPI_%s/PPI_%s_norm', S_effect(k).name, ROIs_names{tmp_ROI});
        if ~isdir(PPI_SPM_dir)
            mkdir(PPI_SPM_dir)
        end
        
        filename = fullfile(D,S_effect(k).name,'Aud/loc/param');
        cd (filename);
        json=dir('*.json');
        json=json.name;
        
        res = get_string_from_json(json, {'RepetitionTime'}, {'num'});
        if res{1}>100
            TR          = res{1}/1000; % millisecond -> second
        else
            TR          = res{1}; % second
        end
        
        scans={};
        cd(fullfile(D, S_effect(k).name,'Aud/loc/swf'))
        vol_name = dir('s5*wts_OC.nii');
        vol_name=fullfile(D, S_effect(k).name,'Vis/loc/swf', vol_name(1).name);
        nb_vol = size(spm_vol(vol_name),1);
        
        matlabbatch{i}.spm.stats.fmri_spec.dir = {PPI_SPM_dir};
        matlabbatch{i}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{i}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        for v=1:nb_vol
            volume=sprintf('%s,%d',vol_name,v);
            scans=[scans;volume];
        end
        
        cd (fullfile(D,S_effect(k).name,'Aud/loc/stats_s5'))
        load(sprintf('PPI_%s_adapt_x_norm', ROIs_names{tmp_ROI}));
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.scans = scans;
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{i}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(1).name = 'PPI-interaction';
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(1).val = PPI.ppi;
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(2).name = sprintf('%s-BOLD', ROIs_names{tmp_ROI});
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(2).val = PPI.Y;
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(3).name = 'Psych_norm';
        
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress(3).val = PPI.P;
        
        multi_reg=fullfile(D,S_effect(k).name,'Aud/loc/param');
        cd (multi_reg);
        mr = dir('multiple*.txt');
        mr = mr.name;
        multi_reg=fullfile(multi_reg,mr);
        matlabbatch{i}.spm.stats.fmri_spec.sess.multi_reg = {''}; %{multi_reg};
        matlabbatch{i}.spm.stats.fmri_spec.sess.hpf = 192;
        matlabbatch{i}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{i}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{i}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{i}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{i}.spm.stats.fmri_spec.mthresh = 0;
        mask = fullfile(D,S_effect(k).name,'anat/brain_extraction_mask.nii');
        matlabbatch{i}.spm.stats.fmri_spec.mask = {mask};
        matlabbatch{i}.spm.stats.fmri_spec.cvi = 'AR(1)';
        i=i+1;
        
        matlabbatch{i}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
        i=i+1;
        
        matlabbatch{i}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'PPI-interaction';
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{i}.spm.stats.con.delete = 0;
        i=i+1;
    end
end
