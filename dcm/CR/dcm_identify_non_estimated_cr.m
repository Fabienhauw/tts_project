% this script is to estimate only not estimated dcm models.
clear;clc;
addpath('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/CR');
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB/spm12'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;
D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'Sujet08', 'FLC', 'FLC (copy)', 'FLC_func'});
S(mask) = [];
subjs={};
nb_subj=numel(S);
path_to_scan = {};
path_to_all_stats = {};
counter = 0;
reg         = {'stgm', 'smg', 'mtg mvpa', 'vwfa'};
stgm=1; smg=2; mtg=3; vwfa=4;
nregions    = length(reg);
nconditions = 2;
sounds=1; words=2;

for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/stats/chaperon/dcm_model_param_modul_endroit');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
end

D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/rosso_controls';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'README'});
S(mask) = [];

nb_subj=nb_subj + numel(S);
for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/stats/chaperon/stc/dcm_model_param_modul_endroit');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
end

all_DCM_non_estimated = {};
for k = 1:35
    disp(path_to_scan{k})
    path_to_stats = path_to_all_stats{k};
    res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models');
    cd(res_path)
    if ~isdir(res_path)
        mkdir(res_path)
    end
    all_DCM = dir('DCM*');
    for dcm = 1 :length(all_DCM)
%         try
%             tmp_dcm = load(all_DCM(dcm).name);
%         catch
%             disp(['error to load: ' fullfile(res_path,all_DCM(dcm).name)])
%             
%             clear DCM
%             
%             SPM = load(fullfile(path_to_stats,'SPM.mat'));
%             SPM = SPM.SPM;
%             
%             reg = {'stgm', 'smg', 'mtg mvpa', 'vwfa'};
%             nregions    = length(reg);
%             
%             stgm=1; smg=2; mtg=3; vwfa=4;
%             
%             f = {
%                 fullfile(path_to_stats,'dcm_minus_Sujet08/VOI_dcm_STGm_perso_1.mat');
%                 fullfile(path_to_stats,'dcm_minus_Sujet08/VOI_dcm_SMG_1.mat');
%                 fullfile(path_to_stats,'dcm_minus_Sujet08/VOI_dcm_mtg_mvpa_1.mat');
%                 fullfile(path_to_stats,'dcm_minus_Sujet08/VOI_dcm_vwfa_sphere_1.mat')
%                 };
%             
%             clear xY
%             tmp_xY = [];
%             for r = 1:length(reg)
%                 XY = load(f{r});
%                 tmp_xY = [tmp_xY;XY];
%             end
%             xY = [tmp_xY.xY];
% 
%             s = struct();
%             s.name       = 'to change';
%             s.u          = include;                 % Conditions
%             s.delays     = repmat(TR/2,1,nregions)';   % Slice timing for each region
%             s.TE         = TE;
%             s.nonlinear  = false;
%             s.two_state  = false;
%             s.stochastic = false;
%             s.centre     = false;
%             s.induced    = 0;
%             s.c          = c;
%             s.d          = d;
%             
%             name_count      = 0;
%             mod_name        = all_DCM(dcm).name;
%             numbers         = regexp(all_DCM(dcm).name, '\d+', 'match');
%             structur        = str2num(numbers{1});
%             modul           = str2num(numbers{2});
%             s.a             = a{1,structur};
%             s.b(:,:,words)  = a{modul+1,structur};
%             s.name          = mod_name;
%             DCM = spm_dcm_specify(SPM,xY,s);            
%         end
        try
            tmp_dcm = load(all_DCM(dcm).name);
        catch
            disp(['still impossible to load: ' fullfile(res_path,all_DCM(dcm).name) ', try to estimate it manually...'])
        end
        if ~isfield(tmp_dcm, 'F')
            all_DCM_non_estimated = [all_DCM_non_estimated;fullfile(res_path,all_DCM(dcm).name)];
            disp(fullfile(res_path,all_DCM(dcm).name))
        end
    end
end

for dcm_mod = 1:length(all_DCM_non_estimated)
    counter = counter + 1;
    matlabbatch{counter}.spm.dcm.estimate.dcms.subj.dcmmat = all_DCM_non_estimated(dcm_mod);
    matlabbatch{counter}.spm.dcm.estimate.output.separate = struct([]);
end
cd('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/CR')

par.run = 1;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '01:30:00';
par.jobname  = 'dcm_cr_estimate';
%%%%%%%% this line below to comment to avoid re writing estimating jobs

% job_ending_rountines(matlabbatch, [], par);
