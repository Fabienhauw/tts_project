%script to automatize conn RS analysis
clear;clc;

addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB'))
addpath(genpath('/network/lustre/isbackps02/home/fabien.hauw/Documents/MATLAB/conn18.b'))


cwd = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/func_conn';
only_right = 1;
modif_freq = 0;

%% definition of all subjects:
subjs={};
group = [];
path_to_scan = {};
repet_time = [];
STRUCTURAL_FILE = {};
FUNCTIONAL_FILE = {};
ROI_AUD_FILE = {};
ONSET_AUD_FILE = {};
GM = {};
WM = {};
CSF = {};
nsessions = {};

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
S = dir(D);
Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet'))));
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));

S = [Syn;Con];

%%

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17

%%
if only_right
    proj_name = 'syn_group_aud_func_conn_rh_s5';
    gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
    mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
else
    proj_name = 'syn_group_aud_func_conn_s5';
    gaucher_appar = {'Control02|Control04|Control07|Control17|Sujet'};
    mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
end

if modif_freq
    proj_name = sprintf('%s_modif_freq', proj_name);
end

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

S_final = [S_syn ; S_con];

for j = 1 : size(S,1)
    if ~isempty(find(~cellfun(@isempty,(regexp({S_final.name},S(j).name)))))
        mask_cov(j,1) = 1;
    else
        mask_cov(j,1) = 0;
    end
end

vector_cov1 = vector_age(mask_cov==1);
vector_cov2 = vector_hand(mask_cov==1);

S = S_final;

aud_coord_names =   {'SMG',    'VWFA',      'lIPS',     'lprecent',   'MFG',  'lpSTG',   'laSTG',  'rpSTG',   'raSTG'};

for k=1:numel(S)
    if ~isempty(regexp(S(k).name, 'Sujet'))
        if isempty(regexp(S(k).name,'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16'))
            subj_group     = 1; % 1 is for RH synesthetes, 2 for RH controls, 3 for LH syn, 4 for LH con
        else
            subj_group     = 3;
        end
    elseif ~isempty(regexp(S(k).name, 'Control'))
        if isempty(regexp(S(k).name,'Control22|Control23|Control24|Control25|Control26'))
            subj_group     = 2; % 1 is for RH synesthetes, 2 for RH controls, 3 for LH syn, 4 for LH con
        else
            subj_group     = 4;
        end
    end
    nbsuj = numel(S);
    subjs = [subjs;S(k).name];
    group = [group; subj_group];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
    nsession = 1;
    nsessions = [nsessions;nsession];
    
    
    % get the TR:
    cd (fullfile(D,S(k).name,'Aud/loc/param'));
    json=dir('*.json');
    json=json.name;
    res = get_string_from_json(json, {'RepetitionTime'}, {'num'});
    TR          = res{1}; if TR > 100, TR = TR/1000; end
    repet_time = [repet_time, TR];
    
    
    % get the structural scan:
    cd (fullfile(D,S(k).name,'anat'));
    anat = dir('wmv*.nii');
    anat = fullfile(D,S(k).name,'anat',anat.name);
    STRUCTURAL_FILE = [STRUCTURAL_FILE;anat];
    
    tmp_gm = dir('wp1*.nii'); tmp_gm = fullfile(D,S(k).name,'anat',tmp_gm.name);
    GM = [GM; tmp_gm];
    
    tmp_wm =  dir('wp2*.nii'); tmp_wm = fullfile(D,S(k).name,'anat',tmp_wm.name);
    WM = [WM; tmp_wm];
    
    tmp_csf = dir('wp3*.nii'); tmp_csf = fullfile(D,S(k).name,'anat',tmp_csf.name);
    CSF = [CSF; tmp_csf];
    
    %get the functional scan:
    cd (fullfile(D,S(k).name,'Aud/loc/swf'));
    func = dir('s5w*.nii');
    func = fullfile(D,S(k).name,'Aud/loc/swf',func.name);
    FUNCTIONAL_FILE = [FUNCTIONAL_FILE;func];
    
    cd (fullfile(D,S(k).name, 'Aud/loc/stats_s5_without_resting'))
    for tmp_roi = 1 : length(aud_coord_names)
        roi_filename = dir(sprintf('VOI*%s_PPI_adj_eoi5_*con16*mask.nii', aud_coord_names{tmp_roi}));
        roi = fullfile(D,S(k).name, 'Aud/loc/stats_s5_without_resting', roi_filename.name);
        ROI_AUD_FILE{tmp_roi}{k,1} = roi;
    end
    
    cd (fullfile(D,S(k).name, 'Aud/loc/cpt_data'))
    onsets_filename = dir('Timedata*without_resting.mat');
    onsets_dur = fullfile(D,S(k).name, 'Aud/loc/cpt_data', onsets_filename.name);
    ONSET_AUD_FILE = [ONSET_AUD_FILE;onsets_dur];
end

NSUBJECTS = length(subjs);
if rem(length(FUNCTIONAL_FILE),NSUBJECTS),error('mismatch number of functional files %n', length(FUNCTIONAL_FILE));end
if rem(length(STRUCTURAL_FILE),NSUBJECTS),error('mismatch number of anatomical files %n', length(FUNCTIONAL_FILE));end
nsessions       = length(FUNCTIONAL_FILE)/NSUBJECTS;
FUNCTIONAL_FILE = reshape(FUNCTIONAL_FILE,[nsessions, NSUBJECTS]);
STRUCTURAL_FILE = {STRUCTURAL_FILE{1:NSUBJECTS}};
disp([num2str(size(FUNCTIONAL_FILE,1)),' session(s)']);
disp([num2str(size(FUNCTIONAL_FILE,2)),' subject(s)']);
TR              = repet_time';

% CONN-SPECIFIC SECTION: RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
% Prepares batch structure
clear batch;
batch.filename  = fullfile(cwd, proj_name);            % New conn_*.mat experiment name

% SETUP & PREPROCESSING step (using default values for most parameters, see help conn_batch to define non-default values)
% CONN Setup                                            % Default options (uses all ROIs in conn/rois/ directory); see conn_batch for additional options
% CONN Setup.preprocessing                               (realignment/coregistration/segmentation/normalization/smoothing)
batch.Setup.isnew       = 1;
batch.Setup.nsubjects   = NSUBJECTS;
batch.Setup.RT          = TR;                                        % TR (seconds)
batch.Setup.functionals = repmat({{}},[NSUBJECTS,1]);       % Point to functional volumes for each subject/session

for nsub = 1:NSUBJECTS
    for nses=1:nsessions
        batch.Setup.functionals{nsub}{nses}{1}=FUNCTIONAL_FILE{nses,nsub};
    end
end %note: each subject's data is defined by three sessions and one single (4d) file per session

batch.Setup.structurals = STRUCTURAL_FILE;                  % Point to anatomical volumes for each subject

%% covariate part;
%% first lvl
batch.Setup.covariates.names{1, 1}    = 'realignment';

for cov = 1 : length(batch.Setup.covariates.names)
    for subj=1:nbsuj
        cd (fullfile(D,S(subj).name,'Aud/loc/param'));
        rp = dir('multiple_reg*.txt'); rp = rp.name;
        %     batch.Setup.l1covariates(1).values{subj}(1) = rp;
        batch.Setup.covariates.files{cov}{subj}{1} = fullfile(D,S(subj).name,'Aud/loc/param',rp);
    end
end

%% second lvl cov


if only_right
    effects = {vector_cov1};
    batch.Setup.subjects.effect_names    = {'age'};
    batch.Setup.subjects.group_names    = {'rh_synesthetes', 'rh_controls'};
else
    effects = {vector_cov1, vector_cov2};
    batch.Setup.subjects.effect_names    = {'age', 'handedness'};
    batch.Setup.subjects.group_names    = {'rh_synesthetes', 'rh_controls', 'lh_syn', 'lh_con'};
end
% effects = {vector_cov2};
% batch.Setup.subjects.effect_names    = {'handedness'};
% batch.Setup.subjects.group_names    = {'rh_synesthetes', 'rh_controls', 'lh_syn', 'lh_con'};

for neffect = 1 : length(batch.Setup.subjects.effect_names)
    batch.Setup.subjects.effects{neffect} = effects{neffect};
end

batch.Setup.subjects.groups         = group;
batch.Setup.covariates.add          = 1;

%% preprocessing part;
batch.Setup.done        = 1;
batch.Setup.overwrite   = 'Yes';
batch.Setup.secondarydatasets(1).functionals_label  = 'mni-space data';
batch.Setup.secondarydatasets(1).functionals_type   = 2;
batch.Setup.secondarydatasets(2).functionals_label  = 'smoothed data';
batch.Setup.secondarydatasets(2).functionals_type   = 4;
for nsub = 1:NSUBJECTS
    for nses=1:nsessions
        batch.Setup.secondarydatasets(1).functionals_explicit{nsub}{nses}{1} = FUNCTIONAL_FILE{nses,nsub};
        batch.Setup.secondarydatasets(2).functionals_explicit{nsub}{nses}{1} = FUNCTIONAL_FILE{nses,nsub};
    end
end

%% roi part;

batch.Setup.rois.names              = {'Grey Matter','White Matter','CSF','networks'};
batch.Setup.rois.multiplelabels     = [0,0,0,1];
batch.Setup.rois.regresscovariates  = [0,1,1,0];
for nses = 1 : nsessions
    for nsub = 1 : NSUBJECTS
        batch.Setup.rois.files{1,nses}{1, nsub}{1, 1}{1, 1}  = GM{nsub,nses};
        batch.Setup.rois.files{2,nses}{1, nsub}{1, 1}{1, 1}  = WM{nsub,nses};
        batch.Setup.rois.files{3,nses}{1, nsub}{1, 1}{1, 1}  = CSF{nsub,nses};
        batch.Setup.rois.files{4,nses}                       = '/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/conn18.b/rois/networks.nii';
        batch.Setup.rois.files{5,nses}{1, nsub}{1, 1}{1, 1} = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/rs/images/gcorr/lantmfg_bin_gcorr_roi_syn>con_510-3_510-2.nii';
        batch.Setup.rois.names{5}                            = 'lantmfg_rs';
        batch.Setup.rois.files{6,nses}{1, nsub}{1, 1}{1, 1} = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/rs/images/gcorr/sma_bin_gcorr_roi_syn>con_510-3_510-2.nii';
        batch.Setup.rois.names{6}                            = 'sma_rs';
        batch.Setup.rois.files{7,nses}{1, nsub}{1, 1}{1, 1} = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/rs/images/gcorr/lantmfg_gcorr_roi_and_cluster_syn>con_10-3_510-2.nii';
        batch.Setup.rois.names{7}                            = 'lantmfg_rs_and_cluster';
        batch.Setup.rois.files{8,nses}{1, nsub}{1, 1}{1, 1} = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/rs/images/gcorr/sma_gcorr_roi_and_cluster_syn>con_10-3_510-2.nii';
        batch.Setup.rois.names{8}                            =  'sma_rs_and_cluster';
        for tmp_roi = 1 : length(aud_coord_names)
            batch.Setup.rois.files{tmp_roi + 8, nses}{1, nsub}{1, 1}{1, 1}  = ROI_AUD_FILE{tmp_roi}{nsub};
            batch.Setup.rois.names{tmp_roi + 8}                             = sprintf('%s roi', aud_coord_names{tmp_roi});
        end
        
    end
    
    batch.Setup.rois.dimensions{1,nses} = 1;
    batch.Setup.rois.dimensions{2,nses} = 16;
    batch.Setup.rois.dimensions{3,nses} = 16;
end

batch.Setup.rois.add = 0; % if set to 1, it does the preproc for Roi

for nses = 1 : nsessions
    for nsub = 1 : NSUBJECTS
        load(ONSET_AUD_FILE{nsub,nses});
        for ncondition = 1 : size(names, 1)
            batch.Setup.conditions.names{ncondition} = names{ncondition};
            batch.Setup.conditions.onsets{ncondition}{nsub}{nses} = onsets{ncondition};
            batch.Setup.conditions.durations{ncondition}{nsub}{nses} = durations{ncondition};
        end
        batch.Setup.conditions.add = 0;
    end
end


batch.Setup.done = 1;
% uncomment the following 3 lines if you prefer to run one step at a time:
conn_batch(batch); % runs Preprocessing and Setup steps only
clear batch;
batch.filename               = fullfile(cwd, proj_name);            % Existing conn_*.mat experiment name

% ..........................................................................................

%% .................................... DENOISING step .....................................

% ..........................................................................................
% CONN Denoising                                    % Default options (uses White Matter+CSF+realignment+scrubbing+conditions as confound regressors); see conn_batch for additional options

if modif_freq
    batch.Denoising.filter      = [0.01, 0.1];                 % frequency filter (band-pass values, in Hz)
else
    batch.Denoising.filter      = [0.008, 0.09];                 % frequency filter (band-pass values, in Hz)
end
batch.Denoising.detrending  = 1;
batch.Denoising.regbp       = 1;
batch.Denoising.done        = 1;
batch.Denoising.overwrite   = 'Yes';

% uncomment the following 3 lines if you prefer to run one step at a time:
conn_batch(batch); % runs Denoising step only
clear batch;
batch.filename=fullfile(cwd, proj_name);            % Existing conn_*.mat experiment name

% ..........................................................................................

%% ............................... FIRST-LEVEL ANALYSIS step ...............................

% ..........................................................................................

% CONN Analysis                                     % Default options (uses all ROIs in conn/rois/ as connectivity sources); see conn_batch for additional options
batch.Analysis.done=1;
batch.Analysis.overwrite='Yes';
conn_batch(batch);% Run all analyses

% CONN Display
% launches conn gui to explore results
% conn
% conn('load',fullfile(cwd,'flc_all_controls.mat'));
% conn gui_results