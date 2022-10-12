clear; clc;
addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))
tts_group.path_to_subject;

wd = pwd;

% res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/syn_vs_con_rh_s5';
res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/syn_vs_con_s5';
Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet'))));
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];
% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17
% gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
gaucher_appar = {'Control02|Control04|Control07|Control17'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con_app = S;
S_con_app(mask_gauch_con) = [];
S = S_con_app;

%%
for k=1:length(S)
    anat_dir = fullfile(D, S(k).name, 'anat');
    cd(anat_dir)
    anat = dir('r_slicerwmv*'); anat = fullfile(anat_dir,anat.name);
    aud_dir = fullfile(D, S(k).name, 'Aud/loc/stats_s5');
    cd(aud_dir);
    func = fullfile(aud_dir,'spmT_0011.nii');
    slicer(...
        {anat,func},... % anatomy and activation volumes (must be same format, reslice if needed)
        'limits',{[],[2 100]},... % when a layer's limit is empty, limits will be adjusted automatically
        'minClusterSize',{0,0},....
        'p-map',{false, false},... % to replace p with 1-p in case of p-map
        'slices',[25],... % slice to show: 0 is the right hemisphere, 40 median sagital.
        'view','ax',... % ax sag cor
        'colormaps',{1,2},... % 1: grey 2: hot 3: cold, type colormaps for options
        'labels',{[],[]},... % no colorbar
        'colormode','black',... % background
        'fontsize',[10 25 5],... % [Title, ColorBarLabel, Coordinates]
        'show',false,... % do not show images
        'noMat',true,... % do not save details as .mat
        'output',sprintf('sujet_%02d',k),... % name of files
        'margins',[0 0 0 0],...
        'showCoordinates',false...
        )
    a       = dir(sprintf('slicer_sujet_*%d.png', k));
    b{k}    = fullfile(aud_dir, a.name);
end

figure(1)
aga = montage(b,'Size',[4 11]); % nb de lignes et de colonnes de la planche

% delete(b{:});

%%
for k=1:length(S)
    anat_dir = fullfile(D, S(k).name, 'anat');
    cd(anat_dir)
    anat = dir('r_slicerwmv*'); anat = fullfile(anat_dir,anat.name);
    aud_dir = fullfile(D, S(k).name, 'Aud/loc/stats_s5');
    cd(aud_dir);
    func = fullfile(aud_dir,'mask.nii');
    
    slicer(...
        {anat,func},... % anatomy and activation volumes (must be same format, reslice if needed)
        'limits',{[],[0 30]},... % when a layer's limit is empty, limits will be adjusted automatically
        'minClusterSize',{0, 0},....
        'p-map',{false, false},... % to replace p with 1-p in case of p-map
        'slices',[25],... % slice to show: 0 is the right hemisphere, 40 median sagital.
        'view','ax',... % ax sag cor
        'colormaps',{1,2},... % 1: grey 2: hot 3: cold, type colormaps for options
        'labels',{[],[]},... % no colorbar
        'colormode','black',... % background
        'fontsize',[10 25 5],... % [Title, ColorBarLabel, Coordinates]
        'show',false,... % do not show images
        'noMat',true,... % do not save details as .mat
        'output',sprintf('sujet_%02d',k),... % name of files
        'margins',[0 0 0 0],...
        'showCoordinates',false...
        )
    
    a       = dir(sprintf('slicer_sujet_*%d.png', k));
    b{k}    = fullfile(aud_dir, a.name);
end

figure(1)
aga = montage(b,'Size',[4 11]); % nb de lignes et de colonnes de la planche

% delete(b{:});