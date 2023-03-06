clear;clc;
Finterest = 27; %index of the F-contrast (EOI) in the SPM.mat;
groupvar = [1, 2, 3, 4, 0, 0, 0, 0, 0]; %contrasts we want to study;
timeaxis = [0:2:15]; %time range and resolution of the resulting plot, in seconds
% these are from all speech conditions > rest
% all_coordinates = [-45 -51 -10; -48 -44 23; -50 -16 50; -50 6 53; -40 -41 46; -70 -28 3; -60 12 -7; 50 -24 16; 65 4 0]; % VWFA, SMG, pre_roland, MFG, IPS, ltemp1, ltemp2, rtemp1, rtemp2;

% these are from phonology contrast
% all_coordinates = [-48 -41 16; -48 -48 -14; -45 19 20; -50 -6 48; -58 -6 -4; -2 4 60]; % SMG, VWFA, IFG, MFG, mlSTG, SMA;

% all_coordinates = [-58 -44 23]; % mIFG;
modeldir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/Sujet01/Aud/loc/stats_s5_without_resting';

%% Single subject analysis

% stanplot_singlemodel(modeldir,Finterest,groupvar,timeaxis,coordinates)


%% Group analysis

for coordinates = 1:size(all_coordinates,1)
    stanplot_group(Finterest,groupvar,timeaxis, all_coordinates(coordinates,:))
end


%% Group analysis for multiple voxels;

for coordinates = 1 : size(all_coordinates,1)
    stanplot_group_multiple_vox(Finterest,groupvar,timeaxis, all_coordinates(coordinates,:))
end


%% Group analysis for best voxel;

for coordinates = 1 : size(all_coordinates,1)
    stanplot_group_best_vox(Finterest,groupvar,timeaxis, all_coordinates(coordinates,:))
end

