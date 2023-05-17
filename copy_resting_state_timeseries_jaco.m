clear;clc;
addpath('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm');
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/VBA-toolbox'))

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

to_dir   = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/rs_for_jaco/TimeSeries';

cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/rs/syn_group_rs_rh_s5_atlas_conn/data')
ROIs = dir('ROI_Subject*');

for tmp_roi = 1 : length(ROIs)
    load(ROIs(tmp_roi).name)
    data = data(4:end);
    names = names(4:end);
    save(fullfile(to_dir, sprintf('TimeSeriesSubject%d', tmp_roi)), 'data', 'names')
end