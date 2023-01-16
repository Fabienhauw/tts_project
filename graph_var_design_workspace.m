filespath         = '/home/fabien.hauw/Desktop/toolbox_matlab/GraphVar_2.03a/workspaces/syn_con_rh_brainnetome/data/CorrMatrix';
cd(filespath);
files_mat = dir('*.mat');
brainSheet        = 'BN_Atlas_246_LUT_reoriented.csv';
corrVar           = '';
fieldName         = {'CorrMatrix'};
filename_end      = 4;
filename_start    = 12;

for i = 1 : length(files_mat)
    files{1,i}             = fullfile(filespath, files_mat(i).name);
    files_rel{i,1}         = fullfile(filespath, files_mat(i).name);
end

partVar           = '';
variableSheet     = 'Variables.csv';

save('/home/fabien.hauw/Desktop/toolbox_matlab/GraphVar_2.03a/workspaces/syn_con_rh_brainnetome/Workspace.mat', 'brainSheet', 'corrVar', 'fieldName', 'filename_end', 'filename_start', 'files', 'files_rel', 'partVar', 'variableSheet');