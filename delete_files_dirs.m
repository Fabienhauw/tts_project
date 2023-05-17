clear;clc;
addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))
D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];
a = 1; b = 48;

% delete files
for k = a:b
    which_dir = fullfile(D, S(k).name, 'Aud/loc/cpt_data');
    cd(which_dir)
    files = dir('*speech_scramble.mat');
    filenames = fullfile(which_dir, {files.name});
    try
        delete( filenames{:} );
    end
end

% for dir
for k = a:b
    which_dir = fullfile(D, S(k).name,'Aud/loc/stats_s5_without_resting/dcm_model_param_modul/lex_cond/all_4_rois_models_4mm');
    try
        cd(which_dir)
        files = dir('stru*9'); % never enter the exact name, or it will also erase repo '.' and '..', so the repo above...
        files = files([files.isdir]);
        filenames = fullfile(which_dir, {files.name});
        for tmp = 1 : length(filenames)
            try
                rmdir(filenames{tmp}, 's')
            end
        end
    end
end