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
for k = a:b
    which_dir = fullfile(D, S(k).name, 'Aud/loc/stats_s5_without_resting');
    cd(which_dir)
    files = dir('*PPI*');
    filenames = fullfile(which_dir, {files.name});
    try
        delete( filenames{:} );
    end
    
    
    files = dir('*sphere*');
    filenames = fullfile(which_dir, {files.name});
    try
        delete( filenames{:} );
    end
    
    
    which_dir = fullfile(D, S(k).name);
    cd(which_dir)
    files = dir('*words_pw*');
    files = files([files.isdir]);
    filenames = fullfile(which_dir, {files.name});
    for tmp = 1 : length(filenames)
        try
            rmdir(filenames{tmp}, 's')
        end
    end
    
    files = dir('*sphere*');
    files = files([files.isdir]);
    filenames = fullfile(which_dir, {files.name});
    for tmp = 1 : length(filenames)
        try
            rmdir(filenames{tmp}, 's')
        end
    end

end