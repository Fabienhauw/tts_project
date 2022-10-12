addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))

wd = pwd;
res_dir_base = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/syn_vs_con_age_cov_s5_paired_tests';

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control07|Control17 
gaucher = {'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control'};
mask_gauch =  cellfun(@isempty,(regexp({S.name},gaucher)));
pairs = {
    'S*18', 'C*13';
    'S*07', 'C*23';
    'S*05', 'C*24';
    'S*13', 'C*03';
    'S*03', 'C*18';
    'S*04', 'C*08';
    'S*01', 'C*19';
    'S*11', 'C*25';
    'S*10', 'C*11';
    'S*06', 'C*01';
    'S*08', 'C*09';
    'S*19', 'C*14';
    'S*17', 'C*20';
    'S*12', 'C*21';
    'S*14', 'C*22';
    'S*21', 'C*10';
    'S*15', 'C*16';
    'S*09', 'C*15';
    'S*16', 'C*26';
    'S*20', 'C*12';
    'S*02', 'C*06';
    'S*22', 'C*05';
    };


cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images');
syn = dir('Sujet*'); con = dir('Control*');
all_subj = [syn;con];

vector_age = [
    25.1013699; 70.8219178; 23.6821918; 24.3342466; 21.109589; 31.7260274; 18.8438356; ...
    36.3424658; 49.7753425; 27.2767123; 26.3616438; 40.8876712; 22.8246575; 43.1013699; 44.309589; ...
    51.2246575; 40.7561644; 18.5835616; 39.939726; 58.9534247; 43.1753425; ... % end of synesthetes
    35.8219178; 23.2054795; 21.9726027; 31.7945205; 30.7589041; 70.3808219; 26.3589041; ... %start of controls
    22.660274; 42.0027397; 45.5589041; 19.6027397; 55.9041096; 19.4958904; 38.3643836; ...
    49.9178082; 46.0739726; 51.8630137; 25.1945205; 23.0547945; 41.4383562; 41.5972603; ...
    29.865753; 26.57534247; 24.67945205; 29.72328767; 58.97534247; ...
    ];


vector_hand = [
    0; 0; 0; 0; 1; 0; 1; 0; 0; 0; 1; 0; 0; 1; 0; 1; 0; 0; 0; 0; 0; ... % end of synesthetes
    zeros(21,1); 1; 1; 1; 1; 1; ... % end of controls
    ]; %0 = right, 1 = left;

cov1 = [];
cov2 = [];

for tmp_p = 1 : length(pairs)
    suj1 = dir(sprintf('%s', pairs{tmp_p,1}));
    suj2 = dir(sprintf('%s', pairs{tmp_p,2}));
    for j = 1 : size(all_subj,1)
        if strcmpi(all_subj(j).name, suj1.name)
            cov1 = [cov1; vector_age(j)];
            cov2 = [cov2; vector_hand(j)];
        end
        if strcmpi(all_subj(j).name, suj2.name)
            cov1 = [cov1; vector_age(j)];
            cov2 = [cov2; vector_hand(j)];
        end
    end
end

clear matlabbatch
i=1;

names = {...
            'words', 'pseudowords', 'numbers', 'normal_speech', 'scramble_speech',...
            'odds', 'motor',...
            'lexicality', '-lexicality','(words + normal_speech + numbers) - 3*pseudowords',...
            'phonology', '-phonology', 'numbers - (words + pseudowords)', ...
            'normal_speech - words', 'words + pseudowords + numbers + normal_speech', ...
            'words - normal_speech', 'words - scramble_speech', 'pseudowords - scramble_speech', '(words+pseudowords) -  scramble_speech', ...
            '(words + normal_speech) -  scramble_speech', '(normal_speech + pseudowords) -  scramble_speech', '(normal_speech + pseudowords + words) -  scramble_speech',...
            '(normal_speech + pseudowords + words + numbers) - 4*scramble_speech',...
            'numbers - words', 'words - numbers', ...
            'EOI', ...
            };

cd(fullfile(D, S(1).name,'Aud/loc/stats_s5'))
all_contrast = dir(sprintf('con*.nii'));
all_contrast = all_contrast(~cellfun(@isempty,(regexp({all_contrast.name}, 'con_\d+.nii'))));

for tmp_con = 1 : length(names)
    res_dir = fullfile(res_dir_base,names{tmp_con});
    if ~isdir(res_dir)
        mkdir(res_dir)
    end
    contrast = all_contrast(tmp_con).name;
    
%     for tmp_p = 1 : length(pairs)
%         cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images');
%         suj1 = dir(sprintf('%s', pairs{tmp_p,1}));
%         suj2 = dir(sprintf('%s', pairs{tmp_p,2}));
%         
%         matlabbatch{i}.spm.stats.factorial_design.dir = {res_dir};
%         
%         scan1 = fullfile(D, suj1.name,'Aud/loc/stats_s5/');
%         scan1 = sprintf('%s%s,1',scan1,contrast);
%         
%         scan2 = fullfile(D, suj2.name,'Aud/loc/stats_s5/');
%         scan2 = sprintf('%s%s,1',scan2,contrast);
%         
%         matlabbatch{i}.spm.stats.factorial_design.des.pt.pair(tmp_p).scans = {
%                                                                   scan1
%                                                                   scan2
%                                                                   };
%     end
%     matlabbatch{i}.spm.stats.factorial_design.des.pt.gmsca = 0;
%     matlabbatch{i}.spm.stats.factorial_design.des.pt.ancova = 0;
%     
%     matlabbatch{i}.spm.stats.factorial_design.cov(1).c = cov1;
%     matlabbatch{i}.spm.stats.factorial_design.cov(1).cname = 'age';
%     matlabbatch{i}.spm.stats.factorial_design.cov(1).iCFI = 1;
%     matlabbatch{i}.spm.stats.factorial_design.cov(1).iCC = 1;
%     
%     matlabbatch{i}.spm.stats.factorial_design.cov(2).c = cov2;
%     matlabbatch{i}.spm.stats.factorial_design.cov(2).cname = 'handedness';
%     matlabbatch{i}.spm.stats.factorial_design.cov(2).iCFI = 1;
%     matlabbatch{i}.spm.stats.factorial_design.cov(2).iCC = 1;
% 
%     matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
%     matlabbatch{i}.spm.stats.factorial_design.masking.tm.tm_none = 1;
%     matlabbatch{i}.spm.stats.factorial_design.masking.im = 1;
%     exp_mask = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks',all_subj(1).name);
%     exp_mask = sprintf('%s_%s_aud_loc_mask_thr_s5.nii',exp_mask, all_subj(end).name);
%     matlabbatch{i}.spm.stats.factorial_design.masking.em = {exp_mask};
%     matlabbatch{i}.spm.stats.factorial_design.globalc.g_omit = 1;
%     matlabbatch{i}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
%     matlabbatch{i}.spm.stats.factorial_design.globalm.glonorm = 1;
%     
%     i = i+1;

%     matlabbatch{i}.spm.stats.fmri_est.spmmat = {fullfile(res_dir, 'SPM.mat')};
%     matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
%     matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
% 
%     i = i+1;

    matlabbatch{i}.spm.stats.con.spmmat = {fullfile(res_dir, 'SPM.mat')};
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'Synesthetes';
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.name = 'Controls';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.weights = [0 1];
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.name = '-Synesthetes';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.weights = -1;
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.name = '-Controls';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{5}.tcon.name = 'syn-con';
    matlabbatch{i}.spm.stats.con.consess{5}.tcon.weights = [1 -1];
    matlabbatch{i}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{6}.tcon.name = 'syn+con';
    matlabbatch{i}.spm.stats.con.consess{6}.tcon.weights = [1 1];
    matlabbatch{i}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{7}.tcon.name = 'con-syn';
    matlabbatch{i}.spm.stats.con.consess{7}.tcon.weights = [-1 1];
    matlabbatch{i}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

    i = i+1;
end


cd '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts'

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'aud_glm';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);