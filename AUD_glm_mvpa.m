% use to make 1/ GLMs design for gcss approach, splitting onsets, durations
% and names in 2 parts; 2/ GLMs design for MVPA analyses, splitting them in
% 8 parts;

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

nb_repeat       = 8;
ncondition1     = 8;
ncondition2     = 0;
ncondition      = ncondition1 + ncondition2;
redo = 0;

% for k = a:b
%     new_spm_1   = fullfile(D,S(k).name,'Aud/loc/mvpa10mm_rois');
%     
%     if exist(new_spm_1) == 0
%         mkdir(new_spm_1)
%     end
%     
%     dinfo = dir(new_spm_1);
%     dinfo([dinfo.isdir]) = [];   %skip directories
%     filenames = fullfile(new_spm_1, {dinfo.name});
%     
%     %% modele specification:
%     if exist ('i', 'var')==0
%         i=1;
%     end
%     % First import the original timedata.
%     time_load=fullfile(D,S(k).name,'Aud/loc/cpt_data/Timedata_');
%     time_load = sprintf('%s%s%s',time_load,S(k).name,'_Aud.mat');
%     load(time_load);
%     tmp_spm_name    = strsplit(time_load, '/'); spm_name{k} = sprintf('mvpa_%s',tmp_spm_name{end});
%     
%     % names
%     
%     names = repmat(names,nb_repeat,1);
%     
%     % onsets
%     
%     for j = 1 : ncondition % for each condition
%         n_interm        = rem(size(onsets{j},1), nb_repeat);
%         stim = 1;
%         for n = 1 : nb_repeat
%             if n<= n_interm
%                 nb_stim     = floor(size(onsets{j},1)/nb_repeat) + 1;
%             else
%                 nb_stim     = floor(size(onsets{j},1)/nb_repeat);
%             end
%             onsets2{(n-1) * ncondition + j, 1} 	= onsets{j}(stim : stim + nb_stim - 1);
%             stim = stim + nb_stim;
%         end
%     end
%     for j = 1 : ncondition*nb_repeat
%         onsets{j}   = onsets2{j};
%     end
%     clear onsets2
%     
%     % durations
%     
%     for j = 1 : ncondition % for each condition
%         n_interm        = rem(size(durations{j},1), nb_repeat);
%         stim = 1;
%         for n = 1 : nb_repeat
%             if n<= n_interm
%                 nb_stim     = floor(size(durations{j},1)/nb_repeat) + 1;
%             else
%                 nb_stim      = floor(size(durations{j},1)/nb_repeat);
%             end
%             durations2{(n-1) * ncondition + j, 1} 	= durations{j}(stim : stim + nb_stim - 1);
%             stim = stim + nb_stim;
%         end
%     end
%     for j = 1 : ncondition*nb_repeat
%         durations{j}   = durations2{j};
%     end
%     clear durations2
%     
%     new_spm     = fullfile(new_spm_1, spm_name{k});
%     save(new_spm,'names','onsets','durations');
% end

%% batch specification for mvpa
clear matlabbatch
i = 1;
for k = a : b
    new_spm_1   = fullfile(D,S(k).name,'Aud/loc/mvpa10mm_rois');
    filename = fullfile(D,S(k).name,'Aud/loc/param');
    cd (filename);
    json=dir('*.json');
    json=json.name;
    
    res = get_string_from_json(json, {'RepetitionTime'}, {'num'});
    if res{1}>100
        TR          = res{1}/1000; % millisecond -> second
    else
        TR          = res{1}; % second
    end
    
    dinfo = dir(new_spm_1);
    dinfo([dinfo.isdir]) = [];   %skip directories
    filenames = fullfile(new_spm_1, {dinfo.name});
    
    if redo & ~isempty (filenames)
        delete( filenames{:} );
    end
    
%     matlabbatch{i}.spm.stats.fmri_spec.dir = {new_spm_1};
%     matlabbatch{i}.spm.stats.fmri_spec.timing.units = 'secs';
%     matlabbatch{i}.spm.stats.fmri_spec.timing.RT = TR;
%     matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t = 16;
%     matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
%     
%     clear scans;
%     scans={};
%     cd(fullfile(D, S(k).name,'Aud/loc/wf'))
%     vol_name = dir('wts_OC.nii');
%     vol_name=fullfile(D, S(k).name,'Aud/loc/wf', vol_name(1).name);
%     nb_vol = size(spm_vol(vol_name),1);
%     
%     for v=1:nb_vol
%         volume=sprintf('%s,%d',vol_name,v);
%         scans=[scans;volume];
%     end
%     
%     matlabbatch{i}.spm.stats.fmri_spec.sess.scans       = scans;
%     matlabbatch{i}.spm.stats.fmri_spec.sess.cond        = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
%     multi_name=fullfile(D,S(k).name,'Aud/loc/mvpa/mvpa_Timedata_');
%     multi_name = sprintf('%s%s%s',multi_name,S(k).name,'_Aud.mat');
%     matlabbatch{i}.spm.stats.fmri_spec.sess.multi       = {multi_name};
%     matlabbatch{i}.spm.stats.fmri_spec.sess.regress     = struct('name', {}, 'val', {});
%     multi_reg = fullfile(D,S(k).name,'Aud/loc/param');
%     cd (multi_reg);
%     mr = 'multiple_regressors.txt';
%     multi_reg = fullfile(multi_reg,mr);
%     matlabbatch{i}.spm.stats.fmri_spec.sess.multi_reg   = {multi_reg};
%     matlabbatch{i}.spm.stats.fmri_spec.sess.hpf         = 128;
%     
%     matlabbatch{i}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
%     matlabbatch{i}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
%     matlabbatch{i}.spm.stats.fmri_spec.volt             = 1;
%     matlabbatch{i}.spm.stats.fmri_spec.global           = 'None';
%     matlabbatch{i}.spm.stats.fmri_spec.mthresh          = 0;
%     anat = fullfile(D,S(k).name,'anat');
%     mask = fullfile(anat,'brain_extraction_mask.nii');
%     matlabbatch{i}.spm.stats.fmri_spec.mask             = {mask};
%     matlabbatch{i}.spm.stats.fmri_spec.cvi              = 'AR(1)';
%     i=i+1;
    
    filename = fullfile(new_spm_1,'SPM.mat');
%     matlabbatch{i}.spm.stats.fmri_est.spmmat            = {filename};
%     matlabbatch{i}.spm.stats.fmri_est.write_residuals   = 0;
%     matlabbatch{i}.spm.stats.fmri_est.method.Classical  = 1;
%     i=i+1;
% 
    matlabbatch{i}.spm.stats.con.spmmat = {filename};
    words           = [repmat([1 0 0 0 0 0 0 0],1,nb_repeat)];
    pseudowords     = [repmat([0 1 0 0 0 0 0 0],1,nb_repeat)];
    numbers         = [repmat([0 0 1 0 0 0 0 0],1,nb_repeat)];
    normal_speech   = [repmat([0 0 0 1 0 0 0 0],1,nb_repeat)];
    scramble_speech = [repmat([0 0 0 0 1 0 0 0],1,nb_repeat)];
    odds            = [repmat([0 0 0 0 0 1 0 0],1,nb_repeat)];
    motor           = [repmat([0 0 0 0 0 0 1 0],1,nb_repeat)];
    resting         = [repmat([0 0 0 0 0 0 0 1],1,nb_repeat)];

    
    alphabetic      =    words + pseudowords;
    lexicality      =    words - pseudowords;
    phonology       =    normal_speech - scramble_speech;
    
    values = {...
        words - resting, pseudowords - resting, alphabetic, lexicality, -lexicality, numbers - resting, ...
        normal_speech - resting, scramble_speech - resting,...
        phonology
        };
    
    names = {...
        'words', 'pseudowords', 'alphabetic', 'lexicality', '-lexicality', 'numbers', 'normal_speech', 'scramble_speech',...
        'phonology', ...
        };
    
    for j=1:length(values)
        matlabbatch{i}.spm.stats.con.consess{j}.tcon.name = names{j};
        matlabbatch{i}.spm.stats.con.consess{j}.tcon.weights = values{j};
        matlabbatch{i}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
    end
    
    matlabbatch{i}.spm.stats.con.delete = 1;
    i=i+1;
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