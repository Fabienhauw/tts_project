% to define a GLM for words/pseudowords run, with normal/reverse speech as a
% modulatory parameter 1/-1 for the DCM model.

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
erase = input('Do you want to erase previous auditive models? [yes/no] ', 's');

redo = 1;
if isequal(erase,'yes')
    redo = 1;
elseif isequal(erase,'no')
    redo = 0;
end

for k = a : b
    if isdir (fullfile(D, S(k).name,'Aud'))
        %% modele specification:
        if exist ('i', 'var')==0
            i=1;
        end
        which_dir=fullfile(D, S(k).name,'Aud/loc/stats_s5_without_resting/dcm_model_param_modul');
%         which_dir=fullfile(D, S(k).name,'Aud/loc/stats/dcm_model_param_modul_endroit');
        if ~isdir(which_dir)
            mkdir(which_dir)
        end
        
        dinfo = dir(which_dir);
        dinfo([dinfo.isdir]) = [];   %skip directories
        filenames = fullfile(which_dir, {dinfo.name});
        if ~exist('redo','var')
            redo=0;
        end
        if redo & ~isempty (filenames)
            delete( filenames{:} );
        end
        
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
        
%         matlabbatch{i}.spm.stats.fmri_spec.dir = {which_dir};
%         matlabbatch{i}.spm.stats.fmri_spec.timing.units = 'secs';
%         matlabbatch{i}.spm.stats.fmri_spec.timing.RT = TR;
%         matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t = 16;
%         matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
%         
%         clear scans;
%         scans={};
%         cd(fullfile(D, S(k).name,'Aud/loc/swf'))
%         vol_name = dir('s5wts_OC.nii');
%         vol_name=fullfile(D, S(k).name,'Aud/loc/swf', vol_name(1).name);
%         nb_vol = size(spm_vol(vol_name),1);
%         
%         for v=1:nb_vol
%             volume=sprintf('%s,%d',vol_name,v);
%             scans=[scans;volume];
%         end
%         
%         matlabbatch{i}.spm.stats.fmri_spec.sess.scans = scans;
%         cd(fullfile(D,S(k).name,'Aud/loc/cpt_data'))
%         multi_name = dir('onsets_dcm_param*.mat');
%         load(multi_name.name)
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.name = names;
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.onset = onsets{1};
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.duration = durations{1};
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.tmod = 0;
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.pmod(1).name = 'speech_scrambled';
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.pmod(1).param = parametric_modul{1};
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.pmod(1).poly = 1;
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.pmod(2).name = 'words_pw';
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.pmod(2).param = parametric_modul{2};
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.pmod(2).poly = 1;
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond.orth = 0; % orthogonalization, try with 1 or 0...
%         matlabbatch{i}.spm.stats.fmri_spec.sess.multi = {''};
%         matlabbatch{i}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
%         multi_reg=fullfile(D,S(k).name,'Aud/loc/param');
%         cd (multi_reg);
%         mr = 'multiple_regressors.txt';
%         multi_reg=fullfile(multi_reg,mr);
%         matlabbatch{i}.spm.stats.fmri_spec.sess.multi_reg = {multi_reg};
%         
%         matlabbatch{i}.spm.stats.fmri_spec.sess.hpf = 128;
%         matlabbatch{i}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
%         matlabbatch{i}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
%         matlabbatch{i}.spm.stats.fmri_spec.volt = 1;
%         matlabbatch{i}.spm.stats.fmri_spec.global = 'None';
%         matlabbatch{i}.spm.stats.fmri_spec.mthresh = 0;
%         anat=fullfile(D,S(k).name,'anat');
%         mask = fullfile(anat,'brain_extraction_mask.nii');
%         matlabbatch{i}.spm.stats.fmri_spec.mask = {mask};%mask
%         matlabbatch{i}.spm.stats.fmri_spec.cvi = 'AR(1)';
%         i=i+1;
        
        filename=fullfile(which_dir,'SPM.mat');
%         matlabbatch{i}.spm.stats.fmri_est.spmmat = {filename};
%         matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
%         matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
%         i=i+1;
%         
        matlabbatch{i}.spm.stats.con.spmmat = {filename};
        
        sounds      =   [1 0 0];
        speech_scr  =   [0 1 0];
        words_pw    =   [0 0 1];

%         sounds      =   [1 0 0 0];
%         speech_scr  =   [0 1 0 0];
%         lists       =   [0 0 1 0];
%         words_pw    =   [0 0 0 1];
        
        values = {...
            sounds, speech_scr, words_pw, ...
            (sounds + speech_scr)/2, (sounds - speech_scr)/2, ...
%             sounds.*words_pw, -sounds.*words_pw, ...
            };
        
        names = {...
            'sounds', 'speech_scr', 'words_pw', ...
            'speech', 'scrambled speech', ...
%             'words', 'pw', ...
            };

%         values = {...
%             sounds, lists, speech_scr, words_pw, ...
%             (sounds + speech_scr)/2, (sounds - speech_scr)/2, ...
%             (lists + words_pw)/2, (lists - words_pw)/2,
%             };
%         
%         names = {...
%             'sounds', 'lists', 'speech_scr', 'words_pw', ...
%             'speech', 'scrambled speech', ...
%             'words', 'pw' ...
%             };
        
        for j=1:length(values)
            matlabbatch{i}.spm.stats.con.consess{j}.tcon.name = names{j};
            matlabbatch{i}.spm.stats.con.consess{j}.tcon.weights = values{j};
            matlabbatch{i}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
        end
        j = j+1;
        matlabbatch{i}.spm.stats.con.consess{j}.fcon.name = 'Effects of interest';
        matlabbatch{i}.spm.stats.con.consess{j}.fcon.weights = eye(3);
%         matlabbatch{i}.spm.stats.con.consess{j}.fcon.weights = eye(4);
        matlabbatch{i}.spm.stats.con.consess{j}.fcon.sessrep = 'none';
        matlabbatch{i}.spm.stats.con.delete = 1;
        i=i+1;
    end
end
% spm_jobman('run', matlabbatch)

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