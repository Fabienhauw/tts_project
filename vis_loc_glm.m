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
erase = input('Do you want to erase previous visual models? [yes/no] ', 's');

redo = 1;
if isequal(erase,'yes')
    redo = 1;
elseif isequal(erase,'no')
    redo = 0;
end

for k = a:b
    if isdir (fullfile(D, S(k).name,'Vis'))
        %% modele specification:
        if exist ('i', 'var')==0
            i=1;
        end
        
        which_dir = fullfile(D, S(k).name,'Vis/loc/stats_s5');
        if ~isdir(which_dir)
            mkdir(which_dir)
        end
        
        dinfo = dir(which_dir);
        dinfo([dinfo.isdir]) = [];   %skip directories
        files = fullfile(which_dir, {dinfo.name});
        if ~exist('redo','var')
            redo=0;
        end
        if redo & ~isempty (files)
            delete (files{:});
        end
        
        filename = fullfile(D,S(k).name,'Vis/loc/param');
        cd (filename);
        json=dir('*.json');
        json=json.name;
        
        res = get_string_from_json(json, {'RepetitionTime'}, {'num'});
        if res{1}>100
            TR          = res{1}/1000; % millisecond -> second
        else
            TR          = res{1}; % second
        end
        
%         matlabbatch{i}.spm.stats.fmri_spec.dir = {fullfile(D, S(k).name,'Vis/loc/stats_s5')};
%         matlabbatch{i}.spm.stats.fmri_spec.timing.units = 'secs';
%         matlabbatch{i}.spm.stats.fmri_spec.timing.RT = TR;
%         matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t = 16;
%         matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
%         
%         scans={};
%         cd(fullfile(D, S(k).name,'Vis/loc/swf'))
%         vol_name = dir('s5*wts_OC.nii');
%         vol_name=fullfile(D, S(k).name,'Vis/loc/swf', vol_name(1).name);
%         nb_vol = size(spm_vol(vol_name),1);
%         
%         for v=1:nb_vol
%             volume=sprintf('%s,%d',vol_name,v);
%             scans=[scans;volume];
%         end
%         matlabbatch{i}.spm.stats.fmri_spec.sess.scans = scans;
%         matlabbatch{i}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
%         multi_name=fullfile(D,S(k).name,'Vis/loc/cpt_data/Timedata_');
%         multi_name = sprintf('%s%s_Vis.mat',multi_name,S(k).name);
%         matlabbatch{i}.spm.stats.fmri_spec.sess.multi = {multi_name};
%         matlabbatch{i}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
%         multi_reg=fullfile(D,S(k).name,'Vis/loc/param');
%         cd (multi_reg);
%         mr = 'multiple_regressors.txt';
%         multi_reg=fullfile(multi_reg,mr);
%         matlabbatch{i}.spm.stats.fmri_spec.sess.multi_reg = {multi_reg};
%         matlabbatch{i}.spm.stats.fmri_spec.sess.hpf = 128;
%         matlabbatch{i}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
%         matlabbatch{i}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
%         matlabbatch{i}.spm.stats.fmri_spec.volt = 1;
%         matlabbatch{i}.spm.stats.fmri_spec.global = 'None';
%         matlabbatch{i}.spm.stats.fmri_spec.mthresh = 0;
%         mask = fullfile(D,S(k).name,'anat/brain_extraction_mask.nii');
%         matlabbatch{i}.spm.stats.fmri_spec.mask = {mask};%mask
%         matlabbatch{i}.spm.stats.fmri_spec.cvi = 'AR(1)';
%         i=i+1;        
        
        filename=fullfile(D, S(k).name,'Vis/loc/stats_s5/SPM.mat');
%         matlabbatch{i}.spm.stats.fmri_est.spmmat = {filename};
%         matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
%         matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
%         i=i+1;
% %         
        matlabbatch{i}.spm.stats.con.spmmat = {filename};
        faces           = [1 0 0 0 0 0 0];
        houses          = [0 1 0 0 0 0 0];
        tools           = [0 0 1 0 0 0 0];
        numbers         = [0 0 0 1 0 0 0];
        words           = [0 0 0 0 1 0 0];
        odds            = [0 0 0 0 0 1 0];
        motor           = [0 0 0 0 0 0 1];
        EOI             = [eye(5)];
        
        images =        houses + faces + tools;
        
        values = {...
            faces, houses, tools, numbers, words, odds, motor, images, ...
            2*words - (faces+houses), 2*faces-(houses+words), 2*houses-(faces+words), ...
            3*words - (faces+houses+tools), 3*faces-(houses+words+tools), 3*houses-(faces+words+tools), 3*tools-(faces+houses+words),...
            3*numbers-(faces+houses+tools), 3*numbers-(faces+houses+words), numbers-words, words - numbers, 2*numbers-(faces+houses), ...
            EOI, ...
            };
        
        names = {...
            'faces', 'houses', 'tools', 'numbers', 'words', 'odds', 'motor', 'images', ...
            'words - (faces+houses)', 'faces-(houses+words)', 'houses-(faces+words)', ...
            'words - (faces+houses+tools)', 'faces-(houses+words+tools)', 'houses-(faces+words+tools)', 'tools -(faces+houses+words)', ...
            'numbers-(faces+houses+tools)', 'numbers-(faces+houses+words)', 'numbers-words', 'words - numbers', '2*numbers-(faces+houses)', ...
            'EOI',...
            };
        
        for j=1:length(values)-1
            matlabbatch{i}.spm.stats.con.consess{j}.tcon.name = names{j};
            matlabbatch{i}.spm.stats.con.consess{j}.tcon.weights = values{j};
            matlabbatch{i}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
        end
        j=length(values);
        matlabbatch{i}.spm.stats.con.consess{j}.fcon.name = names{j};
        matlabbatch{i}.spm.stats.con.consess{j}.fcon.weights = values{j};
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
par.jobname  = 'vis_glm';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);