tts_group.path_to_subject;
tts_group.subj_choice;

erase = input('Do you want to erase previous unframed color models? [yes/no] ', 's');

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
        
        which_dir=fullfile(D, S(k).name,'Vis/unfr_col/stats');
        dinfo = dir(which_dir);
        dinfo([dinfo.isdir]) = [];   %skip directories
        files = fullfile(which_dir, {dinfo.name});
        if ~exist('redo','var')
            redo=0;
        end
        if redo & ~isempty (files)
            delete (files{:});
        end
        
        filename = fullfile(D,S(k).name,'Vis/unfr_col/param');
        cd (filename);
        json=dir('*.json');
        json=json.name;
        
        res = get_string_from_json(json, {'RepetitionTime'}, {'num'});
        if res{1}>100
            TR          = res{1}/1000; % millisecond -> second
        else
            TR          = res{1}; % second
        end
        
        matlabbatch{i}.spm.stats.fmri_spec.dir = {fullfile(D, S(k).name,'Vis/unfr_col/stats')};
        matlabbatch{i}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{i}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        scans={};
        cd(fullfile(D, S(k).name,'Vis/unfr_col/swf'))
        vol_name = dir('s*wts_OC.nii');
        vol_name=fullfile(D, S(k).name,'Vis/unfr_col/swf', vol_name(1).name);
        nb_vol = size(spm_vol(vol_name),1);
        
        for v=1:nb_vol
            volume=sprintf('%s,%d',vol_name,v);
            scans=[scans;volume];
        end
        matlabbatch{i}.spm.stats.fmri_spec.sess.scans = scans;
        matlabbatch{i}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        multi_name=fullfile(D,S(k).name,'Vis/unfr_col/cpt_data/Timedata_');
        multi_name = sprintf('%s%s%s',multi_name,S(k).name,'_Unframed_Col.mat');
        matlabbatch{i}.spm.stats.fmri_spec.sess.multi = {multi_name};
        matlabbatch{i}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        multi_reg=fullfile(D,S(k).name,'Vis/unfr_col/param');
        cd (multi_reg);
        mr = 'multiple_regressors.txt';
        multi_reg=fullfile(multi_reg,mr);
        matlabbatch{i}.spm.stats.fmri_spec.sess.multi_reg = {multi_reg};
        matlabbatch{i}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{i}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{i}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{i}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{i}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{i}.spm.stats.fmri_spec.mthresh = 0;
        mask = fullfile(D,S(k).name,'anat/brain_extraction_mask.nii');
        matlabbatch{i}.spm.stats.fmri_spec.mask = {mask};%mask
        matlabbatch{i}.spm.stats.fmri_spec.cvi = 'AR(1)';
        i=i+1;        
        
        filename=fullfile(D, S(k).name,'Vis/unfr_col/stats/SPM.mat');
        matlabbatch{i}.spm.stats.fmri_est.spmmat = {filename};
        matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
        i=i+1;
        
        matlabbatch{i}.spm.stats.con.spmmat = {filename};
        color           = [1 0];
        uncolor         = [0 1];
        EOI             = [eye(2)];
        
        images =        color + uncolor;
        
        values = {...
            color, uncolor, images, color-uncolor, ...
            EOI, ...
            };
        
        names = {...
            'color', 'uncolor', 'images', 'color - uncolor', 'uncolor - color' ...
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

% cd '/home/fabien/Desktop/Fabien/NeoTopLex'