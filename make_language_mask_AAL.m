clear;clc;
matlabbatch{1}.spm.util.imcalc.input = {'/network/lustre/iss02/cohen/data/Fabien_official/atlas/AAL/AAL3v1_1mm.nii,1'};
matlabbatch{1}.spm.util.imcalc.output = 'left_language_mask_from_aal3';
matlabbatch{1}.spm.util.imcalc.outdir = {'/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/masks'};
matlabbatch{1}.spm.util.imcalc.expression = '1*(i1==1) + 2*(i1==7) + 3*(i1==9) + 4*(i1==15) + 5*(i1==59) + 6*(i1==65) + 7*(i1==67) + 8*(i1==69) + 9*(i1==85) + 10*(i1==89) + 11*(i1==93)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run', matlabbatch)

% then create a new mask with marsbar, based on this image, and with the
% image space to which you want to realign your mask...