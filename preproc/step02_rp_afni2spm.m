%% #matvol

clear
clc

load e

run_dir = e.getSerie('run').removeEmpty.getPath;

% input
afni_dir = gdir(run_dir,'afni');
dfile    = gfile(afni_dir,'dfile_rall.1D');

% output
tedana_dir = gdir(run_dir,'tedana');
output_dir = tedana_dir;

par.redo = 0;

job_rp_afni2spm(dfile, output_dir, par);