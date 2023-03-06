%% Use MarsBaR to make spherical ROIs


%% Set general options

outDir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/dyslexia_feng_fh';
sphereRadius = 4; % mm
% ref_img = spm_vol('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/ANOVA_s5_without_resting_Control01_to_Sujet22_s8/syn_speech>baseline_10-3_510-2.nii');
ref_img = spm_vol('/home/fabien.hauw/Desktop/toolbox_matlab/MRIcron_linux/mricron/Resources/templates/ch2better.nii.gz');

% coordinates are nvoxels rows by 3 columns for X,Y,Z
load dyslex_meta_coord_mean;
% coords = ROIs_coord;
% names = ROIs_names;
%% Error checking: directory exists, MarsBaR is in path
if ~isdir(outDir)
    mkdir(outDir);
end

if ~exist('marsbar')
    error('MarsBaR is not installed or not in your matlab path.');
end


%% Make rois
fprintf('\n');

for i=1:size(coords,1) %needs a variable named coords, with coordinates in it
    thisCoord = coords(i,:);
    
    fprintf('Working on ROI %d/%d...', i, size(coords,1));
    
    roiLabel = sprintf('%g_%g_%g', thisCoord(1), thisCoord(2), thisCoord(3));
%     roiLabel = sprintf('%s_%d_%d_%d', ROIs_names{i}, thisCoord(1), thisCoord(2), thisCoord(3));

    sphereROI = maroi_sphere(struct('centre', thisCoord, 'radius', sphereRadius));
    
    outName = fullfile(outDir, sprintf('%dmmsphere_%s_roi', sphereRadius, roiLabel));    
%     outName = fullfile(outDir, sprintf('%dmmsphere_%s_roi', sphereRadius, roiLabel));

    outName = strrep(outName, '.', '_');
    % save MarsBaR ROI (.mat) file
%     saveroi(sphereROI, [outName '.mat']);
    
    % save the Nifti (.nii) file
%     save_as_image(sphereROI, [outName '.nii']);
    mars_rois2img(sphereROI, [outName '.nii'], ref_img);
    
    fprintf('done.\n');
    
end


fprintf('\nAll done. %d ROIs written to %s.'),
