clear;

%% parameters:
only_right_hand = 1;
volume_option = 2; %1 = ROI, 2 = coordinates;
sphereradius = 10; %% in  mm a 10 mm sphere focuses tightly on the VWFA locationend, useful only for volume option = 2;
voxel_option = 1; % 1 = keep n best voxel; 2 = keep n voxel above a threshold (threshold is applied to the tested contrast, not the reference contrast);
per_vox = 10; % percentage of voxels you want to keep
pvalue = 0.001;
flipsign = 1; %% default is to keep the contrast "as is" with its sign

%%
% D = '/media/biggest_drive/SYNESTHEX/final_images';
D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
S = dir(D);
mask = ismember({S.name}, {'.', '..'});
S(mask) = [];

if only_right_hand
    gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
    mask_gauch =  ~cellfun(@isempty,(regexp({S.name},{'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control'})));
    res_dir = '/home/fabien.hauw/Desktop/Fabien/NeoTopLex/+tts_group/+second_level/RH_';
else
    gaucher_appar = {'Control02|Control04|Control07|Control17|Sujet'};
    mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
    res_dir = '/home/fabien.hauw/Desktop/Fabien/NeoTopLex/+tts_group/+second_level/';
end
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con_app = S;
S_con_app(mask_gauch_con) = [];

S_droit = S;
S_droit(mask_gauch) = [];

S_effect = [S_droit ; S_con_app];
S = S_effect;


%% ROI definition
ROIfiles = {};
ROIfiles = {
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/svc_roi/sphere_15--3_-1_68.nii',...
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/svc_roi/sphere_15--3_-55_59.nii',...
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/svc_roi/sphere_15--33_5_26.nii',...
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/svc_roi/sphere_15--48_-43_32.nii',...
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/svc_roi/sphere_8--44_-50_-14.nii',...
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/svc_roi/sphere_10--42_-46_-10.nii',...
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/svc_roi/sphere_15--42_-46_-10.nii',...
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/masks/smg_mask_from_brainnetome_and_flc_activations_mars.nii',...
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/masks/ifg_mask_from_brainnetome_and_flc_activations_mars.nii',...
    '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/masks/vwfa_mask_from_brainnetome_and_flc_activations_mars.nii',...
    };
%
ROInames = {'SMA','PRECUN','IFG', 'SMG', 'VWFA', 'VWFA2', 'VWFA3', 'SMG_atlas', 'IFG_atlas', 'VWFA_atlas'};

if volume_option == 1
    all_xyzmm = ROInames;
elseif volume_option == 2
    
    % different coordinates, according to where they come from
%     AUD normal > scr, RH synesthetes + controls,  
%                  SMG,      pSTS,        aSTS,    vis VWFA,      VWFA,         aVWFA,       iTPol,       sTPol,      sIFG,      mIFG,      iIFG.
%     all_xyzmm = [-58 -44 23; -52 -38 0; -58 -6 -2; -48 -54 -20; -52 -51 -20; -40 -36 -27; -40 -11 -47; -40 -14 -34; -50 -8 50; -45 22 23; -50 26 -2]; 
% --> no diff    

%     AUD "big" phonology, RH synesthetes + controls
%                    SMG,      pSTS,        aSTS,    vis VWFA,      VWFA,         aVWFA,       iTPol,       sTPol,      sIFG,      mIFG,      iIFG.
%         all_xyzmm = [-48 -41 18; -52 -48 6; -58 -6 -4; -48 -54 -20; -42 -51 -14; -40 -38 -24; -38 -8 -44; -38 -18 -32; -50 -8 43; -42 12 26; -45 29 -2];
% --> p=0.03 in SMG

%     AUD words + PW, RH synesthetes + controls
%                    rSMG,      rpMTS,        aSTS,    vis VWFA,      VWFA,         aVWFA,       iTPol,       sTPol,      sIFG,      mIFG,      iIFG.
%         all_xyzmm = [48 -44 28; 60 -51 -7];
% --> no diff   

%     AUD words - PW, RH synesthetes + controls
%                    lvisuomot,    lSMA,      VWFA,    rvisuomot.
        all_xyzmm = [-28 -68 36; -30 2 46; -40 -34 -22; 28 -66 30];
% --> p=0.0439 in lvisuomot.

%     MVPA phonology, RH synesthetes + controls
%                   left aSTG,  right mSTG,   SMA,   l visual motor.
%         all_xyzmm = [-58 -6 -4; -52 -48 6; -2 22 66; -35 -68 48];
% --> no diff    

end

%% paths and scans

for k = 1:length(S)
    if ~isempty(regexp(S(k).name, 'Sujet'))
        group_S(k)    = 1; % 1 is for synesthetes, 2 for controls
    elseif ~isempty(regexp(S(k).name, 'Control'))
        group_S(k)    = 2; % 1 is for synesthetes, 2 for controls
    end
end

path_to_stats = {};
subjs={};
for k=1:numel(S)
    subjs                   = [subjs;S(k).name];
    path_to_subj            = fullfile(D, S(k).name, 'Aud/loc/mvpa10mm');
    path_to_stats           = [path_to_stats;path_to_subj];
end

cd(path_to_stats{k})
test_con = dir('results*');

for tmp_comp = 1 : length(test_con)
    str_path = strsplit(test_con(tmp_comp).name, '_');
    all_comp_inv{tmp_comp} = strjoin({str_path{1}, str_path{4}, str_path{3}, str_path{2}}, '_');
    all_comp{tmp_comp} = test_con(tmp_comp).name;
end

experiments = struct(...
    'test_path',    path_to_stats,...  % subpath to the test data SPM.mat inside each subject
    'data',         subjs...
    );


res_dir = sprintf('%scomparison_best_voxels/', res_dir);

if volume_option == 1
    res_dir = fullfile(res_dir, 'ROI');
elseif volume_option == 2
    res_dir = fullfile(res_dir, 'coord');
end

res_dir = fullfile(res_dir, 'decoding_10mm_aud');

if ~isdir(res_dir)
    mkdir(res_dir)
end

for zz=1:length(test_con)
    current_con = test_con(zz).name;
    for qq=1:size(all_xyzmm,1)
        
        if volume_option == 1
            xyzmm = sprintf('%s',ROInames{qq});
        elseif volume_option == 2
            xyzmm = all_xyzmm(qq,:);
        end
        
        %% start of computation of the results
        totsub = size(experiments,1);
        spmfiles={};
        spmfiles_test={};
        
        for nsub=1:totsub
            spmfiles_test{nsub}     = fullfile(experiments(nsub).test_path,'SPM.mat');
        end
        
        %% define the search volume
        if (volume_option == 1)
            load(spmfiles_test{1}); %%% load one SPM.mat just to get the transformation matrix iM
            ROIheader   = spm_vol(ROIfiles{qq});
            ROIvol      = spm_read_vols(ROIheader);
            ROInumbers  = unique(round(ROIvol(ROIvol>0)));
        elseif (volume_option == 2)
            if size(xyzmm,1)==1
                xyzmm = xyzmm';
            end
            load(spmfiles_test{1}); %%% load one SPM.mat just to get the transformation matrix iM
            iM = SPM.xVol.iM ;
            xyzvox = iM( 1:3, : ) * [ xyzmm ; 1 ] ;
            spherecoords = xyzmm';
            
            %% the following code is just to get the voxel sizes
            if ~isdir(fullfile(experiments(nsub).test_path, all_comp_inv{zz}))
                select_tfile = fullfile(experiments(nsub).test_path, all_comp{zz}, 'res_accuracy_minus_chance.nii');
%                 select_tfile = fullfile(experiments(nsub).test_path, all_comp{zz}, 'perm/res_accuracy_minus_chance_1-p_value.nii');
            else
                select_tfile = fullfile(experiments(nsub).test_path, all_comp_inv{zz}, 'res_accuracy_minus_chance.nii');
%                 select_tfile = fullfile(experiments(nsub).test_path, all_comp_inv{zz}, 'perm/res_accuracy_minus_chance_1-p_value.nii');
            end
            select_theader = spm_vol(select_tfile);
            select_tvol=spm_read_vols(select_theader);
            
            ROInames{1} =  'sphere';
            %%%% define sphere
            disp('computing sphere voxels....');
            vox = 1;
            for x=1:size(select_tvol,1)
                for y=1:size(select_tvol,2)
                    for z=1:size(select_tvol,3)
                        voxx = x * select_theader.mat(1,1) + select_theader.mat(1,4);
                        voxy = y * select_theader.mat(2,2) + select_theader.mat(2,4);
                        voxz = z * select_theader.mat(3,3) + select_theader.mat(3,4);
                        if (norm([voxx voxy voxz]-spherecoords)<=sphereradius)
                            spherevol(x,y,z)=1;
                        else
                            spherevol(x,y,z)=0;
                        end
                    end
                end
            end
        end
        
        fprintf( 'FIND BEST VOXELS %d/%d\n\n', qq + (size(all_xyzmm,1) * (zz - 1)), size(all_xyzmm,1) * length(test_con)) ;
        h = waitbar( 0, 'Finding best voxels...' ) ;
        
        for i_subj = 1:totsub  %%%% loop across subjects
            
            %%% read all the test images
            clear test_convol;
            
            if ~isdir(fullfile(experiments(i_subj).test_path, all_comp_inv{zz}))
                test_confile = fullfile(experiments(i_subj).test_path, all_comp{zz}, 'res_accuracy_minus_chance.nii');
%                 test_confile = fullfile(experiments(i_subj).test_path, all_comp{zz}, 'perm/res_accuracy_minus_chance_1-p_value.nii');
            else
                test_confile = fullfile(experiments(i_subj).test_path, all_comp_inv{zz}, 'res_accuracy_minus_chance.nii');
%                 test_confile = fullfile(experiments(i_subj).test_path, all_comp_inv{zz}, 'perm/res_accuracy_minus_chance_1-p_value.nii');
            end
            test_conheader      = spm_vol(test_confile);
            test_convol         = spm_read_vols(test_conheader);
            
            
            waitbar( (qq+(i_subj-1)*length(xyzmm))/(totsub*length(xyzmm)), h ) ;
            if (volume_option == 1)
                anal{qq}.volume_str = sprintf('Region %d (%s)', qq, ROInames{qq});
                %                     searchvol = double(round(ROIvol)==ROInumbers(qq));
                searchvol = double(round(ROIvol)==1);
            elseif (volume_option == 2)
                anal{qq}.volume_str = sprintf('%d mm sphere centered at [%d %d %d]',sphereradius,xyzmm);
                searchvol = spherevol;
            end
            
            
            
            test_convol(searchvol(:)==0) = Inf; % every voxel outside of the mask is set to 'Inf'
            tvol = flipsign * test_convol(:) .* (searchvol(:)>0);
            
            if (voxel_option == 1)
                [tvalues,xyz] = sort(tvol,'descend');
                xyz(isnan(tvalues)) = [];
                tvalues(isnan(tvalues)) = [];
                nvox = round(length(find(searchvol>0))*(per_vox/100));
                xyz = xyz(1:nvox);
                
                % then
                [tvalues_thr,xyz_thr] = sort(tvol,'descend');
                xyz_thr(isnan(tvalues_thr)) = [];
                tvalues_thr(isnan(tvalues_thr)) = [];
                load(spmfiles_test{i_subj});
                criticalt = tinv(1-pvalue,SPM.xX.erdf);
                xyz_thr(tvalues_thr<criticalt) = [];
                tvalues_thr(tvalues_thr<criticalt) = [];
                nvox_thr = length(find(xyz_thr>0));
            elseif (voxel_option == 2)
                load(spmfiles_test{i_subj});
                criticalt = tinv(1-pvalue,SPM.xX.erdf);
                [tvalues,xyz] = sort(tvol,'descend'); % comment this line and next 4 if you don't want the values to be sorted
                xyz(isnan(tvalues)) = [];
                tvalues(isnan(tvalues)) = [];
                xyz(tvalues<criticalt) = [];
                tvalues(tvalues<criticalt) = [];
%                 xyz = find(tvol>=criticalt); % uncomment if you don't want the values to be sorted
                nvox = length(find(xyz>0));
                nvox_thr = length(find(xyz>0));
            end
            
            anal{qq}.nbvoxels(i_subj) = nvox;
            anal{qq}.nbvoxels_thr(i_subj) = nvox_thr;
            if nvox>0
                [x,y,z]=ind2sub(size(test_convol),xyz);
                anal{qq}.voxels{i_subj}             = [ x y z ]';  %%% store the voxels for this subject
                anal{qq}.xyz{i_subj}                = xyz;  %%% this is the easiest storage: direct indexes to the data matrices, can be used as a data selector
                anal{qq}.coordsmm(1, i_subj)        = mean(x * select_theader.mat(1,1) + select_theader.mat(1,4)); %%% mean coordinates of the identified voxels
                anal{qq}.coordsmm(2, i_subj)        = mean(y * select_theader.mat(2,2) + select_theader.mat(2,4));
                anal{qq}.coordsmm(3, i_subj)        = mean(z * select_theader.mat(3,3) + select_theader.mat(3,4));
                anal{qq}.bvcoordsmm(1, i_subj)      = (x(1) * select_theader.mat(1,1) + select_theader.mat(1,4)); %%% mean coordinates of the identified voxels
                anal{qq}.bvcoordsmm(2, i_subj)      = (y(1) * select_theader.mat(2,2) + select_theader.mat(2,4));
                anal{qq}.bvcoordsmm(3, i_subj)      = (z(1) * select_theader.mat(3,3) + select_theader.mat(3,4));
                anal{qq}.coords(1, i_subj)          = mean(x);
                anal{qq}.coords(2, i_subj)          = mean(y);
                anal{qq}.coords(3, i_subj)          = mean(z);
                indiv_coords(i_subj,:)              = anal{qq}.coordsmm(:, i_subj);
                anal{qq}.all_activation{i_subj}     = test_convol(xyz);
                anal{qq}.activation(i_subj)         = mean(test_convol(xyz));
                anal{qq}.group(i_subj)              = group_S(i_subj);
            end
            
            if voxel_option == 1
                anal{qq}.voxel_str = sprintf('The %d best voxels for aud localizer %s.', nvox, current_con);
            elseif voxel_option == 2
                anal{qq}.voxel_str = sprintf('All voxels for aud localizer %s with p-value < %7.5f.', current_con, pvalue);
            end
            
            % to write a .nii map with only voxels kepts:
            test_convol_copy = test_convol;
            test_convol_copy(:) = 0;
            test_convol_copy(xyz) = 1;
            if volume_option==1
                test_conheader.fname = sprintf('%s_best_vox_in %s.nii',test_conheader.fname(1:end-4), xyzmm);
                test_conheader.descrip = sprintf('%s thresholded to 10 prc best voxels in region %s.',test_conheader.descrip, xyzmm);
            elseif volume_option ==2
                test_conheader.fname = sprintf('%s_best_vox_sph_%d_%d_%d.nii', test_conheader.fname(1:end-4),xyzmm);
                test_conheader.descrip = sprintf('%s thresholded to 10 prc best voxels in 10 mm rad sphere centered on %d %d %d',test_conheader.descrip, xyzmm);
            end
            
            spm_write_vol(test_conheader,test_convol_copy);
            
            
        end %%% subject loop
        close(h)
        
    end % coord loop
    if volume_option == 1
        save(sprintf('%s/best_voxels_in_aud_%s.mat',...
            res_dir, current_con), 'anal');
    elseif volume_option == 2
        save(sprintf('%s/best_voxels_in_aud_%s.mat',...
            res_dir, current_con), 'anal');
    end
    
end % contrast loop

% return

%% now the part to plot
cd(res_dir)

do_plot = 0;

results = dir('best*aud*');
count_fig = 0;

load('/home/fabien.hauw/Desktop/Fabien/NeoTopLex/+tts_group/cpt_data/pca_score.mat')

for tmp_anal = 1 : size(results, 1)
    clear anal
    load(results(tmp_anal).name)
    clear matcon Group_Id Subj_Id vox_activ
    for j = 1 : size(anal,2) % j = ROI;
        count = 0;
        for k = 1 : size(anal{j}.all_activation,2) % k = subjs;
            for l = 1 : size(anal{j}.all_activation{k},1)
                count = count + 1;
                matcon{j}{count, 1}     = anal{j}.group(k);
                Group_Id{j}{count, 1}   = anal{j}.group(k);
                matcon{j}{count, 2}     = k;
                Subj_Id{j}{count, 1}    = k;
                matcon{j}{count, 3}     = anal{j}.all_activation{k}(l);
                vox_activ{j}{count, 1}  = anal{j}.all_activation{k}(l);
            end
        end
    end
    
    
%     matcon_sma      = table(Subj_Id{1}, Group_Id{1}, vox_activ{1}); matcon_sma.Properties.VariableNames = {'Subj_Id', 'Group_Id', 'vox_activ'};
%     matcon_precun   = table(Subj_Id{2}, Group_Id{2}, vox_activ{2}); matcon_precun.Properties.VariableNames = {'Subj_Id', 'Group_Id', 'vox_activ'};
%     matcon_ifg      = table(Subj_Id{3}, Group_Id{3}, vox_activ{3}); matcon_ifg.Properties.VariableNames = {'Subj_Id', 'Group_Id', 'vox_activ'};
%     matcon_smg      = table(Subj_Id{4}, Group_Id{4}, vox_activ{4}); matcon_smg.Properties.VariableNames = {'Subj_Id', 'Group_Id', 'vox_activ'};
%     matcon_vwfa     = table(Subj_Id{5}, Group_Id{5}, vox_activ{5}); matcon_vwfa.Properties.VariableNames = {'Subj_Id', 'Group_Id', 'vox_activ'};
%     matcon_vwfa2    = table(Subj_Id{6}, Group_Id{6}, vox_activ{6}); matcon_vwfa2.Properties.VariableNames = {'Subj_Id', 'Group_Id', 'vox_activ'};
%     matcon_vwfa3    = table(Subj_Id{7}, Group_Id{7}, vox_activ{7}); matcon_vwfa3.Properties.VariableNames = {'Subj_Id', 'Group_Id', 'vox_activ'};
%     
%     writetable(matcon_sma, sprintf('%s/detail_best_voxels_in_sma_aud_con%04d.txt',...
%         res_dir, current_con));
%     writetable(matcon_precun, sprintf('%s/detail_best_voxels_in_precun_aud_con%04d.txt',...
%         res_dir, current_con));
%     writetable(matcon_ifg, sprintf('%s/detail_best_voxels_in_ifg_aud_con%04d.txt',...
%         res_dir, current_con));
%     writetable(matcon_smg, sprintf('%s/detail_best_voxels_in_smg_aud_con%04d.txt',...
%         res_dir, current_con));
%     writetable(matcon_vwfa, sprintf('%s/detail_best_voxels_in_vwfa_aud_con%04d.txt',...
%         res_dir, current_con));
%     writetable(matcon_vwfa2, sprintf('%s/detail_best_voxels_in_vwfa2_aud_con%04d.txt',...
%         res_dir, current_con));
%     writetable(matcon_vwfa3, sprintf('%s/detail_best_voxels_in_vwfa3_aud_con%04d.txt',...
%         res_dir, current_con));
    
    
    clear h
    for s = 1 : size(anal,2)
        tmp_mask = cellfun(@isempty,(anal{1,s}.all_activation(:)));
%         anal{1,s}.activation(tmp_mask) = '';
        anal{1,s}.coordsmm(:,tmp_mask) = '';
        anal{1,s}.bvcoordsmm(:,tmp_mask) = '';
        anal{1,s}.group_mask = anal{1,s}.group(~tmp_mask);
        
        % this is to uncomment to nmake a linear regression based on PCA
        % score for each subject, and make comparisons with the residuals
%         test = fitlm([anal{1,s}.activation(anal{1,s}.group_mask == 1), anal{1,s}.activation(anal{1,s}.group_mask == 2)], score_pca)
%         test = test.Residuals.Raw;
%         
%         [h{tmp_anal,s},p{tmp_anal,s},ci{tmp_anal,s},stats{tmp_anal,s}] = ...
%             ttest2(...
%             test(1:length(anal{1,s}.activation(anal{1,s}.group_mask == 1))),...
%         test(length(anal{1,s}.activation(anal{1,s}.group_mask == 1))+1:end));
        
        [h{tmp_anal,s},p{tmp_anal,s},ci{tmp_anal,s},stats{tmp_anal,s}] = ...
            ttest2(...
            anal{1,s}.activation(anal{1,s}.group_mask == 1), ...
            anal{1,s}.activation(anal{1,s}.group_mask == 2));
        tval{tmp_anal,s} = stats{tmp_anal,s}.tstat;
        
        [h_coord1{tmp_anal,s},p_coord1{tmp_anal,s},ci_coord1{tmp_anal,s},stats_coord1{tmp_anal,s}] = ...
            ttest2(...
            anal{1,s}.coordsmm(1,anal{1,s}.group_mask == 1), ...
            anal{1,s}.coordsmm(1,anal{1,s}.group_mask == 2));
        tval_coord1{tmp_anal,s} = stats_coord1{tmp_anal,s}.tstat;
        
        [h_coord2{tmp_anal,s},p_coord2{tmp_anal,s},ci_coord2{tmp_anal,s},stats_coord2{tmp_anal,s}] = ...
            ttest2(...
            anal{1,s}.coordsmm(2,anal{1,s}.group_mask == 1), ...
            anal{1,s}.coordsmm(2,anal{1,s}.group_mask == 2));
        tval_coord2{tmp_anal,s} = stats_coord2{tmp_anal,s}.tstat;
        
        [h_coord3{tmp_anal,s},p_coord3{tmp_anal,s},ci_coord3{tmp_anal,s},stats_coord3{tmp_anal,s}] = ...
            ttest2(...
            anal{1,s}.coordsmm(3,anal{1,s}.group_mask == 1), ...
            anal{1,s}.coordsmm(3,anal{1,s}.group_mask == 2));
        tval_coord3{tmp_anal,s} = stats_coord3{tmp_anal,s}.tstat;
        
        [h_coord_bv_1{tmp_anal,s},p_coord_bv_1{tmp_anal,s},ci_coord_bv_1{tmp_anal,s},stats_coord_bv_1{tmp_anal,s}] = ...
            ttest2(...
            anal{1,s}.bvcoordsmm(1,anal{1,s}.group_mask == 1), ...
            anal{1,s}.bvcoordsmm(1,anal{1,s}.group_mask == 2));
        tval_coord_bv_1{tmp_anal,s} = stats_coord_bv_1{tmp_anal,s}.tstat;
        
        [h_coord_bv_2{tmp_anal,s},p_coord_bv_2{tmp_anal,s},ci_coord_bv_2{tmp_anal,s},stats_coord_bv_2{tmp_anal,s}] = ...
            ttest2(...
            anal{1,s}.bvcoordsmm(2,anal{1,s}.group_mask == 1), ...
            anal{1,s}.bvcoordsmm(2,anal{1,s}.group_mask == 2));
        tval_coord_bv_2{tmp_anal,s} = stats_coord_bv_2{tmp_anal,s}.tstat;
        
        [h_coord_bv_3{tmp_anal,s},p_coord_bv_3{tmp_anal,s},ci_coord_bv_3{tmp_anal,s},stats_coord_bv_3{tmp_anal,s}] = ...
            ttest2(...
            anal{1,s}.bvcoordsmm(3,anal{1,s}.group_mask == 1), ...
            anal{1,s}.bvcoordsmm(3,anal{1,s}.group_mask == 2));
        tval_coord_bv_3{tmp_anal,s} = stats_coord_bv_3{tmp_anal,s}.tstat;
    
        [h_nvox{tmp_anal,s},p_nvox{tmp_anal,s},ci_nvox{tmp_anal,s},stats_nvox{tmp_anal,s}] = ...
            ttest2(...
            anal{1,s}.nbvoxels(anal{1,s}.group == 1), ...
            anal{1,s}.nbvoxels(anal{1,s}.group == 2));
        tval_nvox{tmp_anal,s} = stats_nvox{tmp_anal,s}.tstat;
        
        [h_nvox_thr{tmp_anal,s},p_nvox_thr{tmp_anal,s},ci_nvox_thr{tmp_anal,s},stats_nvox_thr{tmp_anal,s}] = ...
            ttest2(...
            anal{1,s}.nbvoxels_thr(anal{1,s}.group == 1), ...
            anal{1,s}.nbvoxels_thr(anal{1,s}.group == 2));
        tval_nvox_thr{tmp_anal,s} = stats_nvox_thr{tmp_anal,s}.tstat;
    end
    
    [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
    nbgraph = size(anal,2);
    
    if do_plot
        for tmp_grap = 1 : nbgraph
            count_fig = count_fig + 1;
            figure(count_fig)
            tmp_matcon = matcon{tmp_grap};
            group = [tmp_matcon{:,1}]';
            vox_activ = [tmp_matcon{:,3}]';
            syn_group = group == 1;
            con_group = group == 2;
            vox_activ_syn = vox_activ(syn_group);
            vox_activ_con = vox_activ(con_group);
            
            h1 = raincloud_plot(vox_activ_syn, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35, 'box_col_match', 1);
            h2 = raincloud_plot(vox_activ_con, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75, 'box_col_match', 1);
            % red is the control group;
            
            title(sprintf('plot of the %d prc best voxels in %s', per_vox, anal{tmp_grap}.volume_str));
        end
    end
end

fprintf('\np-values for:              SMG,     pSTS,  aSTS,  vis VWFA, VWFA,   aVWFA,  iTPol,  sTPol,  sIFG,  mIFG,  iIFG.')
% fprintf('\np-values for:              left aSTG,  right mSTG,   SMA,   l visual motor.')
for tmp_stat = 1 : size(results, 1)
    tmp_con = regexp(results(tmp_stat).name, '[^_]+_vs_+[^_]', 'match');
    fprintf('\np-values for con %s: ', tmp_con{1})
    fprintf('%4.4f  ', p{tmp_stat,:});
    fprintf('\nt-values for con %s: ', tmp_con{1})
    fprintf('%4.4f  ', tval{tmp_stat,:});fprintf('\n')
end
p_coord1
p_coord2
p_coord3
p_coord_bv_1
p_coord_bv_2
p_coord_bv_3
p_nvox_thr