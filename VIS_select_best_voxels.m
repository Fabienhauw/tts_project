clear;

%% parameters:
select_con = [9]; % 9 = words - (faces + houses) in the visual GLM; %select_con = number of the contrast used to select the voxels (in SPM.mat of the select_path)
test_con = [16]; % 12 = w-(f+h+t), 13 = f-(w+h+t), 14 = h-(f+w+t), 15 = t-(f+h+w), 16 = n-(f+h+t) %test_con = number of the contrast used to test the voxels (in SPM.mat of the select_path)
diff_ref = 0; % 1 = use a different contrast for selecting the best_voxels and for comparing signals;
only_right_hand = 1;
volume_option = 2; %1 = ROI, 2 = coordinates;
sphereradius = 10; %% in  mm a 10 mm sphere focuses tightly on the VWFA locationend, useful only for volume option = 2;
voxel_option = 1; % 1 = keep n best voxel; 2 = keep n voxel above a threshold;
per_vox = 10; % percentage of voxels you want to keep
pvalue = 0.01;
flipsign = 1; %% default is to keep the contrast "as is" with its sign

%%
tts_group.path_to_subject;

if only_right_hand
    gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
    mask_gauch =  ~cellfun(@isempty,(regexp({S.name},{'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control'})));
    res_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Vis/loc/RH_comparison_best_voxels';
else
    gaucher_appar = {'Control02|Control04|Control07|Control17|Sujet'};
    mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
    res_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Vis/loc/comparison_best_voxels';
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
    };
%
ROInames = {'SMA','PRECUN'};

if volume_option == 1
    all_xyzmm = ROInames;
elseif volume_option == 2
    %different coordinates, according to where they come from
    if test_con == 12
        peaks = 'vis_12';
        coord_names = {'lOcc', 'rOcc', 'lIPS', 'rIPS', 'SMA', 'mIFG', 'iIFG', 'VWFA', 'lSTS', 'rSTS'};
        all_xyzmm = [-20 -94 -4; 20 -88 -4; -30 -48 43; 38 -54 46; 0 12 53; -40 6 30; -50 36 13; -48 -54 -20; -68 -44 6; 58 -31 0]; %% VISUAL WORDS-(h+f+t), RH
    elseif test_con == 13
        peaks = 'vis_13';
        coord_names = {'FFA', 'OFA', 'lFFA'};
        all_xyzmm = [42 -56 -20; 52 -64 16; -38 -51 -20]; %% VISUAL FACES-(w+h+t), rh, no difference.
    elseif test_con == 14
        peaks = 'vis_14';
        coord_names = {'post rPPA', 'post lPPA', 'ant rPPA', 'ant lPPA'};
        all_xyzmm = [22 -76 -12; -18 -78 -12; 28 -54 -7; -28 -46 -7]; %% VISUAL HOUSES-(w+f+t), rh, no difference.
    elseif test_con == 15
        peaks = 'vis_15';
        coord_names = {'lPPA', 'rPPA', 'lOTA', 'rOTA'};
        all_xyzmm = [-25 -54 12; 30 -48 -12; -45 -68 -2; 48 -64 -7]; %% VISUAL TOOLS-(h+f+w), rh, no difference.
    elseif test_con == 16
        peaks = 'vis_16';
        coord_names = {'rIPS', 'lIPS'};
        all_xyzmm = [32 -61 43; -45 -41 43]; %% VISUAL NUMBERS-(h+f+t), rh, no difference.
    end
    
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
path_to_visual_stats = {};
path_to_scan = {};
path_to_stats_mask = {};
subjs={};
for k=1:numel(S)
    subjs                   = [subjs;S(k).name];
    path_to_subj            = fullfile(D, S(k).name);
    path_to_scan            = [path_to_scan; path_to_subj];
    path_to_visual_sub      = fullfile(D, S(k).name, 'Vis/loc/stats_s5');
    path_to_subj            = fullfile(D, S(k).name, 'Vis/loc/stats_s5');
    path_to_mask            = fullfile(D, S(k).name, 'Vis/loc/stats_s5/mask.nii');
    path_to_visual_stats    = [path_to_visual_stats;path_to_visual_sub];
    path_to_stats           = [path_to_stats;path_to_subj];
    path_to_stats_mask      = [path_to_stats_mask;path_to_mask];
end

experiments = struct(...
    'select_path',  path_to_visual_stats,...  % subpath to the localizer SPM.mat inside each subject
    'test_path',    path_to_stats,...  % subpath to the test data SPM.mat inside each subject
    'test_mask',    path_to_stats_mask,... % subpath to the test data SPM.mat mask for each subject
    'data',         subjs...
    );
if ~diff_ref
    for tmp_path = 1 : size(experiments,1)
        experiments(tmp_path).select_path = experiments(tmp_path).test_path;
    end
end

str_path = strsplit(experiments(1).select_path, '/');
if sum(~cellfun(@isempty,(regexp(str_path, 'Aud'))))
    ref_mod = 'auditive';
elseif sum(~cellfun(@isempty,(regexp(str_path, 'Vis'))))
    ref_mod = 'visual';
end

if ~isdir(res_dir)
    mkdir(res_dir)
end

for zz=1:length(test_con)
    current_con = test_con(zz);
    if ~diff_ref
        select_con = current_con;
    end
    for qq=1:size(all_xyzmm,1)
        
        if volume_option == 1
            xyzmm = sprintf('%s',ROInames{qq});
        elseif volume_option == 2
            xyzmm = all_xyzmm(qq,:);
        end
        
        %% start of computation of the results
        totsub = size(experiments,1);
        spmfiles={};
        spmfiles_select={};
        spmfiles_test={};
        
        for nsub=1:totsub
            spmfiles_select{nsub}   = fullfile(experiments(nsub).select_path,'SPM.mat');
            spmfiles_test{nsub}     = fullfile(experiments(nsub).test_path,'SPM.mat');
        end
        
        %% define the search volume
        if (volume_option == 1)
            load(spmfiles_select{1}); %%% load one SPM.mat just to get the transformation matrix iM
            ROIheader   = spm_vol(ROIfiles{qq});
            ROIvol      = spm_read_vols(ROIheader);
            ROInumbers  = unique(round(ROIvol(ROIvol>0)));
        elseif (volume_option == 2)
            if size(xyzmm,1)==1
                xyzmm = xyzmm';
            end
            load(spmfiles_select{1}); %%% load one SPM.mat just to get the transformation matrix iM
            iM = SPM.xVol.iM ;
            xyzvox = iM( 1:3, : ) * [ xyzmm ; 1 ] ;
            spherecoords = xyzmm';
            
            %% the following code is just to get the voxel sizes
            select_tfile = fullfile(experiments(nsub).select_path,sprintf('spmT_%04d.nii',select_con));
            %             select_tfile = fullfile(experiments(nsub).select_path,sprintf('con_%04d.nii',select_con));
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
                            spherevol_gen(x,y,z)=1;
                        else
                            spherevol_gen(x,y,z)=0;
                        end
                    end
                end
            end
        end
        
        fprintf( 'FIND BEST VOXELS %d/%d\n\n', qq + (size(all_xyzmm,1) * (zz - 1)), size(all_xyzmm,1) * length(test_con)) ;
        h = waitbar( 0, 'Finding best voxels...' ) ;
        
        for i_subj = 1:totsub  %%%% loop across subjects
            clear tvol
            clear searchvol
            
            %%% extract localizer t-test to optimize for this subject:
            clear select_tvol
            select_tfile    = fullfile(experiments(i_subj).select_path,sprintf('spmT_%04d.nii',select_con));
            %             select_tfile    = fullfile(experiments(i_subj).select_path,sprintf('con_%04d.nii',select_con));
            select_theader  = spm_vol(select_tfile);
            select_tvol     = spm_read_vols(select_theader);
            
            %%% read all the test images
            clear test_convol;
            
            t                   = test_con;
            test_confile        = fullfile(experiments(i_subj).test_path,sprintf('spmT_%04d.nii',current_con));
            %             test_confile        = fullfile(experiments(i_subj).test_path,sprintf('con_%04d.nii',current_con));
            test_conheader      = spm_vol(test_confile);
            test_convol         = spm_read_vols(test_conheader);
            
            %%% read the mask.nii image
            mask_confile        = experiments(i_subj).test_mask;
            mask_conheader      = spm_vol(mask_confile);
            mask_convol         = spm_read_vols(mask_conheader);
            
            for x=1:size(spherevol_gen,1)
                for y=1:size(spherevol_gen,2)
                    for z=1:size(spherevol_gen,3)
                        if mask_convol(x,y,z)==1 & spherevol_gen(x,y,z)==1 % to exclude voxels which are outside the statistical map mask.nii
                            spherevol(x,y,z)=1;
                        else
                            spherevol(x,y,z)=0;
                        end
                    end
                end
            end
            
            waitbar( (qq+(i_subj-1)*length(xyzmm))/(totsub*length(xyzmm)), h ) ;
            if (volume_option == 1)
                anal(qq).volume_str = sprintf('Region %d (%s)', qq, ROInames{qq});
                %                     searchvol = double(round(ROIvol)==ROInumbers(qq));
                searchvol = double(round(ROIvol)==1);
            elseif (volume_option == 2)
                anal(qq).volume_str = sprintf('%d mm sphere centered at [%d %d %d]',sphereradius,xyzmm);
                searchvol = spherevol;
            end
            anal(qq).roi_name = coord_names{qq};
            anal(qq).peaks_origin = peaks;
            
            select_tvol(searchvol(:)==0) = Inf; % every voxel outside of the mask is set to 'Inf'
            tvol = flipsign * select_tvol(:) .* (searchvol(:)>0);
            
            %write an image of the sphere where you selected your voxels;
            test_convol_copy = test_convol;
            test_convol_copy(searchvol(:)==0) = 0;
            test_convol_copy(searchvol(:)~=0) = 1;
            test_conheader_copy = test_conheader;
            test_conheader_copy.fname = sprintf('%s_%dmm_sphere_in_%d_%d_%d_%s.nii',test_conheader.fname(1:end-4), sphereradius, xyzmm, ref_mod);
            spm_write_vol(test_conheader_copy,test_convol_copy);
            
            
            if (voxel_option == 1)
                [tvalues,xyz] = sort(tvol,'descend');
                xyz(isnan(tvalues)) = [];
                tvalues(isnan(tvalues)) = [];
                nvox = round(length(find(searchvol>0))*(per_vox/100)); % determinate how many vox is X% of the ROI (searchvol >0)
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
            
            anal(qq).nbvoxels(i_subj) = nvox;
            anal(qq).nbvoxels_thr(i_subj) = nvox_thr;
            if nvox>0
                [x,y,z]=ind2sub(size(select_tvol),xyz);
                anal(qq).voxels{i_subj}             = [ x y z ]';  %%% store the voxels for this subject
                anal(qq).xyz{i_subj}                = xyz;  %%% this is the easiest storage: direct indexes to the data matrices, can be used as a data selector
                anal(qq).coordsmm(1, i_subj)        = mean(x * select_theader.mat(1,1) + select_theader.mat(1,4)); %%% mean coordinates of the identified voxels
                anal(qq).coordsmm(2, i_subj)        = mean(y * select_theader.mat(2,2) + select_theader.mat(2,4));
                anal(qq).coordsmm(3, i_subj)        = mean(z * select_theader.mat(3,3) + select_theader.mat(3,4));
                anal(qq).bvcoordsmm(1, i_subj)      = (x(1) * select_theader.mat(1,1) + select_theader.mat(1,4)); %%% mean coordinates of the identified voxels
                anal(qq).bvcoordsmm(2, i_subj)      = (y(1) * select_theader.mat(2,2) + select_theader.mat(2,4));
                anal(qq).bvcoordsmm(3, i_subj)      = (z(1) * select_theader.mat(3,3) + select_theader.mat(3,4));
                anal(qq).coords(1, i_subj)          = mean(x);
                anal(qq).coords(2, i_subj)          = mean(y);
                anal(qq).coords(3, i_subj)          = mean(z);
                indiv_coords(i_subj,:)              = anal(qq).coordsmm(:, i_subj);
                anal(qq).all_activation{i_subj}     = test_convol(xyz);
                anal(qq).activation(i_subj)         = mean(test_convol(xyz));
                anal(qq).group(i_subj)              = group_S(i_subj);
            end
            anal(qq).ref_con_path{i_subj,1} = fullfile(experiments(i_subj).select_path, sprintf('spmT_%04d.nii',select_con));
            anal(qq).test_con_path{i_subj,1} = fullfile(experiments(i_subj).test_path, sprintf('spmT_%04d.nii',current_con));
            
            if voxel_option == 1
                anal(qq).voxel_str = sprintf('The %d best voxels for localizer contrast %d based on %s contrast %d',nvox,current_con, ref_mod, select_con);
            elseif voxel_option == 2
                anal(qq).voxel_str = sprintf('All voxels for localizer contrast %d with p-value < %7.5f based on %s contrast %d', current_con, pvalue, ref_mod, select_con);
            end
            
            % to write a .nii map with only voxels kepts:
            test_convol_copy = test_convol;
            test_convol_copy(:) = 0;
            test_convol_copy(xyz) = 1;
            if volume_option==1
                test_conheader.fname = sprintf('%s_best_vox_in %s_based_on_%s_con_%d.nii',test_conheader.fname(1:end-4), xyzmm, ref_mod, select_con);
                test_conheader.descrip = sprintf('%s thresholded to 10 prc best voxels in region %s based on visual con %d',test_conheader.descrip, xyzmm, select_con);
            elseif volume_option ==2
                test_conheader.fname = sprintf('%s_best_vox_sph_%d_%d_%d_based_on_%s_con_%d.nii', test_conheader.fname(1:end-4),xyzmm, ref_mod, select_con);
                test_conheader.descrip = sprintf('%s thresholded to 10 prc best voxels in 10 mm rad sphere centered on %d %d %d',test_conheader.descrip, xyzmm);
            end
            
            spm_write_vol(test_conheader,test_convol_copy);
            
            
        end %%% subject loop
        close(h)
        
    end % coord loop
    
    if diff_ref
        save_name = sprintf('best_voxels_in_vis_spmT%04d_based_on_%s_spmT%04d_%s_peaks', current_con, ref_mod, select_con, peaks)
    else
        save_name = sprintf('best_voxels_in_vis_spmT%04d_%s_peaks', current_con, peaks)
    end
    
    if volume_option == 1
        save_name = sprintf('best_voxels_in_vis_spmT%04d_ROI.mat', current_con);
    elseif volume_option == 2
        save_name = sprintf('%s_coord.mat',save_name);
    end
    save(fullfile(res_dir, save_name), 'anal')
end % contrast loop

% return

%% now the part to plot
cd(res_dir)

do_plot = 0;

results = dir('best*vis*');
count_fig = 0;
for tmp_anal = 1 : size(results, 1)
    clear anal
    load(results(tmp_anal).name)
    clear matcon Group_Id Subj_Id vox_activ
    for j = 1 : size(anal,2) % j = ROI;
        count = 0;
        for k = 1 : size(anal(j).all_activation,2) % k = subjs;
            for l = 1 : size(anal(j).all_activation{k},1)
                count = count + 1;
                matcon{j}{count, 1}     = anal(j).group(k);
                Group_Id{j}{count, 1}   = anal(j).group(k);
                matcon{j}{count, 2}     = k;
                Subj_Id{j}{count, 1}    = k;
                matcon{j}{count, 3}     = anal(j).all_activation{k}(l);
                vox_activ{j}{count, 1}  = anal(j).all_activation{k}(l);
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
        tmp_mask = cellfun(@isempty,(anal(s).all_activation(:)));
        anal(s).activation(tmp_mask) = '';
        anal(s).coordsmm(:,tmp_mask) = '';
        anal(s).bvcoordsmm(:,tmp_mask) = '';
        anal(s).group_mask = anal(s).group(~tmp_mask);
        
        [h{tmp_anal,s},p{tmp_anal,s},ci{tmp_anal,s},stats{tmp_anal,s}] = ...
            ttest2(...
            anal(s).activation(anal(s).group_mask == 1), ...
            anal(s).activation(anal(s).group_mask == 2));
        tval{tmp_anal,s} = stats{tmp_anal,s}.tstat;
        anal(s).p_val = p{tmp_anal,s};
        anal(s).t_val = tval{tmp_anal,s};
        
        [h_coord1{tmp_anal,s},p_coord1{tmp_anal,s},ci_coord1{tmp_anal,s},stats_coord1{tmp_anal,s}] = ...
            ttest2(...
            anal(s).coordsmm(1,anal(s).group_mask == 1), ...
            anal(s).coordsmm(1,anal(s).group_mask == 2));
        tval_coord1{tmp_anal,s} = stats_coord1{tmp_anal,s}.tstat;
        anal(s).p_coord1 = p_coord1{tmp_anal,s};
        anal(s).tval_coord1 = tval_coord1{tmp_anal,s};
        
        [h_coord2{tmp_anal,s},p_coord2{tmp_anal,s},ci_coord2{tmp_anal,s},stats_coord2{tmp_anal,s}] = ...
            ttest2(...
            anal(s).coordsmm(2,anal(s).group_mask == 1), ...
            anal(s).coordsmm(2,anal(s).group_mask == 2));
        tval_coord2{tmp_anal,s} = stats_coord2{tmp_anal,s}.tstat;
        anal(s).p_coord2 = p_coord2{tmp_anal,s};
        anal(s).tval_coord2 = tval_coord2{tmp_anal,s};
        
        [h_coord3{tmp_anal,s},p_coord3{tmp_anal,s},ci_coord3{tmp_anal,s},stats_coord3{tmp_anal,s}] = ...
            ttest2(...
            anal(s).coordsmm(3,anal(s).group_mask == 1), ...
            anal(s).coordsmm(3,anal(s).group_mask == 2));
        tval_coord3{tmp_anal,s} = stats_coord3{tmp_anal,s}.tstat;
        anal(s).p_coord3 = p_coord3{tmp_anal,s};
        anal(s).tval_coord3 = tval_coord3{tmp_anal,s};
        
        [h_coord_bv_1{tmp_anal,s},p_coord_bv_1{tmp_anal,s},ci_coord_bv_1{tmp_anal,s},stats_coord_bv_1{tmp_anal,s}] = ...
            ttest2(...
            anal(s).bvcoordsmm(1,anal(s).group_mask == 1), ...
            anal(s).bvcoordsmm(1,anal(s).group_mask == 2));
        tval_coord_bv_1{tmp_anal,s} = stats_coord_bv_1{tmp_anal,s}.tstat;
        anal(s).p_coord_bv_1 = p_coord_bv_1{tmp_anal,s};
        anal(s).tval_coord_bv_1 = tval_coord_bv_1{tmp_anal,s};
        
        [h_coord_bv_2{tmp_anal,s},p_coord_bv_2{tmp_anal,s},ci_coord_bv_2{tmp_anal,s},stats_coord_bv_2{tmp_anal,s}] = ...
            ttest2(...
            anal(s).bvcoordsmm(2,anal(s).group_mask == 1), ...
            anal(s).bvcoordsmm(2,anal(s).group_mask == 2));
        tval_coord_bv_2{tmp_anal,s} = stats_coord_bv_2{tmp_anal,s}.tstat;
        anal(s).p_coord_bv_2 = p_coord_bv_2{tmp_anal,s};
        anal(s).tval_coord_bv_2 = tval_coord_bv_2{tmp_anal,s};
        
        [h_coord_bv_3{tmp_anal,s},p_coord_bv_3{tmp_anal,s},ci_coord_bv_3{tmp_anal,s},stats_coord_bv_3{tmp_anal,s}] = ...
            ttest2(...
            anal(s).bvcoordsmm(3,anal(s).group_mask == 1), ...
            anal(s).bvcoordsmm(3,anal(s).group_mask == 2));
        tval_coord_bv_3{tmp_anal,s} = stats_coord_bv_3{tmp_anal,s}.tstat;
        anal(s).p_coord_bv_3 = p_coord_bv_3{tmp_anal,s};
        anal(s).tval_coord_bv_3 = tval_coord_bv_3{tmp_anal,s};
        
        [h_nvox{tmp_anal,s},p_nvox{tmp_anal,s},ci_nvox{tmp_anal,s},stats_nvox{tmp_anal,s}] = ...
            ttest2(...
            anal(s).nbvoxels(anal(s).group == 1), ...
            anal(s).nbvoxels(anal(s).group == 2));
        tval_nvox{tmp_anal,s} = stats_nvox{tmp_anal,s}.tstat;
        anal(s).p_nvox = p_nvox{tmp_anal,s};
        anal(s).tval_nvox = tval_nvox{tmp_anal,s};
        
        [h_nvox_thr{tmp_anal,s},p_nvox_thr{tmp_anal,s},ci_nvox_thr{tmp_anal,s},stats_nvox_thr{tmp_anal,s}] = ...
            ttest2(...
            anal(s).nbvoxels_thr(anal(s).group == 1), ...
            anal(s).nbvoxels_thr(anal(s).group == 2));
        tval_nvox_thr{tmp_anal,s} = stats_nvox_thr{tmp_anal,s}.tstat;
        anal(s).p_nvox_thr = p_nvox_thr{tmp_anal,s};
        anal(s).tval_nvox_thr = tval_nvox_thr{tmp_anal,s};
    end
    
    save(fullfile(res_dir, results(tmp_anal).name), 'anal')
    
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
            
            title(sprintf('plot of the %d prc best voxels in %s', per_vox, anal(tmp_grap).volume_str));
        end
    end
end

for tmp_stat = 1 : size(results, 1)
    tmp_con = regexp(results(tmp_stat).name, 'spmT\d\d\d\d', 'match');
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