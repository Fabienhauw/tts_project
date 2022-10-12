clear;
clc;

tts_group.path_to_subject;
only_right_hand = 1;

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

%% parameters:
select_con = [9]; % 9 = words - (faces + houses) in the visual GLM; %select_con = number of the contrast used to select the voxels (in SPM.mat of the select_path)
test_con = [11;22]; % 11 = normal - scrambled in the auditive script; 22 = (normal + numbers + words + PW) - scr speech; %test_con = number of the contrast used to select the voxels (in SPM.mat of the select_path)
diff_ref = 0; % 1 = use a different contrast for selecting the best_voxels and for comparing signals;
volume_option = 2; %1 = ROI, 2 = coordinates;
sphereradius = 10; %% in  mm a 10 mm sphere focuses tightly on the VWFA locationend, useful only for volume option = 2;
voxel_option = 1; % 1 = keep n best voxel; 2 = keep n voxel above a threshold;
per_vox = 10; % percentage of voxels you want to keep
pvalue = 0.001;
flipsign = 1; %% default is to keep the contrast "as is" with its sign

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
% ROInames = {'SMA','PRECUN','IFG', 'SMG', 'VWFA', 'VWFA2', 'VWFA3'};
if volume_option == 1
    all_xyzmm = ROInames;
elseif volume_option == 2
    %different coordinates, according to where they come
%     all_xyzmm = [-58 -44 23; -52 -51 -20; -45 22 23]; %% AUD normal > scr, RH synesthetes,  SMG, VWFA, IFG.
%     all_xyzmm = [-58 -44 23; -40 -24 -22; -50 -8 46]; %% AUD normal > scr, all synesthetes,  SMG, VWFA, IFG, no difference.
%     all_xyzmm = [-50 -46 26; -45 -51 -12; -40 4 33]; %% VISUAL w-(h+f), all synesthetes, SMG, VWFA, IFG, no difference.
%     all_xyzmm = [-50 -46 26; -48 -54 -20; -40 6 30]; %% VISUAL w-(h+f), RH synesthetes, SMG, VWFA, IFG, no difference.
    all_xyzmm = [-48 -41 18; -42 -51 -14; -50 -8 43]; %% AUD "big" phonology, RH synesthetes, SMG, VWFA, IFG, no difference.
%     all_xyzmm = [-48 -41 18; -42 -51 -14; -50 -8 46]; %% AUD "big" phonology, all synesthetes, SMG, VWFA, IFG, no difference.
end

%% paths and scans
scans ={};
for k = 1:length(S)
    vol_name = fullfile(D, S(k).name,'Aud/loc/stats_s4/spmT_0011.nii');
    scans = [scans;vol_name];
    if ~isempty(regexp(S(k).name, 'Sujet'))
        group_S(k)    = 1; % 1 is for synesthetes, 2 for controls
    elseif ~isempty(regexp(S(k).name, 'Control'))
        group_S(k)    = 2; % 1 is for synesthetes, 2 for controls
    end
end

fff = scans;

% check file number of fff
loopn_fff=size(fff,1);

%
pouf=1;

path_to_stats = {};
path_to_visual_stats = {};
path_to_scan = {};
subjs={};
for k=1:numel(S)
    subjs                   = [subjs;S(k).name];
    path_to_subj            = fullfile(D, S(k).name);
    path_to_scan            = [path_to_scan; path_to_subj];
    path_to_visual_sub      = fullfile(D, S(k).name, 'Vis/loc/stats_s5');
    path_to_subj            = fullfile(D, S(k).name, 'Aud/loc/stats_s5');
    path_to_visual_stats    = [path_to_visual_stats;path_to_visual_sub];
    path_to_stats           = [path_to_stats;path_to_subj];
end

experiments = struct(...
    'select_path',  path_to_visual_stats,...  % subpath to the localizer SPM.mat inside each subject
    'test_path',    path_to_stats,...  % subpath to the test data SPM.mat inside each subject
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

if diff_ref
    res_dir = sprintf('%scomparison_best_voxels_diff_reference/', res_dir);
else
    res_dir = sprintf('%scomparison_best_voxels/', res_dir);
end

if volume_option == 1
    res_dir = fullfile(res_dir, 'ROI');
elseif volume_option == 2
    res_dir = fullfile(res_dir, 'coord');
end

if ~isdir(res_dir)
    mkdir(res_dir)
end

for zz=1:length(test_con)
    current_con = test_con(zz);
    if ~diff_ref
        select_con = current_con;
    end
    for qq=1:length(all_xyzmm)
        
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
            select_theader = spm_vol(select_tfile);
            select_tvol=spm_read_vols(select_theader);
            
            ROInames{1} =  'sphere';
            %%%% define sphere
            disp('computing sphere voxels.... ');
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
        
        fprintf( '\n FIND BEST VOXELS\n' ) ;
        h = waitbar( 0, 'Finding best voxels...' ) ;
        
        for i_subj = 1:totsub  %%%% loop across subjects
            
            %%% extract localizer t-test to optimize for this subject:
            clear select_tvol
            select_tfile    = fullfile(experiments(i_subj).select_path,sprintf('spmT_%04d.nii',select_con));
            select_theader  = spm_vol(select_tfile);
            select_tvol     = spm_read_vols(select_theader);
            
            %%% read all the test images
            clear test_convol;
            
            t                   = test_con;
            test_confile        = fullfile(experiments(i_subj).test_path,sprintf('spmT_%04d.nii',current_con));
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
            
            
            
            select_tvol(searchvol(:)==0) = Inf; % every voxel outside of the mask is set to 'Inf'
            tvol = flipsign * select_tvol(:) .* (searchvol(:)>0);
            
            if (voxel_option == 1)
                [tvalues,xyz] = sort(tvol,'descend');
                xyz(isnan(tvalues)) = [];
                tvalues(isnan(tvalues)) = [];
                nvox = round(length(find(searchvol>0))*(per_vox/100));
                xyz = xyz(1:nvox);
            elseif (voxel_option == 2)
                load(spmfiles_select{i_subj});
                criticalt = tinv(1-pvalue,SPM.xX.erdf);
                xyz = find(tvol>=criticalt);
                nvox = length(find(xyz>0));
            end
            
            anal{qq}.nbvoxels(i_subj) = nvox;
            if nvox>0
                [x,y,z]=ind2sub(size(select_tvol),xyz);
                anal{qq}.voxels{i_subj}         = [ x y z ]';  %%% store the voxels for this subject
                anal{qq}.xyz{i_subj}            = xyz;  %%% this is the easiest storage: direct indexes to the data matrices, can be used as a data selector
                anal{qq}.coordsmm(1, i_subj)    = mean(x * select_theader.mat(1,1) + select_theader.mat(1,4)); %%% mean coordinates of the identified voxels
                anal{qq}.coordsmm(2, i_subj)    = mean(y * select_theader.mat(2,2) + select_theader.mat(2,4));
                anal{qq}.coordsmm(3, i_subj)    = mean(z * select_theader.mat(3,3) + select_theader.mat(3,4));
                anal{qq}.coords(1, i_subj)    = mean(x);
                anal{qq}.coords(2, i_subj)    = mean(y);
                anal{qq}.coords(3, i_subj)    = mean(z);
                indiv_coords(i_subj,:)          = anal{qq}.coordsmm(:, i_subj);
                anal{qq}.all_activation{i_subj} = test_convol(xyz);
                anal{qq}.activation(i_subj)     = mean(test_convol(xyz));
                anal{qq}.group(i_subj)          = group_S(i_subj);
            end
            
            if voxel_option == 1
                anal{qq}.voxel_str = sprintf('The %d best voxels for localizer contrast %d based on %s contrast %d',nvox,current_con, ref_mod, select_con);
            elseif voxel_option == 2
                anal{qq}.voxel_str = sprintf('All voxels for localizer contrast %d with p-value < %7.5f based on %s contrast %d', current_con, pvalue, ref_mod, select_con);
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
        if volume_option==1
            aga{pouf,1}=xyzmm;
            aga{pouf,2}=current_con;
            aga{pouf,3}=nvox;
        elseif volume_option==2
            aga(pouf,[1:3])=xyzmm;
            aga(pouf,4)=current_con;
            aga(pouf,5)=nvox;
        end
        pouf=pouf+1;
    end % coord loop
    if volume_option == 1
        save(sprintf('%s/best_voxels_in_aud_spmT%04d.mat',...
            res_dir, current_con), 'anal');
    elseif volume_option == 2
        save(sprintf('%s/best_voxels_in_aud_spmT%04d.mat',...
            res_dir, current_con), 'anal');
    end
    
end % contrast loop

return

%% now the part to plot
cd(res_dir)

do_plot = 0;

results = dir('best*');
count_fig = 0;
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
        anal{1,s}.activation(tmp_mask) = '';
        anal{1,s}.group_mask = anal{1,s}.group(~tmp_mask);
        [h{tmp_anal,s},p{tmp_anal,s},ci{tmp_anal,s},stats{tmp_anal,s}] = ...
            ttest2(anal{1,s}.activation(anal{1,s}.group_mask == 1), anal{1,s}.activation(anal{1,s}.group_mask == 2));
        tval{tmp_anal,s} = stats{tmp_anal,s}.tstat;
        
        [h_coord1{tmp_anal,s},p_coord1{tmp_anal,s},ci_coord1{tmp_anal,s},stats_coord1{tmp_anal,s}] = ...
            ttest2(anal{1,s}.coordsmm(1,anal{1,s}.group_mask == 1), anal{1,s}.coordsmm(1,anal{1,s}.group_mask == 2));
        tval_coord1{tmp_anal,s} = stats_coord1{tmp_anal,s}.tstat;
        
        [h_coord2{tmp_anal,s},p_coord2{tmp_anal,s},ci_coord2{tmp_anal,s},stats_coord2{tmp_anal,s}] = ...
            ttest2(anal{1,s}.coordsmm(2,anal{1,s}.group_mask == 1), anal{1,s}.coordsmm(2,anal{1,s}.group_mask == 2));
        tval_coord2{tmp_anal,s} = stats_coord2{tmp_anal,s}.tstat;
        
        [h_coord3{tmp_anal,s},p_coord3{tmp_anal,s},ci_coord3{tmp_anal,s},stats_coord3{tmp_anal,s}] = ...
            ttest2(anal{1,s}.coordsmm(3,anal{1,s}.group_mask == 1), anal{1,s}.coordsmm(3,anal{1,s}.group_mask == 2));
        tval_coord3{tmp_anal,s} = stats_coord3{tmp_anal,s}.tstat;
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

fprintf('\np-values for con %d: ', test_con(1))
fprintf('%4.4f  ', p{1,:});
fprintf('\nt-values for con %d: ', test_con(1))
fprintf('%4.4f  ', tval{1,:});fprintf('\n')
fprintf('\np-values for con %d: ', test_con(2))
fprintf('%4.4f  ', p{2,:});
fprintf('\nt-values for con %d: ', test_con(2))
fprintf('%4.4f  ', tval{2,:}); fprintf('\n')
p_coord1
p_coord2
p_coord3
