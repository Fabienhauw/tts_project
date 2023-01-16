%to plot specified contrasts activation for one subj in one voxel. Also
%useful to get t and p values in the voxels and contrasts you want.
clc;clear;
xyzmm = [-45 -51 -10; -48 -44 23; -40 -41 46; -50 -16 50; -50 6 53; -70 -28 3; -60 12 -7; 50 -24 16; 65 4 0];

roi_names = {'VWFA','lSMG','lIPS', 'lprecent', 'lMFG', 'lSTG', 'lSTS', 'rSTG', 'rSTS'};

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17 

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
% gaucher_appar = {'Control02|Control04|Control05|Control07|Control17|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con_app = S;
S_con_app(mask_gauch_con) = [];
con_group = repmat({'control'}, 1, length(S_con_app));

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_droit = S;
S_droit(mask_gauch) = [];
syn_group = repmat({'syn'}, 1, length(S_droit));

categ = [syn_group, con_group];
nsub(1) = length(S_droit); % synesthetes
nsub(2) = length(S_con_app); % controls
ngroups = 2;

S_effect = [S_droit ; S_con_app];
vis_spec_con = [1:5];
% spec_con = [19;20;5]; %this is the t-con of interests;
% Ic=22; % this is the 'f contrast' with grouped conditions: this is useful
% if you have multiple beta to merge (for example high freq words). Pay
% attention to the order of "spec_con" and contrasts in "Ic".

for i_subj = 1 : length(S_effect)
    path_to_vis_stats{i_subj} = fullfile(D, S_effect(i_subj).name, 'Vis/loc/stats_s5_without_resting');
end

for tmp_con = 1 : length(vis_spec_con)
    vis_contrasts{tmp_con} = sprintf('con_000%d.nii', vis_spec_con(tmp_con));
end
vis_con_names = {'Faces', 'Houses', 'Tools', 'Numbers', 'Words'};

%% other parameters
% odd number of voxels defining the edge of the cube, centered on the chosen voxel, in which data will be averaged
cubesize=1;

%% display of exeriments as subplots
% grid of subplots
plot_cols=1;
%%% plot position and dimensions
x0=10;
y0=10;
width=1000;
height=600;

% save plots as image files
save_plots=0;
plot_filename='topdown_barchart'; % generic file name

% equate means across groups (within experiment)
equate_means=0; % leave at 0 or debug the plotting part below


number_coord = size(xyzmm,1);
for numcoord = 1 : number_coord %numcoord is for multiple coordinates...
    
    % groups dirs within rootdir:
    group_dir    = {'Synesthetes', 'Controls'};
    group_labels = group_dir; % labels for the plots (may be different from dir names)
    
    % visual run
    experiment_labels={'Visual'};
    
    % get coordinates
    load(fullfile(path_to_vis_stats{i_subj}, 'SPM.mat' )) ;
    iM = SPM.xVol.iM ;
    XYZ = iM( 1:3, : ) * [ xyzmm(numcoord,:)' ; 1 ] ;
    
    for i_subj = 1 : length(path_to_vis_stats)
        clear SPM
        cd (fullfile(path_to_vis_stats{i_subj}))
        load(fullfile(path_to_vis_stats{i_subj}, 'SPM.mat' )) ;

        if ~isfield(SPM,'VCbeta') % xSPM.STAT ~= 'P'
            %-Parameter estimates:   beta = xX.pKX*xX.K*y;
            %-Residual mean square: ResMS = sum(R.^2)/xX.trRV
            %----------------------------------------------------------------------
            beta_vis  = spm_get_data(SPM.Vbeta, XYZ);
            ResMS_vis = spm_get_data(SPM.VResMS,XYZ);
            Bcov_vis  = ResMS_vis*SPM.xX.Bcov;
        end
        
        %     X=SPM.xX.xKXs.X; % load the whitened/filtered design matrix
        %     c=SPM.xCon(Ic).c'; % load the weights of the contrast
        %     Res_header=spm_vol('ResMS.nii');
        %     vol_Res=spm_read_vols(Res_header);
        %     vol_SD = sqrt(vol_Res * (c * inv(X'*X) * c'));
        
        Ci    = 1.6449;  % = spm_invNcdf(1 - 0.05);
        for sc=1:length(vis_spec_con)
            cbeta_vis{i_subj}(sc, numcoord) = SPM.xCon(vis_spec_con(sc)).c'*beta_vis;
%             SE_vis{i_subj}(sc, numcoord)    = sqrt(diag(SPM.xCon(vis_spec_con(sc)).c'*Bcov_vis*SPM.xCon(vis_spec_con(sc)).c));
%             CI_vis{i_subj}(sc, numcoord)    = Ci*SE_vis{i_subj}(sc,numcoord);
        end
        cbeta_vis_mean{i_subj}(1, numcoord) = mean(cbeta_vis{i_subj}(:, numcoord));
        
        % normalization for each subject
        for sc=1:length(vis_spec_con)
            ncbeta_vis{i_subj}(sc, numcoord) = cbeta_vis{i_subj}(sc, numcoord) - cbeta_vis_mean{i_subj}(1, numcoord);
        end
        
        %     this part works if you give a f contrast for Ic (because only one
        %     contrast...)
        %     cbeta(:, numcoord) = SPM.xCon(Ic).c'*beta;
        %     SE(:, numcoord)    = sqrt(diag(SPM.xCon(Ic).c'*Bcov*SPM.xCon(Ic).c));
        %     CI(:, numcoord)    = Ci*SE(:,numcoord);
    end % end of subject loop
end

for g = 1 : ngroups
    if g>1
        startnumber = sum(nsub(1:g-1));
    else
        startnumber = 0;
    end
    for numcoord = 1 : number_coord
        for sc=1:length(vis_spec_con)
            cbeta_tot_vis{g}(sc, numcoord) = 0;
            num_moy_vis{g}(sc, numcoord) = 0;
            for i_subj = startnumber+(1:nsub(g))
                if ~isnan(cbeta_vis{i_subj}(sc, numcoord))
                    num_moy_vis{g}(sc, numcoord) = num_moy_vis{g}(sc, numcoord) + 1;
                    ncbeta_tot_vis{g}(sc, numcoord, i_subj) = ncbeta_vis{i_subj}(sc, numcoord);
                    cbeta_tot_vis{g}(sc, numcoord) = cbeta_tot_vis{g}(sc, numcoord) + cbeta_vis{i_subj}(sc, numcoord);
                end
            end
            
            %then you can calculate the mean for each group
            cbeta_tot_vis{g}(sc, numcoord) = cbeta_tot_vis{g}(sc, numcoord) / num_moy_vis{g}(sc, numcoord);
            
            %then you can calculate the std for each group, based on the
            %normalized cbeta
            std_vis{g}(sc, numcoord) = std(ncbeta_tot_vis{g}(sc, numcoord, :));
            CI_vis{g}(sc, numcoord) = std_vis{g}(sc, numcoord)/sqrt(num_moy_vis{g}(sc, numcoord));
        end
    end
end

grouporder = {};
for numcoord = 1 : number_coord
    grouporder = [grouporder sprintf('%s \n(%d %d %d)', roi_names{numcoord}, xyzmm(numcoord,1),xyzmm(numcoord,2),xyzmm(numcoord,3))];
end

% Col    = [ 1 0 0 ; 0 1 0 ; 0 0 1; 0 1 1; 1 0 1];
Col    = [ 1 1 0 ; 0.317 0.745 0.576 ; 0 0.133 0.576 ; 1 0.5 0];

final_fig = figure(31);
final_fig.WindowState = 'maximized';

% h = bar(cbeta');
all_y = [];
for numcoord=1:number_coord
    
    ax = subplot(1,number_coord,numcoord,'Parent',31);
    cla(ax)
    hold(ax,'on')
    % estimates
    %--------------------------------------------------------------
    clear y
    for g = 1 : ngroups
        for sc=1:length(vis_spec_con)
            y(sc,g) = cbeta_tot_vis{g}(sc, numcoord);
            ci(sc,g) = CI_vis{g}(sc, numcoord);
        end
    end
    all_y = [all_y y];
    h     = bar(y);
    [nbars, ngroup] = size(y);
    x = nan(nbars, ngroup);
    set(h(1), 'LineWidth',1, 'FaceColor', Col(4,:))
    set(h(2), 'LineWidth',1, 'FaceColor', Col(1,:))
    
    set(gca,'XTick',[1:length(xyzmm)],'XTickLabel',vis_con_names,...
        'FontSize', 6,'fontweight','n','LineWidth',1); % axes (fontweight can be b)
    xlabel(grouporder(numcoord), 'FontSize', 15)
    ylabel('Contrast value')
    legend(group_labels,'Location','northeast','Orientation','vertical')
    
    for i = 1:ngroup
        x(:,i) = h(i).XEndPoints;
    end
    
%     for g = 1 : ngroups
%         for sc=1:length(vis_spec_con)
%             line([x(sc,g) x(sc,g)],([ci(sc,g) -ci(sc,g)] + y(sc,g)),...
%                 'LineWidth',1, 'Color', [0.3 0.3 0.3])
%         end
%     end
    errorbar(x,y,ci,'k','linestyle','none');%Adding the errorbars
    
    hold on
end

for numcoord=1:number_coord
    figure(31);
    subplot(1,number_coord,numcoord);
    mi = min(all_y(:));
    ma = max(all_y(:));
    ylim([mi - (ma-mi)*0.05     ma + (ma-mi)*0.05 ]);
end

set (final_fig,'Color','w'); % background

        
hold off
if save_plots
    saveas(fn, [plot_filename '-' num2str(fn) '.png'],'png')
end



