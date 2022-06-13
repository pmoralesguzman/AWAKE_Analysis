%________________________________________________________________________
% Script to manually (not using Plotty) plot the density profiles of the
% experiment together with the ones from simulations
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 02/06/2021
%________________________________________________________________________

% close all;
clear;

% file location variables
datadirs = {'gm10','g0','gp10'};
datadirs_sim = {'gm10','g0','gp10'};

% datadirs = {'gm20'};
% datadirs_sim = {'gm20'};

dataformat = 'mat';
useAvg = false;
dump = 100;

% saving data
save_flag = true;
save_format = {'pdf','png'};

% plasma properties
plasmaden = 1.81e14;

% choose fields to plot
wakefields_direction = ''; % trans, long

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm for transverse, in ps for longitudinal)
trans_range = [0 0.16];
xi_range = [-20,618]; %ps

% choose property to plot
property_plot = 'density'; % density, wakefields, both

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose scaling
plot_scale = 'linear';

% choose if make pause or not
make_pause = false;

% figure number
fig_number = 4;

letters = {'a) $g = -0.93$\,\%/m','b) $g = -1$\,\%/m','c) $g = +0.03$\,\%/m',...
    'd) $g = 0$\,\%/m','e) $g = +0.87$\,\%/m','f) $g = +1$\,\%/m','g)','h)','i)','j)','k)','l)'};

% letters = {'a) g = 0.03 %/m',...
%     'b) g = 0 %/m'};

% directory to save the plots
plots_dir = ['gradsim_convergence/','field_density_longprofile_all/'];
% plots_dir = ['gradsim/','field_density_longprofile/'];

P = Plotty('datadir',datadirs{1},'dataformat',dataformat,...
    'useAvg',useAvg,...
    'plot_scale',plot_scale,...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'fig_number',fig_number,'title_flag',false);

fig_handle = figure(fig_number);

tile_handle = tiledlayout(fig_handle,length(datadirs),1,...
    'TileSpacing','compact','Padding','compact');

t1 = tiledlayout(tile_handle,4,1,'TileSpacing','none');
t1.Layout.Tile = 1;
t2 = tiledlayout(tile_handle,4,1,'TileSpacing','none');
t2.Layout.Tile = 2;
t3 = tiledlayout(tile_handle,4,1,'TileSpacing','none');
t3.Layout.Tile = 3;

% plot_position = [0.2857    0.0444    0.3979    0.9514];
% fig_handle.Units = 'normalized';
plot_position = [1,1,8.6*1.5,18];
fig_handle.Units = 'centimeters';
fig_handle.OuterPosition = plot_position; %[0.1 0.3 0.8 0.5]
sigma_rexp = 0.536; %mm

hletter = 0.84; % normalized
fontsize = 0.3; % cm
fontsize_den = 0.25; %cm

EDA = ExperimentalDataAnalyser('plasmaden',1.81e14);

%% Plot part

for d = 1:length(datadirs)
    
    EDA.datadir = datadirs{d};
    EDA.loadSCdata(); % P.SC
    z_plot = EDA.SCI_z; %ps
    r_plot = EDA.SCI_r*10; %mm
    
    %first manual trimming
    ind_r = (r_plot > -trans_range(2)*10) & (r_plot < trans_range(2)*10);
    density_plot = EDA.SCI(ind_r,:);
    r_plot = r_plot(ind_r);
    
    % stablish the opaque index
    % the standard deviation is used as a measure the avoid
    % noisy peaks that sets a wrong scale for the opaqueness
    % get 3 times the std deviation with no weights
    meanstd_density = 4*std(density_plot,[],'all');
    max_opaqueness = 1;
    ind_opaque = max_opaqueness*density_plot;
    ind_opaque(density_plot > meanstd_density) = max_opaqueness*meanstd_density;
    if d == 1
        max_indopaque = max(ind_opaque,[],'all');
    end
    ind_opaque = ind_opaque/max_indopaque;
    switch d
        case 1
            ax_density_exp(d) = nexttile(t1);
        case 2
            ax_density_exp(d) = nexttile(t2);
        case 3
            ax_density_exp(d) = nexttile(t3);
    end
    
    
    % set limits to the axis
    ax_density_exp(d).XLim = xi_range; %[min(P.SCxaxis),max(P.SCxaxis)];
    ax_density_exp(d).YLim = [min(r_plot),max(r_plot)];
    ax_density_exp(d).XDir = 'reverse';
    ax_density_exp(d).Box = 'on';
    %     ax_density_exp(d).YTickMode = 'manual';
    
    
    %----------------- PLOT EXP 2D ----------------------------------------
    imagesc(ax_density_exp(d),'XData',z_plot,'YData',r_plot,'CData',double(density_plot>0),'alphadata',ind_opaque);
    grad = [1 1 1; 0 0 0];
    colormap(ax_density_exp(d),grad);
    text(0.025,hletter,letters{2*(d-1)+1},'Units','normalized',...
        'FontUnits','centimeters','FontSize',fontsize,'interpreter','latex')
    
    if P.title_flag
        title(['propagation dist. = ',num2str(P.propagation_distance/100,3),' m',''])
    end
    
    yline(sigma_rexp,'r','LineWidth',1)
    yline(-sigma_rexp,'r','LineWidth',1)
    
    if d == 1
        ylabel('$x$ (mm)','FontUnits','centimeters','FontSize',fontsize,'Interpreter','latex')
    end
    
    ax_density_exp(d).XTickLabel = [];
    ax_density_exp(d).YTick = [-1,0,1];
    ax_density_exp(d).FontUnits = 'centimeters';
    ax_density_exp(d).FontSize = fontsize_den;
    
    
    %----------------- PLOT EXP LINEOUT -----------------------------------
        
    switch d
        case 1
            axlongprofile_exp(d) = nexttile(t1);
        case 2
            axlongprofile_exp(d) = nexttile(t2);
        case 3
            axlongprofile_exp(d) = nexttile(t3);
    end
%     axlongprofile_exp(d) = nexttile;
    % manual triming of experimental data
    
    ind_r = (r_plot >= -sigma_rexp) & (r_plot < sigma_rexp);
    density_trim = density_plot(ind_r,:);
    %end of manual triming
    
    long_profile = sum(density_trim);
    plongprofile_exp = plot(z_plot,long_profile,'k');
    ylim([0 5e4]);
    xlim(xi_range);
    
    axlongprofile_exp(d).XDir = 'reverse';    
    axlongprofile_exp(d).XTickLabel = [];
    axlongprofile_exp(d).YTickLabel = [];
    axlongprofile_exp(d).YTick = [];
    axlongprofile_exp(d).FontSize = fontsize;
    
    if d == 1
        ylabel({'density','(arb. units)'},'FontUnits','centimeters','FontSize',fontsize_den,'Interpreter','latex')
    end
    
    %---------------- PLOT SIM 2D -----------------------------------------
    
    timetemp = load([datadirs_sim{d},'/n130_m1350_slit.mat']);
    simtimeprofile = timetemp.densitymatrix;
    z_plot = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2))-2.6;
    r_plot = linspace(0,timetemp.ylims(2),size(simtimeprofile,1))*10;
    
    %first manual trimming
    ind_r = (r_plot > trans_range(1)*10) & (r_plot < trans_range(2)*10);
    density_plot_half = simtimeprofile(ind_r,:);
    r_plot = r_plot(ind_r);
    
    density_plot = [flip(density_plot_half);density_plot_half];
    meanstd_density = 3*std(density_plot,[],'all');
    max_opaqueness = 1;
    ind_opaque = max_opaqueness*density_plot;
    ind_opaque(density_plot > meanstd_density) = max_opaqueness*meanstd_density;
%     if d == 1
        max_indopaque_exp = max(ind_opaque,[],'all');
%     end

    ind_opaque = ind_opaque/max_indopaque_exp;
    
    switch d
        case 1
            ax_density_sim(d) = nexttile(t1);
        case 2
            ax_density_sim(d) = nexttile(t2);
        case 3
            ax_density_sim(d) = nexttile(t3);
    end
%     ax_density_sim(d) = nexttile;

    % set limits to the axis
    ax_density_sim(d).XLim = xi_range; %[min(P.SCxaxis),max(P.SCxaxis)];
    ax_density_sim(d).YLim = [-max(r_plot),max(r_plot)];
    ax_density_sim(d).XDir = 'reverse';
    ax_density_sim(d).Box = 'on';
    
    imagesc(ax_density_sim(d),'XData',z_plot,'YData',[-max(r_plot),max(r_plot)],'CData',double(density_plot>0),'alphadata',ind_opaque);
    grad = [1 1 1; 0 0 0];
    colormap(ax_density_sim(d),grad);
    text(0.025,hletter,letters{2*(d-1)+2},'Units','normalized',...
        'FontUnits','centimeters','FontSize',fontsize,'interpreter','latex')
    
    yline(sigma_rexp,'r','LineWidth',1)
    yline(-sigma_rexp,'r','LineWidth',1)
    
    ax_density_sim(d).XTickLabel = [];
%     ax_density_sim(d).XTick = [];
    ax_density_sim(d).YTick = [-1,0,1];
    ax_density_sim(d).FontUnits = 'centimeters';
    ax_density_sim(d).FontSize = fontsize_den;
    
    
    %---------------- PLOT SIM LINEOUT ------------------------------------
    switch d
        case 1
            axlongprofile_sim(d) = nexttile(t1);
        case 2
            axlongprofile_sim(d) = nexttile(t2);
        case 3
            axlongprofile_sim(d) = nexttile(t3);
    end
%     axlongprofile_sim(d) = nexttile;
    
    % manual triming of experimental data
    ind_r = (r_plot < sigma_rexp*10);
    density_trim = density_plot(ind_r,:);
    r_plot = r_plot(ind_r);
    
    %end of manual triming
    
    long_profile = sum(density_trim);
    plongprofile_sim = plot(z_plot,long_profile,'k');
    axlongprofile_sim(d).XDir = 'reverse';
    xlim(xi_range); ylim([0 2.9e7])
    
    if d ~= 3
%         axlongprofile_sim(d).XTick = [];
        axlongprofile_sim(d).XTickLabel = [];
        axlongprofile_sim(d).XTickLabel = {''};
    end
    axlongprofile_sim(d).YTickLabel = [];
    axlongprofile_sim(d).YTick = [];
end

xlabel(tile_handle,'$t$ (ps)','FontSize',9,'Interpreter','latex'); % 1 point = 1/72 inches = 0.0353 cm; 9 point = 0.388
drawnow;
pause(1);

P.fig_handle = fig_handle;
P.plot_name = 'allprofiles';
P.save_plot();
