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
fig_number = 3;

letters = {'a) $g = -0.93$\,\%/m','a) $g = -1$\,\%/m','c) $g = +0.03$\,\%/m',...
    'b) $g = 0$\,\%/m','e) $g = +0.87$\,\%/m','c) $g = +1$\,\%/m','g)','h)','i)','j)','k)','l)'};

% letters = {'a) g = 0.03 %/m',...
%     'b) g = 0 %/m'};

% directory to save the plots
plots_dir = ['gradsim_eps/','density_longprofile_all/'];
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

tile_handle = tiledlayout(fig_handle,2*length(datadirs),1);
tile_handle.TileSpacing = 'none';
tile_handle.Padding = 'compact';
% plot_position = [0.2857    0.0444    0.3979    0.9514];
% fig_handle.Units = 'normalized';
plot_position = [15.5340    8.2785   12.8999   11.9944];
fig_handle.Units = 'centimeters';
fig_handle.OuterPosition = plot_position; %[0.1 0.3 0.8 0.5]
sigma_rexp = 0.536; %mm

hletter = 0.84; % normalized
fontsize = 0.3; % cm
fontsize_den = 0.25; %cm


%% Plot part

for d = 1:length(datadirs)
    
    
    %---------------- PLOT SIM 2D -----------------------------------------
    
    timetemp = load([datadirs_sim{d},'/n100_m1350_slit.mat']);
    simtimeprofile = timetemp.densitymatrix;
    z_plot = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2));
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
    if d == 1
        max_indopaque_exp = max(ind_opaque,[],'all');
    end
    
    ind_opaque = ind_opaque/max_indopaque_exp;
    ax_density_sim(d) = nexttile;
    
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
    
    if d == 1
        ylabel('x (mm)','FontUnits','centimeters','FontSize',fontsize)
    end
    
    ax_density_sim(d).XTickLabel = [];
    ax_density_sim(d).XTick = [];
    ax_density_sim(d).YTick = [-1,0,1];
    ax_density_sim(d).FontUnits = 'centimeters';
    ax_density_sim(d).FontSize = fontsize_den;
    
    
    %---------------- PLOT SIM LINEOUT ------------------------------------
    
    axlongprofile_sim(d) = nexttile;
    
    % manual triming of experimental data
    ind_r = (r_plot < sigma_rexp*10);
    density_trim = density_plot(ind_r,:);
    r_plot = r_plot(ind_r);
    
    %end of manual triming
    
    long_profile = sum(density_trim);
    plongprofile_sim = plot(z_plot,long_profile,'k');
    axlongprofile_sim(d).XDir = 'reverse';
    xlim(xi_range); ylim([0 2.9e7])
    
    if d == 1
        ylabel({'density','(arb. units)'},'FontUnits','centimeters','FontSize',fontsize_den)
    end
    
    if d ~= 3
        axlongprofile_sim(d).XTick = [];
        axlongprofile_sim(d).XTickLabel = [];
    end
    axlongprofile_sim(d).YTickLabel = [];
    axlongprofile_sim(d).YTick = [];
end

xlabel(tile_handle,'delay (ps)','FontSize',9); % 1 point = 1/72 inches = 0.0353 cm; 9 point = 0.388
drawnow;
pause(1);

P.fig_handle = fig_handle;
P.plot_name = 'allprofiles';
P.save_plot();
