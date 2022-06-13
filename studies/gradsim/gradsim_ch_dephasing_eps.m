%________________________________________________________________________
% gradsim paper
% Script joining the dephasing, charge and mean transverse field amplitude
% for a 3x3 tile plot with all info. It does not calculate anything, bith
% rather takes the data from
% - gradsim_ch_mean_trans_fields_in_xi
% - gradsim_ch_charge_fraction_in_xi
% - gradsim_ch_x0shift_in_xi
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 18/02/2020
%________________________________________________________________________

clear;
% close all;

plots_dir = ['gradsim_eps/dephasing'];
plot_name = ['dephasing','xi7','allgrads'];

% load color order for 9 gradients
load('color_red_to_blue.mat'); % loaded vars: ccrb

% load x0 shift data
load('loading_files/gradsim_cache/gradsim_dephasing.mat'); % loaded vars: dephasing_z, dephasing_lines

% plotting parameters
fontsize_annotation = 9; % points (1 point = 1/72 inches = 0.0353 cm; 9 point = 0.388)
fontsize_label = 0.4; % cm
dephasing_first = 0;
letterbox_x = 0.34;
letterbox_y1 = -6.8;
letterbox_y2 = 82.7;
letterbox_y3 = 1.37;


% cell plotting parameters
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
leg = {'$-2$\,\%/m','$-1.5$\,\%/m','$-1$\,\%/m','$-0.5$\,\%/m',...
    '\ $0$\,\%/m','$+0.5$\,\%/m','$+1$\,\%/m','$+1.5$\,\%/m','$+2$\,\%/m'};
letters = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
line_style = {':','--','-.','-','-','-','-.','--',':'};

% study parameters
dephasing_xi = [14,7,1]; % cm

% initialize counters
i_letter = 1;

% initialize classes
P = Plotty('plasmaden',1.81e14,'plots_dir',plots_dir,'plot_name',plot_name);

%% start the plotting, with dephasing
fig_join = figure(2);
% fig_join.Units = 'normalized';
% fig_join.OuterPosition = [100 100 1200 800];
fig_join.Units = 'centimeters';
fig_join.Position = [1,1,8.6,8.6*3/4]*1.5;
colororder(ccrb);
% tt = tiledlayout(1,1);
% tt.TileSpacing = 'compact';
% tt.Padding = 'compact';

% dephasing (x0 shift)

for xi = 2:2
    
%     ax_x0(xi) = nexttile;
    ax_x0(xi) = gca;
    ax_x0(xi).FontUnits = 'centimeters';
    ax_x0(xi).FontSize = fontsize_label;
    
    hold on
    for d = 1:length(datadirs)
        p1 = plot(squeeze(dephasing_z(xi,d,:)),squeeze(dephasing_lines(xi,d,:)),...
            line_style{d},'LineWidth',2);
    end
    
    yline(-1,'--','LineWidth',1,'color',[0 0 0])
    yline(-2,'--','LineWidth',1,'color',[0 0 0])
    if xi < 3
        xline(4,'--','LineWidth',1,'color',[0 0 0]);
    else
        plot([4 4],[1 -2],'--','LineWidth',1,'color',[0 0 0])
    end
    hold off

    
    xlim([0 10])
    ylim([-4 0.5]*2)
    
end % xi

ylabel('x_0 (\lambda_{pe0}/2)')

% Create textarrow
annotation(fig_join,'textarrow',[0.636956092492287 0.693171188026193],...
    [0.385324980194506 0.489795918367347],...
    'Color',[0.79296875 0.09375 0.11328125],...
    'String','$g < 0$\,\%/m',...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',12);

% Create textarrow
annotation(fig_join,'textarrow',[0.42843779232928 0.604303086997194],...
    [0.445839874411303 0.665620094191523],'String','$g = 0$\,\%/m',...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',12);

% Create textarrow
annotation(fig_join,'textarrow',[0.31711880261927 0.481758652946679],...
    [0.58712715855573 0.759811616954474],...
    'Color',[0.12890625 0.44140625 0.70703125],...
    'String','$g > 0$\,\%/m',...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',12);

xlabel('z (m)')

P.fig_handle = fig_join;
P.save_plot();









