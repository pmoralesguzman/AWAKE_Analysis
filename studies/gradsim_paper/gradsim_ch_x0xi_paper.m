%________________________________________________________________________
% gradsim paper
% Script to follow the zero crossing of fields from waterfall plots. It
% also creates the waterfall plots if not yet made.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 18/02/2021
%________________________________________________________________________
clear;
close all;

plots_dir = ['gradsim_paper/ch/x0xi/'];
plot_name = ['dephasing','1471','allgrads0'];

% load color order for 9 gradients
load('color_purple_to_green.mat'); % cc
load('color_red_to_blue.mat'); % ccrb

% cell plotting parameters
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
line_style = {':','--','-.','-','-','-','-.','--',':'};

% plotting parameters
fontsize_annotation = 12;
fontsize_label = 14;

% study parameters
plasmaden = 1.81e14;
dump_list = 0:1:100;
useAvg = false;
dataformat = 'mat';
dephasing_xi = [14,7,1]; % cm

% initialize classes
P = Plotty('plasmaden',1.81e14,'plots_dir',plots_dir,'plot_name',plot_name);
OPA = OsirisPhaseAnalysis('datadir',datadirs{1},...
    'property','fields','wakefields_direction','long',...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,'dephasing_xi',dephasing_xi(1),...
    'force_waterfall',false);

OPA.dephasing();

% theoretical phase
T = linspace(0.01/P.c_m,10/P.c_m,1000);
dT = T(2) - T(1); 
z = T*100*P.c_m;
vph_xi1 = P.c_cm - 0.5.*(0.25./z).^(1/3).*(0.0376649*P.e_mass_kg/(2*P.p_mass_kg*426.439))^(1/3).*P.c_cm;
ph = 0;
for t = 1:length(T)
    ph(t+1) = ph(t) + dT*vph_xi1(t);
end
ph(end) = [];
% figure
% plot(z,vph_xi1);
% figure
% plot(z(2:end),(ph(2:end) - P.c_cm*T(1:end-1))/0.25);

% initialize variables
dephasing_z = zeros(length(dephasing_xi),length(datadirs),length(OPA.dephasing_line));
dephasing_lines = zeros(length(dephasing_xi),length(datadirs),length(OPA.dephasing_line));

%% start the script
fig_dephase = figure(1);
fig_dephase.OuterPosition = [100 100 1200 400];
colororder(ccrb);
tt = tiledlayout(1,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

for xi = 1:length(dephasing_xi)
    
    for d = 1:length(datadirs)
        
        OPA.datadir = datadirs{d};
        OPA.dephasing_xi = dephasing_xi(xi);
        OPA.dephasing_first = [];
        
        OPA.dephasing();
        
        dephasing_z(xi,d,:) = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),...
            length(OPA.dephasing_line));
        dephasing_lines(xi,d,:) = OPA.dephasing_line*2; % HALF PLASMA WAVELENGTH !!!
        OPA.progress_dump('dephasing lines',d,length(datadirs))
        
    end % datadirs
    
    ax_x0(xi) = nexttile;
    ax_x0(xi).FontSize = fontsize_label;
    hold on
    
    
    
    for d = 1:length(datadirs)
        p1 = plot(squeeze(dephasing_z(xi,d,:)),squeeze(dephasing_lines(xi,d,:)),line_style{d},'LineWidth',2);
%         plot(z(2:end)/100,(ph(2:end) - P.c_cm*T(1:end-1))/(0.5*P.plasma_wavelength));
    end % datadirs
    yline(-1,'--','LineWidth',1,'color',[0 0.4470 0.7410])
    yline(-2,'--','LineWidth',1,'color',[0 0.4470 0.7410])
    %     yline(0.0,'--','LineWidth',1,'color',[0 0.4470 0.7410])
    xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
    hold off
    
    switch dephasing_xi(xi)
        case 14
            position_word = '(back of the bunch)';
        case 7
            position_word = '(middle of the bunch)';
        case 1
            position_word = '(front of the bunch)';
    end % dephasing xi (title)
    title(['\xi_0 = ',num2str(dephasing_xi(xi)),' cm ',position_word]);
    
    xlim([0 10])
    
    if OPA.dephasing_first == 40
        ylim([-6.5 1.5]*2)
    else
        ylim([-6.5 0.5]*2)
    end % if OPA dephasing first
    
end % xi

legend(ax_x0(3),leg,'location','southeast','FontSize',fontsize_annotation,'NumColumns',2)
ylabel(tt,'zero-crossing position shift (\lambda_p/2)','FontSize',fontsize_label)
xlabel(tt,'z (m)','FontSize',fontsize_label)

% Create textarrow
annotation(fig_dephase,'textarrow',[0.159431137724551 0.221556886227545],...
    [0.6 0.787412587412587],'Color',[33,113,181]/256,...
    'String',{'positive','gradient values'},...
    'Interpreter','latex',...
    'HeadStyle','none','FontSize',fontsize_annotation);

% Create textarrow
annotation(fig_dephase,'textarrow',[0.264221556886228 0.306137724550898],...
    [0.300699300699301 0.488111888111888],...
    'Color',[203,24,29]/256,...
    'String',{'negative','gradient values'},...
    'Interpreter','latex',...
    'HeadStyle','none','FontSize',fontsize_annotation);

% Create textarrow
annotation(fig_dephase,'textarrow',[0.166167664670659 0.25374251497006],...
    [0.404195804195804 0.692307692307692],'String','g = 0 \%/m',...
    'Interpreter','latex',...
    'HeadStyle','none','FontSize',fontsize_annotation);

P.fig_handle = fig_dephase;
P.save_plot();

save('loading_files/gradsim_cache/gradsim_dephasing.mat','dephasing_z','dephasing_lines');


