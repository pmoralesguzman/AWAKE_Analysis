%---------------------------------------------------------------------
% Plot waterfall, energy gain from Ez along xi (point-like bunch),
% maximum energy gain and position, and phase evolution.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/05/2021
%---------------------------------------------------------------------
clear;
close all;

load('color_red_to_blue.mat');

%% input parameters

gpm = 'a';

switch gpm
    case 'p'
        d_start = 5;
        datadirs = {'g0','gp5','gp10','gp15','gp20'};
        leg = {'0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
    case 'a'
        d_start = 1;
        datadirs = {'gm20','gm15','gm10','gm5',...
            'g0','gp5','gp10','gp15','gp20'};
        leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
end
ccrb = ccrb(d_start:9,:);

plasmaden = 1.81e14;
dump_list = 0:1:100;
useAvg = 0;
dataformat = 'mat';

xi_range = [5.7 6.3];
lineout_point = 5;
plots_dir = ['gradsim_paper/ene/wz_ene'];


OPA = OsirisPhaseAnalysis('property','fields','wakefields_direction','long',...
    'plasmaden',plasmaden,'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,'force_waterfall',false,...
    'lineout_point',lineout_point);
P = Plotty('plasmaden',plasmaden,'plots_dir',plots_dir);


for d = 1:length(datadirs)
    OPA.datadir = datadirs{d};
    
    % waterfall plot
    OPA.build_waterfall();
    ind_xi = fliplr((OPA.waterfall_xi < xi_range(2)) & (OPA.waterfall_xi > xi_range(1)));
    waterfall_mat = OPA.waterfall_mat(:,ind_xi);
    waterfall_z = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),size(OPA.waterfall_mat,1));
    
    % field integral
    field_integral = cumtrapz(waterfall_z,-(waterfall_mat));
    
    %% max energy gain
    [max_energygain,pos_maxenergygain] = max(field_integral);
    
    [maxmax_energygain,best_position] = max(max_energygain);
    energy_gain(d,:) = field_integral(:,best_position);
    field_amplitude(d,:) = -flip(waterfall_mat(:,best_position));
    max_field_amplitude(d,:) = max(flip(waterfall_mat),[],2);
    
end

% energy gain
fig_eg = figure(1);
colororder(ccrb);
line_style = {':','--','-.','-','-','-','-.','--',':'};
hold on
for d = 1:length(datadirs)
    plot(waterfall_z,energy_gain(d,:),line_style{d + d_start - 1},'LineWidth',2);
end
hold off

xlim([0 10])
ylim([-80 2000])

xlabel('distance from plasma end (m)')
ylabel ('maximum energy gain (MeV)')

legend(leg,'Location','best')
P.fig_handle = fig_eg;
switch gpm
    case 'p'
        P.plot_name = 'wz_int_positive';
    case 'a'
        P.plot_name = 'wz_int';
end
P.save_plot();

% field amplitude
fig_fa = figure(2);
colororder(ccrb);
hold on
for d = 1:length(datadirs)
    plot(waterfall_z,max_field_amplitude(d,:),line_style{d + d_start - 1},'LineWidth',2);
end
hold off
xlim([0 10])

xlabel('z(m)')
ylabel ('E_z (MV)')

% legend(leg,'Location','best')
P.fig_handle = fig_fa;
switch gpm
    case 'p'
        P.plot_name = 'field_positive';
    case 'a'
        P.plot_name = 'field';
end

P.save_plot();

