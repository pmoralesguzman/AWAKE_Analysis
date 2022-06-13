%________________________________________________________________________
% Script to calculate, from the waterfall data, the dephasing of the fields
% close to some selected xi.
% Especially developed for the seeding position studies.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/07/2020
%________________________________________________________________________

clear;
close all;

load('color_red_to_blue.mat');

datadir = 'gp5';

plasmaden = 1.81e14;
property = 'fields';
dump_list = 1:1:100;
useAvg = false;
dataformat = 'mat';

wakefields_direction = 'long';
lineout_point = 5;

trans_range = [0 0.02]; % if density
dephasing_xis = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5]; % cm
label_leg = num2cell(dephasing_xis);
dephasing_search = '0x'; % 0x, max
force_waterfall = false;

% plot parameters
save_flag = 0;
save_format = {'png','eps'};
plots_dir = 'gradsim/dephasing/';
legs = {'1','2','3','4','5'};
x_range = [-inf inf]; %OPA.dephasing_xi + [0 21];

load_microbunches = 0;

% initialize some plotty
P = Plotty('plasmaden',plasmaden,...
    'save_flag',save_flag,...
    'plots_dir',plots_dir,'save_format',save_format,...
    'wakefields_direction',wakefields_direction);

OPA = OsirisPhaseAnalysis('datadir',datadir,...
    'property',property,'species','proton_beam',...
    'wakefields_direction',wakefields_direction,...
    'trans_range',trans_range,...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,...
    'dephasing_search',dephasing_search,...
    'force_waterfall',force_waterfall,...
    'lineout_point',lineout_point);

for ph = 1:length(dephasing_xis)
    OPA.dephasing_xi = dephasing_xis(ph);
    OPA.dephasing_first = [];
    OPA.dephasing();
    
    % waterfall plot datas
    P.waterfall_xi = OPA.waterfall_xi;
    P.waterfall_z = OPA.waterfall_z;
    P.waterfall_mat = OPA.waterfall_mat;
    P.property = OPA.property;
    
    P.waterfall_plot();
    xlim(x_range);
        
    % add the dephasing line (zero-crossing or max) to the waterfall plot
    hold on
    plot((-OPA.dephasing_line)*OPA.plasma_wavelength + OPA.simulation_window - OPA.dephasing_first,...
        linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
        'LineWidth',2,'color',[0.5,0.5,0.5])
    % load line from the microbunches for double plotting
    if load_microbunches
        phase_filename = ['save_files/phases/phase_',datadir,'_',num2str(OPA.dephasing_xi),'.mat'];
        bunches = load(phase_filename);
        bunches_x = bunches.phase_x_plot;
        bunches_y = bunches.phase_y_plot;
        plot(bunches_x,bunches_y,'LineWidth',2,'color','k','LineStyle',':')
    end
    
    title(['\xi = ',num2str(OPA.dephasing_xi),' cm behind seed. pos. (',datadir,')'])
    hold off
    
    drawnow;
    P.plot_name = ['gradsim_field_',datadir,'_',num2str(OPA.dephasing_xi)]; % PLOT NAME
    P.save_plot();
    dephasing_lines{ph} = OPA.dephasing_line;
end

% just the dephasing line, holds the figure for the next datadir
fig2 = figure(100);
colororder(ccrb);
for xi = 1:length(dephasing_xis)
    hold on
    p1 = plot(linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
        dephasing_lines{xi},'LineWidth',2);
    hold off
end % for xi

% labels for the dephasing lines
ylabel('phase (\lambda_p)')
xlabel('z (m)')
legend(legs{:},'Location','best')
xlim([0 10])
%     ylim([-6 1]);
title(['\xi = ',num2str(OPA.dephasing_xi),' cm behind seed. pos.'])

P.plot_name = ['gradsim_dephasing_fields_',datadir,num2str(OPA.dephasing_xi)]; % PLOT NAME
P.fig_handle = fig2;
P.save_plot();

