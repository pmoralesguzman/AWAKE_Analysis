
%________________________________________________________________________
% Plot the 2D wakefields together with the proton bunch, or each
% individually
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 16/03/2021
%________________________________________________________________________

% close all;

% file location variables
datadir = 'gm20';%'gradsim_shortbunch';
dataformat = 'mat';
useAvg = 0;
dump_list = 45:1:55;

% saving data
save_flag = 1;
save_format = {'png'};

% plasma properties
plasmaden = 1.81e14;

% choose property to plot
property_plot = 'both'; % density, wakefields, both
include_lineout = 'density_profile'; % 'no','both','field_lineout','density_profile'

% choose fields to plot
wakefields_direction = 'long'; % trans, long
lineout_point = 2;

% choose species density to plot
species = 'proton_beam';    

% choose limits (in cm, must denormalize)
trans_range = [0 0.24];
% xi_range = [4.8 4.4];
xi_range = [14 0.38];
plot_density_lims = [-inf inf];

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


% directory to save the plots
plots_dir = ['field_density_longprofile/',datadir,'/',...
    property_plot,'/',wakefields_direction,'/',...
    'xi',num2str(round(xi_range(1))),'xi',...
    num2str(round(xi_range(2)))];

P = Plotty('datadir',datadir,'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list,...
    'plot_scale',plot_scale,...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'lineout_point',lineout_point,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'plot_density_lims',plot_density_lims,...
    'make_pause',make_pause,'fig_number',fig_number);

figure(fig_number);
P.plot_field_density('trans_lines_position',[-0.01,0.01],'include_lineout',include_lineout);

% xline(0.1,'LineWidth',2)
hold on
% scatter((P.dtime+P.simulation_window-locs),(pks),50,'r','filled')
hold off

P.plot_name = 'FreqSmoothZoom';
P.save_plot();

