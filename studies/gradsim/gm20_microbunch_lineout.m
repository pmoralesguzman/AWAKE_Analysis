
%________________________________________________________________________
% Plot the lineouts to show that the microbunches come from defocused
% charge
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 03/05/2021
%________________________________________________________________________

% close all;

% file location variables
datadir = 'gm20'; %'gradsim_shortbunch';
dataformat = 'h5';
useAvg = 1;
dump_list = 38:1:38;

% saving data
save_flag = 1;
save_format = {'png'};

% plasma properties
plasmaden = 1.81e14;

% choose property to plot
property_plot = 'density'; % density, wakefields, both
include_lineout = 'both'; % 'no','both','field_lineout','density_profile'

% choose fields to plot
wakefields_direction = 'trans'; % trans, long
lineout_point = 13;

% choose species density to plot
species = 'proton_beam';    

% choose limits (in cm, must denormalize)
trans_range = [0 0.018];
% xi_range = [10 7];
xi_range = [2 0];
plot_density_lims = [0 6e11];
plot_field_lims = [-100 100];

% create movie or not
create_movie = 0;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose scaling
plot_scale = 'linear';

% choose if make pause or not
make_pause = false;

% figure number 
fig_number = 3;


% directory to save the plots
plots_dir = ['gradsim/field_density/bunch_moving/',datadir];

P = Plotty('datadir',datadir,'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list,...
    'plot_scale',plot_scale,...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'lineout_point',lineout_point,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'plot_density_lims',plot_density_lims,'plot_field_lims',plot_field_lims,...
    'make_pause',make_pause,'fig_number',fig_number);

P.plot_lineout();
P.plot_handle.LineWidth = 2;

P.dump_list = 60;
hold on
P.plot_lineout();
hold off

P.plot_handle.Color = [0.03,0.6641,0.0781];
P.plot_handle.LineWidth = 2;

P.dump_list = 80;
hold on
P.plot_lineout();
hold off

P.plot_handle.Color = [0.0641,0.0039,0.6781];
P.plot_handle.LineWidth = 2;


P.dump_list = 100;
hold on
P.plot_lineout();
hold off

P.plot_handle.Color = [0.6641,0.0039,0.0781];
P.plot_handle.LineWidth = 2;

legend('z = 3.8 m','z = 6 m','z = 8 m','z = 10 m')

P.plot_name = 'lineout_normalized_4';
P.save_plot();

