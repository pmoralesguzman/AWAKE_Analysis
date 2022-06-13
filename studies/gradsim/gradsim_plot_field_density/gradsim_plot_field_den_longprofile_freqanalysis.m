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
% Last update: 05/08/2020
%________________________________________________________________________

clear;
close all;

% file location variables
datadir = 'g0';
dataformat = 'h5';
useAvg = true;
dump_list = 3:1:3;

% saving data
save_flag = true;
save_format = {'png'};

% plasma properties
plasmaden = 1.81e14;

% choose fields to plot
wakefields_direction = 'long'; % trans, long
lineout_point = 3;

% choose species density to plot
species = 'proton_beam'; 

% choose limits (in cm, must denormalize)
trans_range = [0 0.1];
% xi_range = [4.8 4.4];
xi_range = [14.2 0.0];

% choose property to plot
property_plot = 'both'; % density, wakefields, both

% create movie or not
create_movie = true;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose scaling
plot_scale = 'linear';

% choose if make pause or not
make_pause = false;

% figure number 
fig_number = 3;


% directory to save the plots
plots_dir = ['gradsim_test/2D/',datadir,'/',...
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
    'make_pause',make_pause,'fig_number',fig_number);


figure(fig_number);
P.plot_field_density('include_lineout','both','trans_lines_position',[-0.0536 0.0536]);

% P.plot_name = [datadir,property_plot,n];
% P.save_plot();

