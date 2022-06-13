
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
datadir = 'gp5'; %'gradsim_shortbunch';
dataformat = 'h5';
useAvg = 1;
dump_list = 100:1:100;

% saving data
save_flag = 1;
save_format = {'png'};

% plasma properties
plasmaden = 1.81e14;

% choose property to plot
property_plot = 'both'; % density, wakefields, both
include_lineout = 'no'; % 'no','both','field_lineout','density_profile'

% choose fields to plot
wakefields_direction = 'long'; % trans, long
lineout_point = 13;

% choose species density to plot
species = 'proton_beam';    

% choose limits (in cm, must denormalize)
trans_range = [0 0.16];
% xi_range = [10 7];
xi_range = [10 0];
plot_density_lims = [0 2e10];
plot_field_lims = [-400 400];

% create movie or not
create_movie = 1;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose scaling
plot_scale = 'linear';

% choose if make pause or not
make_pause = false;

% figure number 
fig_number = 3;


% directory to save the plots
plots_dir = ['tCSC/field_density/',datadir];

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

figure(fig_number);
P.plot_field_density('ylinepos',[-0.018,0.018],'ylineflag',1,...
    'include_lineout',include_lineout);




