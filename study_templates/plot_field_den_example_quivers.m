
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
% Last update: 06/07/2022
%________________________________________________________________________

% close all;
clear;

% file location variables
datadir = 'r2o_2'; 
extradatadir = ''; 

dataformat = 'h5';
useAvg = 0;
dump_list = [1:5:100];


% saving plot
save_flag = 0;
save_format = {'png'}; % fig, pdf

% plasma properties
plasmaden = 2e14; % ! 2e14, 7e14 

% choose property to plot
property_plot = {'density','wakefields'}; % density, wakefields, both
species_to_plot = {'proton_beam'}; %proton_beam, electron_seed, electrons
include_phasespace = 0; % 0,1
field_geometry = 'cartesian'; % cartesian, cylindrical
on_axis = 'lineout'; % how to get the density profile
% int = charge, sum = density, lineout = density

% choose fields to plot
wakefields_direction = 'trans'; % trans, long
lineout_point = '0.02'; % cm

% choose limits (in cm, must denormalize)
trans_range = [0 0.11]; % [smaller,larger]
xi_range = [2 0]; % [larger, smaller]
plot_density_lims = [-inf inf]; % [smaller,larger] 
plot_field_lims = [-inf inf]; % [smaller,larger]

% create movie or not
create_movie = 1;

% choose between normalized and denormalized units
denormalize_flag = 1; % true, false

% choose if make pause or not
make_pause = 1;

% figure number
fig_number = 1;

% directory to save the plots
plots_dir = ['example/fielden/',datadir,'x'];

code = 'osiris';

include_quivers = 1;

switch datadir

    otherwise
        extitles{1} = '';
        exlegends{1} = '';

end

switch extradatadir

    otherwise
        extitles{2} = '';
        exlegends{2} = '';
end


P = Plotty('datadir',datadir,'extradatadir',extradatadir,'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list,...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,...
    'lineout_point',lineout_point,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'plot_density_lims',plot_density_lims,'plot_field_lims',plot_field_lims,...
    'make_pause',make_pause,'fig_number',fig_number,...
    'extitles',extitles,'exlegends',exlegends,...
    'field_geometry',field_geometry,'on_axis',on_axis,...
    'species_to_plot',species_to_plot,'code',code);%,...
    %'plot_name',datadir);

figure(fig_number);
P.plot_field_density('ylinepos',[-0.02,0.02],'yline_flag',1,...
    'xlinepos',[0.423726],'xline_flag',0,...
    'include_density_profile',1,'include_field_lineout',1,...
    'include_phasespace_profile',0,'include_quivers',include_quivers); %#ok<NBRAK> 




