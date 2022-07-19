
%________________________________________________________________________
% Plot transvere fields in 2D and lineout
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress 
%
% P. I. Morales Guzman
% Last update: 19/07/2022
%________________________________________________________________________

% close all;
clear;

% file location variables
datadir = 'fdr_26e'; 
extradatadir = ''; 
movie_dir = 'fdr_26e';

dataformat = 'mat';
useAvg = 0;
dump_list = [0:1:0];


% saving plot
save_flag = 0;
save_format = {'png'}; % fig, pdf

% plasma properties
plasmaden = 2e14; % ! 2e14, 7e14 

% choose property to plot
property_plot = {'wakefields'}; % density, wakefields, both
species_to_plot = {'electrons','electron_seed'}; %proton_beam, electron_seed, electrons
include_phasespace = 0; % 0,1
field_geometry = 1; % cartesian, cylindrical
integration_type = 'sum'; % how to get the density profile
% trapz,simpsons,sum = charge, sum_density = density, just_lineout = density

% choose fields to plot
wakefields_direction = ''; % trans, long

lineout_point = '0.02'; % cm (fie the field lineout)

% choose limits (in cm, must denormalize)
trans_range = [0 0.11]; % [smaller,larger]
xi_range = [10 0]; % [larger, smaller]
ylinepos = [-0.02,0.02]; % [lower,upper] also limits for the charge/density longitudinal profile
yline_flag = 1;
xlinepos = 0;
xline_flag = 0;
plot_density_lims = [-inf inf]; % [smaller,larger] 
plot_field_lims = [-inf inf]; % [smaller,larger]

% includes
include_density_profile = 0;
include_field_lineout = 1;
include_phasespace_profile = 0;

% create movie or not
create_movie = 0;

% choose between normalized and denormalized units
denormalize_flag = 1; % true, false

% choose if make pause or not
make_pause = 0;

% figure number
fig_number = 1;

% directory to save the plots
plots_dir = ['example/fielden/',datadir,'x'];

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
    'movie_dir',movie_dir,'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,...
    'lineout_point',lineout_point,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...    
    'plot_density_lims',plot_density_lims,'plot_field_lims',plot_field_lims,...
    'make_pause',make_pause,'fig_number',fig_number,...
    'extitles',extitles,'exlegends',exlegends,...
    'integration_type',integration_type,...
    'species_to_plot',species_to_plot,'direction','r');%,...
    %'plot_name',datadir);

P.plot_field_density('ylinepos',ylinepos,'yline_flag',yline_flag,...
    'xlinepos',xlinepos,'xline_flag',xline_flag,...
    'include_density_profile',include_density_profile,...
    'include_field_lineout',include_field_lineout,...
    'include_phasespace_profile',include_phasespace_profile); 




