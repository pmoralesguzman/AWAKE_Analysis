
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
datadir = 'fdr_26'; 
extradatadir = ''; 
movie_dir = 'example_density';

dataformat = 'mat';
useAvg = 0;
dump_list = [0:1:200];


% saving plot
save_flag = 1;
save_format = {'png'}; % fig, pdf

% plasma properties
plasmaden = 2e14; % ! 2e14, 7e14 

% choose property to plot
property_plot = {'density'}; % density, wakefields, both
species_to_plot = {'electrons','electron_seed'}; %proton_beam, electron_seed, electrons
transparency_flag = 0;
include_phasespace = 0; % 0,1
field_geometry = 'cartesian'; % cartesian, cylindrical
integration_type = 'just_lineout'; % how to get the density profile
% trapz,simpsons,sum = charge, sum_density = density, just_lineout = density

% choose fields to plot
wakefields_direction = 'trans'; % trans, long
lineout_point = '0.02'; % cm (fie the field lineout)

% choose limits (in cm, must denormalize)
trans_range = [0 0.09]; % [smaller,larger]
xi_range = [10 0]; % [larger, smaller]
ylinepos = [0,0.005,0.01,0.015,0.02]; % [lower,upper] also limits for the charge/density longitudinal profile
yline_flag = 0;
xlinepos = 0;
xline_flag = 0;
plot_density_lims = [-inf inf]; % [smaller,larger] 
plot_field_lims = [-inf inf]; % [smaller,larger]

% includes
include_density_profile = 1;
include_field_lineout = 0;
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
    'create_movie',create_movie,'movie_dir',movie_dir,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'field_geometry',field_geometry,'integration_type',integration_type,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,...
    'lineout_point',lineout_point,...
    'plot_density_lims',plot_density_lims,'plot_field_lims',plot_field_lims,...
    'make_pause',make_pause,'fig_number',fig_number,...
    'extitles',extitles,'exlegends',exlegends,...
    'species_to_plot',species_to_plot);%,...
    %'plot_name',datadir);

P.plot_field_density('ylinepos',ylinepos,'yline_flag',yline_flag,...
    'xlinepos',xlinepos,'xline_flag',xline_flag,...
    'include_density_profile',include_density_profile,...
    'include_field_lineout',include_field_lineout,...
    'include_phasespace_profile',include_phasespace_profile,...
    'transparency_flag',transparency_flag); 




