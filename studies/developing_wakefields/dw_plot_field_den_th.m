
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
clear;

% file location variables
datadir = 'protonbm'; %
extradatadir = ''; %

dataformat = 'mat';
useAvg = 1;
dump_list = 0:0;

% saving data
save_flag = 1;
save_format = {'png'};

% plasma properties
plasmaden = 2e14; % !!!!!!

% choose property to plot
property_plot = 'both'; % density, wakefields, both
include_lineout = 'both'; % 'no','both','field_lineout','density_profile'

% choose fields to plot
wakefields_direction = 'trans'; % trans, long
lineout_point = '0.02';

% choose species density to plot
species = 'proton_beam';    

% choose limits (in cm, must denormalize)
trans_range = [0 0.11];
% xi_range = [9.5 6.5];
xi_range = [19 0];
% plot_density_lims = [0 2.5e11];
% plot_field_lims = [-40 40];

plot_density_lims = [-inf inf];
plot_field_lims = [-inf inf];

% create movie or not
create_movie = 0;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose scaling
plot_scale = 'linear';

% choose if make pause or not
make_pause = false;

% figure number 
fig_number = 1;


% directory to save the plots
plots_dir = ['DW/fielden/',datadir,''];

switch datadir
    case 'DWdc3_lcode_x1'
        extitles{1} = '100 % density';
        exlegends{1} = '100 % density';
    case 'DWdc3_lcode_x666'
        extitles{1} = '15 % density';
        exlegends{1} = '15 % density';
    case 'DWdc3_lcode_x666_pi'
        extitles{1} = '15 % density, pi';
        exlegends{1} = '15 % density, pi';
    case 'DWdc_lcode_x666'
        extitles{1} = '15 % density, no den. feat.';
        exlegends{1} = '15 % density, no den. feat.';
    case 'electron_bunch'
        extitles{1} = 'electron bunch 550 pC';
        exlegends{1} = 'electron bunch 550 pC';
    case 'proton_electron'
        extitles{1} = 'electron bunch 550 pC';
        exlegends{1} = 'electron bunch 550 pC';

end

switch extradatadir
    case 'DWdc3_lcode_x1'
        extitles{2} = '100 % density';
        exlegends{2} = '100 % density';
    case 'DW_lcode_nofront'
        extitles{2} = 'no front';
        exlegends{2} = 'no front';
    case 'DWdc3_lcode_x666'
        extitles{2} = '15 % density';
        exlegends{2} = '15 % density';
    case 'DWdc_lcode_x666_pi'
        extitles{2} = '15 % density, no den. feat., pi';
        exlegends{2} = '15 % density, no den. feat., pi';
    otherwise
        warning('extradir title otherwise')
        extitles{2} = '';
        exlegends{2} = '';
end



P = Plotty('datadir',datadir,'extradatadir',extradatadir,'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list,...
    'plot_scale',plot_scale,...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'lineout_point',lineout_point,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'plot_density_lims',plot_density_lims,'plot_field_lims',plot_field_lims,...
    'make_pause',make_pause,'fig_number',fig_number,...
    'extitles',extitles,'exlegends',exlegends);%,...
    %'plot_name',datadir);

figure(fig_number);
P.plot_field_density('ylinepos',[-0.02,0.02],'ylineflag',1,...
    'xlinepos',[0],'xlineflag',1,...
    'include_lineout',include_lineout); %#ok<NBRAK> 




