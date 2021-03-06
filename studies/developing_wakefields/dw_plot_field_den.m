
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
% datadir = 'DW_lcode_x5_pi_rz2_th_m2_w3'; %''; DW_lcode_justelectron DW_justelectron
% extradatadir = 'DW_lcode_x5_pi_rz2_th_m2'; %'DWdc3_lcode_x1_pi'; DW_lcode_nofront

datadir = 'DW_lcode_x20_pi'; %''; DW_lcode_justelectron DW_justelectron
extradatadir = ''; %'DWdc3_lcode_x1_pi'; DW_lcode_nofront

dataformat = 'mat';
useAvg = 0;
dump_list = 0:2:100;

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
species = 'proton_beam';    %proton_beam  electron_seed

% choose limits (in cm, must denormalize)
trans_range = [0 0.12];
xi_range = [15 0];
plot_density_lims = [0 15e9];
plot_field_lims = [-210 210];
% plot_field_lims = [-inf inf];

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
% plots_dir = ['DW/fielden/',datadir,'x'];
plots_dir = ['DW/fielden/','x20_pi_p2'];

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
        extitles{1} = '15 % density';
        exlegends{1} = '15 % density';
    case 'DWdc_lcode_nofront_nodf'
        extitles{1} = 'baseline';
        exlegends{1} = 'baseline';
    case 'DW_justelectron'
        extitles{1} = 'osiris';
        exlegends{1} = 'osiris';
    case {'DW_lcode_x5','DW_lcode_x5_pi'}
        extitles{1} = 'baseline';
        exlegends{1} = 'baseline';
    case 'DW_lcode_x5_th'
        extitles{1} = 'time HD';
        exlegends{1} = 'time HD';
    case 'DW_lcode_x20'
        extitles{1} = 'front 20\%';
        exlegends{1} = 'front 20\%';
    otherwise
        extitles{1} = '';
        exlegends{1} = '';

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
    case 'DWdc_lcode_x666'
        extitles{2} = '15 % density';
        exlegends{2} = '15 % density';
    case 'DWdc_lcode_x666_pi'
        extitles{2} = '15 % density, pi';
        exlegends{2} = '15 % density, pi';
    case 'DW_lcode_jet'
        extitles{2} = 'jet';
        exlegends{2} = 'jet';
    case 'DW_lcode_x20_pi'
        extitles{2} = 'front 20\% dephased $\pi$';
        exlegends{2} = 'front 20\% dephased $\pi$';
    otherwise
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
    'extitles',extitles,'exlegends',exlegends,'include_phasespace',1);%,...
    %'plot_name',datadir);

figure(fig_number);
P.plot_field_density('ylinepos',[-0.02,0.02],'ylineflag',1,...
    'xlinepos',[0.423726],'xlineflag',0,...
    'include_lineout',include_lineout); %#ok<NBRAK> 




