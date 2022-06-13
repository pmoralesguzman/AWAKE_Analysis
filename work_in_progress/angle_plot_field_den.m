
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

datadir = 'r2l_302_c_noe'; %''; DW_lcode_justelectron DW_justelectron
extradatadir = ''; %'DWdc3_lcode_x1_pi'; DW_lcode_nofront

dataformat = 'mat';
useAvg = 0;
dump_list = 0:1:0;
% dump_list = 0:1:10;

% saving data
save_flag = 1;
save_format = {'png'};

% plasma properties
plasmaden = 2e14; % !!!!!!

% choose property to plot
property_plot = 'both'; % density, wakefields, both
include_lineout = 'both'; % 'no','both','field_lineout','density_profile'
include_phasespace = 1;

% choose fields to plot
wakefields_direction = 'trans'; % trans, long
lineout_point = '0.02';

% choose species density to plot
species = 'proton_beam';    %proton_beam  electron_seed

% choose limits (in cm, must denormalize)
trans_range = [0 0.11];
xi_range = [30 0];
plot_density_lims = [0 inf];
plot_field_lims = [-inf inf];

% create movie or not
create_movie = 1;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose scaling
plot_scale = 'linear';

% choose if make pause or not
make_pause = 1;

% figure number
fig_number = 1;


% directory to save the plots
plots_dir = ['angle/fielden/',datadir,''];
% plots_dir = ['r2plot/fielden_summary/','302ce550_b_vs_pi'];

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
    'plot_scale',plot_scale,...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'lineout_point',lineout_point,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'plot_density_lims',plot_density_lims,'plot_field_lims',plot_field_lims,...
    'make_pause',make_pause,'fig_number',fig_number,...
    'extitles',extitles,'exlegends',exlegends,'include_phasespace',include_phasespace);%,...
    %'plot_name',datadir);

figure(fig_number);
P.plot_field_density('ylinepos',[-0.02,0.02],'ylineflag',1,...
    'xlinepos',[0.423726],'xlineflag',0,...
    'include_lineout',include_lineout,...
    'include_phasespace_profile',1,...
    'include_density_profile',0,'include_field_lineout',1); %#ok<NBRAK> 




