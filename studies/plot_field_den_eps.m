
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

datadir = 'r2l_302_c_e550_l'; %''; DW_lcode_justelectron DW_justelectron
extradatadir = ''; %r2l_302_c_pi_e550_l, 'DWdc3_lcode_x1_pi'; DW_lcode_nofront

dataformat = 'mat';
useAvg = 0;
dump_list = 0:1:10;
% dump_list = 15:1:15;

% saving plot
save_flag = 1;
save_format = {'png'}; % eps, fig

% plasma properties
plasmaden = 2e14; % !!!!!! 2e14, 7e14 

% choose property to plot
property_plot = {'density','wakefields'}; % density, wakefields, both
species_to_plot = {'proton_beam','electron_seed'}; 
include_lineout = 'both'; % 'no','both','field_lineout','density_profile'
include_phasespace = 0; 
field_geometry = 'cartesian'; % cylindrical
on_axis = 'int'; 

% choose fields to plot
wakefields_direction = 'trans'; % trans, long
lineout_point = '0.02'; % cm

% choose species density to plot
species = {'electrons'};    %proton_beam, electrons,  electron_seed

% choose limits (in cm, must denormalize)
trans_range = [0 0.11];
xi_range = [10 0];
plot_density_lims = [0 1e10];
plot_field_lims = [-210 210];

% create movie or not
create_movie = 1;

% choose between normalized and denormalized units
denormalize_flag = 1; % true, false


% choose if make pause or not
make_pause = 0;

% figure number
fig_number = 1;


% directory to save the plots
plots_dir = ['eps2022/fielden/',datadir,''];

switch datadir
    case 'r2l_202'
        extitles{1} = 'e- $1\,\sigma$ behind density cut';
        exlegends{1} = 'e- $1\,\sigma$ behind density cut';
    case 'r2l_2'
        extitles{1} = 'e- 6 ps behind density cut';
        exlegends{1} = 'e- 6 ps behind density cut';
    case 'r2l_2_2s'
        extitles{1} = 'e- 6 ps behind density cut';
        exlegends{1} = 'e- 6 ps behind density cut';
    case 'r2l_202_noe'
        extitles{1} = 'no e-';
        exlegends{1} = 'no e-';
    case 'r2l_202_c'
        extitles{1} = 'e- $0.5\,\sigma$ behind density cut';
        exlegends{1} = 'e- $0.5\,\sigma$ behind density cut';
    case 'r2l_226_c1s_e550'
        extitles{1} = 'e- 1\,cm behind density cut';
        exlegends{1} = 'e- 1\,cm behind density cut';

    case 'r2l_302'
        extitles{1} = 'e- $1.5\,\sigma$ behind density cut';
        exlegends{1} = 'e- $1.5\,\sigma$ behind density cut';
    case 'r2l_401'
        extitles{1} = 'e- $2\,\sigma$ behind density cut';
        exlegends{1} = 'e- $2\,\sigma$ behind density cut';
    case 'r2l_2_longp'
        extitles{1} = ['density cut (3 $\sigma_z$) e- at $\xi = 1$ cm'];
        exlegends{1} = ['density cut (3 $\sigma_z$) e- at $\xi = 1$ cm'];
    case 'r2l_302_c_e550'
        extitles{1} = ['density cut (2 $\sigma_z$) e- at $\xi = 4$ cm'];
        exlegends{1} = ['density cut (2 $\sigma_z$) e- at $\xi = 4$ cm'];
    case 'r2l_302_c_e550_l'
        extitles{1} = ['e- at $\xi = 4$ cm'];
        exlegends{1} = ['e- at $\xi = 4$ cm'];
    otherwise
        extitles{1} = '';
        exlegends{1} = '';

end

switch extradatadir
    case 'r2l_202_pi'
        extitles{2} = 'e- $1\,\sigma + \pi$ behind density cut';
        exlegends{2} = 'e- $1\,\sigma + \pi$ behind density cut';
    case 'r2l_202_c_pi'
        extitles{2} = 'e- $0.5\,\sigma + \pi$ behind density cut';
        exlegends{2} = 'e- $0.5\,\sigma + \pi$ behind density cut';
    case 'r2l_226_c1s_pi_e550'
        extitles{2} = 'e- 1\,cm + $\pi$ behind density cut';
        exlegends{2} = 'e- 1\,cm + $\pi$ behind density cut';
    case 'r2l_2_pi'
        extitles{2} = 'e- 6 ps $+ \pi$ behind density cut';
        exlegends{2} = 'e- 6 ps $+ \pi$  behind density cut';
    case 'r2l_2_pi_2s'
        extitles{2} = 'e- 6 ps $+ \pi$ behind density cut';
        exlegends{2} = 'e- 6 ps $+ \pi$  behind density cut';
    case 'r2l_302_pi'
        extitles{2} = 'e- $1.5\,\sigma + \pi$ behind density cut';
        exlegends{2} = 'e- $1.5\,\sigma + \pi$ behind density cut';
    case 'r2l_401_pi'
        extitles{2} = 'e- $2\,\sigma + \pi$ behind density cut';
        exlegends{2} = 'e- $2\,\sigma + \pi$ behind density cut';
    case 'r2l_2_longp_noe'
        extitles{2} = ['density cut (3 $\sigma_z$) no e-'];
        exlegends{2} = ['density cut (3 $\sigma_z$) no e-'];
    case 'r2l_302_c_pi_e550'
        extitles{2} = ['density cut (2 $\sigma_z$) e- at $\xi = 4 + \lambda_{pe}/2$ cm'];
        exlegends{2} = ['density cut (2 $\sigma_z$) e- at $\xi = 4 + \lambda_{pe}/2$ cm'];
    case 'r2l_302_c_pi_e550_l'
        extitles{2} = ['e- at $\xi = 4 + \lambda_{pe}/2$ cm'];
        exlegends{2} = ['e- at $\xi = 4 + \lambda_{pe}/2$ cm'];
    otherwise
        extitles{2} = '';
        exlegends{2} = '';
end


P = Plotty('datadir',datadir,'extradatadir',extradatadir,'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list,...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'lineout_point',lineout_point,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'plot_density_lims',plot_density_lims,'plot_field_lims',plot_field_lims,...
    'make_pause',make_pause,'fig_number',fig_number,...
    'extitles',extitles,'exlegends',exlegends,...
    'field_geometry',field_geometry,'on_axis',on_axis,...
    'species_to_plot',species_to_plot);%,...
    %'plot_name',datadir);

figure(fig_number);
P.plot_field_density('ylinepos',[-0.02,0.02],'ylineflag',1,...
    'xlinepos',[0.423726],'xlineflag',0,...
    'include_density_profile',1,'include_field_lineout',1); %#ok<NBRAK> 




