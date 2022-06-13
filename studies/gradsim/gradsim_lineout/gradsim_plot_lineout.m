%________________________________________________________________________
% All gradients
% Plot lineout standard file
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 22/02/2021
%________________________________________________________________________


% data directory
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'}; % data directory
plots_dir = ['gradsim/lineout/']; % save directory for the plot

% parameters
plasma_density = 1.81e14;

% properties
property_plot = 'wakefields';
species = 'proton_beam';
field = 'e';
direction = 'z';
wakefields_direction = 'long';

% simulation parameters
dump_list = 0:1:100;
dataformat = 'mat';
useAvg = false;

% plot parameters
lineout_point = 5;
xi_range = [15,0];


P = Plotty('plasmaden',plasma_density,'property_plot',property_plot,...
    'species',species,'field',field,'direction',direction,...
    'wakefields_direction',wakefields_direction,...
    'lineout_point',lineout_point,'xi_range',xi_range,...
    'dump_list',dump_list,'dataformat',dataformat,'useAvg',useAvg,...
    'create_movie',true);

for d = 1:length(datadirs)
    
    P.datadir = datadirs{d};
    P.plots_dir = ['gradsim/lineout/',datadirs{d}]; % save directory for the plot
    P.plot_lineout();
    
end




