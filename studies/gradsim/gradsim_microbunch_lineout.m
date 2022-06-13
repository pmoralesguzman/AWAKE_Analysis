%________________________________________________________________________
% Plot the 2D wakefields together with the proton bunch, or each
% individually.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 11/02/2020
%________________________________________________________________________

close all;
co =[0 0 0;0.9336 0.2305 0.1719];

% file location variables
% datadirs = {'g0','gradsim_shortbunch'};
% legs = {'g = 0 %/m',...
%     'short bunch'};
datadirs = {'gm20','gm20'};
legs = {'z = 3.8 m','z = 10 m'};
% datadirs = {'gradsim_shortbunch'};
% legs = {'short bunch'};
linestyles = {'-','-.'};
dataformat = 'mat';
useAvg = 0;
dump_list = [38,100];

% saving data
save_flag = 0;
save_format = {'png','eps','fig'};

% plasma properties
plasmaden = 1.81e14;

% choose fields to plot
wakefields_direction = 'long'; % trans, long

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 0.018];
xi_range = [2 0];
lineout_point = [25,25,50];

% choose property to plot
property_plot = 'density'; % density, wakefields, both

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose if make pause or not
make_pause = false;

% directory to save the plots
plots_dir = ['gradsim_paper/'];

P = Plotty('datadir',datadirs{1},'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list(1),...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'lineout_point',lineout_point(1));

for n = 1:length(dump_list)
    
    P.dump_list = dump_list(n);
    
    hold on
    P.plot_lineout();
    P.plot_handle.LineStyle = linestyles{n};
    P.plot_handle.Color = co(n,:);
    P.plot_handle.LineWidth = 2;
    P.plot_handle.LineWidth = 2;
    hold off
    
    title('')
    legend(legs{:},'location','best','FontSize',12)
    
    P.plot_name = ['microbunchgm20','n',num2str(n)];
    P.fig_handle = gcf;
    ylabel('density (arb. unit)','FontSize',14)
    xlabel('\xi (cm)','FontSize',14)
    P.fig_handle.Children(2).FontSize = 14;
    ylim([0 4.5e11])
%     title(['propagation dist. = ' num2str(P.propagation_distance,3),' cm']);
    P.save_flag = 1;
    P.save_plot();
end

