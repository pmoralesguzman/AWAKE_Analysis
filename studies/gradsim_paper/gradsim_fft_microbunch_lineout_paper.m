%________________________________________________________________________
% Plot the lineout of the microbunches near the axis
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 02/07/2021
%________________________________________________________________________

clear;
close all;
co =[0 0 0;0.9336 0.2305 0.1719];

% file location variables
datadir = 'gm20';
legs = {'z = 3.8 m','z = 10 m'};

linestyles = {'-','-.'};
dataformat = 'h5';
useAvg = 1;
dump_list = [38,100];

% saving data
save_flag = 1;
save_format = {'fig'};
plot_name = 'microbunchgm20';

% plasma properties
plasmaden = 1.81e14;

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 0.018];
xi_range = [2 0];

% choose property to plot
property_plot = 'density'; % density, wakefields, both

% choose between normalized and denormalized units
denormalize_flag = true; % true, false


% directory to save the plots
plots_dir = 'gradsim_paper/fft/lineout';

P = Plotty('datadir',datadir,'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list(1),...
    'save_flag',0,'save_format',save_format,'plots_dir',plots_dir,...
    'plot_name',plot_name,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag);

hold on
for n = 1:length(dump_list)
    
    P.dump_list = dump_list(n);
    P.plot_lineout();
    P.plot_handle.LineStyle = linestyles{n};
    P.plot_handle.Color = co(n,:);
    P.plot_handle.LineWidth = 2;
    
end
hold off


ylabel('density (arb. unit)','FontSize',14)
xlabel('\xi (cm)','FontSize',14)
P.fig_handle.Children.FontSize = 14;

P.save_flag = 1; % needs to be here
P.plot_name = plot_name; % needs to be here
P.fig_handle = gcf;
P.save_plot();

ylim([0 4.5e11])
legend(legs{:},'location','best','FontSize',12)

