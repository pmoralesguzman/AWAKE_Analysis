%________________________________________________________________________
% Plot the lineout of the microbunches near the axis together with the
% envelope
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
co =[0 0 0;0.9336 0.2305 0.1719;0.9336 0.2305 0.1719];

datadir = 'gp20';
dataformat = 'h5';
useAvg = 1;
dump_list = [38,100];

legs = {'z = 3.8 m','z = 10 m','z = 10 m (envelope)'};
linestyles = {'-','-.'};

% saving data
save_flag = 1;
save_format = {'fig'};

% plasma properties
plasmaden = 1.81e14;

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 0.018];
xi_range = [14 0];

% choose property to plot
property_plot = 'density'; % density, wakefields, both

% choose between normalized and denormalized units
denormalize_flag = true; % true, false


% directory to save the plots
plots_dir = 'gradsim_paper/fft/lineout/';
plot_name = ['microbunch',datadir];

P = Plotty('datadir',datadir,'dataformat',dataformat,'useAvg',useAvg,...
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
    
    if n == 2 % <=
        [yupper,ylower] = envelope(P.lineout,250,'peak'); % 290 for both
        plot(P.z_out,yupper,'LineWidth',2,'Color',co(n,:))
    end
end
hold off

ylabel('density (arb. unit)','FontSize',14)
xlabel('\xi (cm)','FontSize',14)
P.fig_handle.Children.FontSize = 14;
ylim([0 1.8*4.5e11])
legend(legs{:},'location','best','FontSize',12)
title(['propagation dist. = ' num2str(P.propagation_distance,3),' cm']);

P.save_flag = 1; % needs to be here
P.plot_name = plot_name; % needs to be here
P.fig_handle = gcf;
P.save_plot();

