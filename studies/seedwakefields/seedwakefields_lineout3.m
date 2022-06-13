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
datadirs = {'g0z2','g0z2-n'};


linestyles = {'-','-.'};
dataformat = 'h5';
useAvg = 1;
dump_list = [3];

% saving data
save_flag = 1;
save_format = {'png','fig'};
plot_name = 'seedwakefields';

% plasma properties
plasmaden = 1.81e14;

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 0.018];
xi_range = [19 0];
% xi_range = [19 15];

% choose property to plot
property_plot = 'wakefields'; % density, wakefields, both
wakefields_direction = 'long';

% choose between normalized and denormalized units
denormalize_flag = true; % true, false


% directory to save the plots
plots_dir = 'seedwakefields';

P = Plotty('dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list(1),...
    'save_flag',0,'save_format',save_format,'plots_dir',plots_dir,...
    'plot_name',plot_name,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'species',species,'lineout_point',3,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag);

xx = 1;
hold on
for d = 1:length(datadirs)
    P.datadir = datadirs{d};
    
    for n = 1:length(dump_list)
        
        P.dump_list = dump_list(n);
        P.plot_lineout();
        P.plot_handle.LineStyle = linestyles{xx};
        P.plot_handle.Color = co(xx,:);
        xx = xx + 1;
        P.plot_handle.LineWidth = 2;
        z(n) = P.dtime; % position of the front
        
    end

end
hold off


ylabel('E_z on axis','FontSize',14)
xlabel('\xi (cm)','FontSize',14)
P.fig_handle.Children.FontSize = 14;
P.fig_handle.Children(1).Title = [];

if length(dump_list) > 1
    legs = {['RIF at z = ',num2str(z(1),2),' cm'],['RIF at z = ',num2str(z(2),2),' cm']};
elseif length(datadirs) > 1
    legs = {['RIF at \xi = 128 ps (ahead of bunch center)'],...
        ['RIF at \xi = -128 ps (behind bunch center)']};
end
legend(legs{:},'location','best','FontSize',12)


P.save_flag = 1; % needs to be here
P.plot_name = plot_name; % needs to be here
P.fig_handle = gcf;

P.save_plot();

% ylim([0 4.5e11])

