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
clear;

% file location variables
datadirs = {'TFHbl','TFHblsim'};
% legs = {'plasma rad. = 1 mm','plasma rad. = 2 mm','theory'};
legs = {'plasma rad. = 1 mm','theory'};
linestyles = {'-','-.','--'};
co =[0 0 0;0.9336 0.2305 0.1719;0.0 0.5 0.15];


dataformat = 'h5';
useAvg = 1;
dump_list = 1;

% saving data
save_flag = 0;
save_format = {'png','eps','fig'};

% plasma properties
plasmaden = 2e14;

% choose fields to plot
wakefields_direction = 'long'; % trans, long

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 0.2];
xi_range = [21 0];
lineout_point = '0.01';

% choose property to plot
property_plot = 'wakefields'; % density, wakefields, both

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = 1; % true, false

% choose if make pause or not
make_pause = false;

% directory to save the plots
plots_dir = ['TFH/','lineoutcompare/',...
    property_plot,'/',wakefields_direction,'/'];

P = Plotty('datadir',datadirs{1},'dataformat','h5',...
    'useAvg',useAvg,'dump_list',dump_list(1),...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'lineout_point',lineout_point);

for n = 1:length(dump_list)
    
    P.dump_list = dump_list(n);
    
    for d = 1:length(datadirs)
        P.datadir = datadirs{d};
        P.lineout_point = lineout_point;
        if d == 2
            P.useAvg = 1;
            P.dataformat = 'mat';
            %P.useAvg = 0;
            P.dump_list = 0;
        end
        hold on
        P.plot_lineout();
        P.plot_handle.LineStyle = linestyles{d};
        P.plot_handle.Color = co(d,:);
        P.plot_handle.LineWidth = 2;
        hold off
    end
    
    legend(legs{:},'location','southwest','interpreter','latex')
    %      legend(legs{:},'location','best')
    
    P.plot_name = ['lineout_comp','',num2str(n),lineout_point(end-1:end)];
    P.fig_handle = gcf;
    ylabel('E_z (MV/m)')
    %     title(['propagation dist. = ' num2str(P.propagation_distance,3),' cm']);
    title(['lineout at r = 0.',lineout_point(end-1:end),' mm'])
    P.save_flag = 1;
    P.save_plot();
end

