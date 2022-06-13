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
load("colororder_defaultblack.mat"); %corder_default

% file location variables

datadirs = {'DWdc_th1','DWdc_th2','DWdc_th3'};

legs ={'th, bunch','th, den. feat.','th, both','sum'};

linestyles = {'-',':','--','-.','-','-','-','-.','--',':'};


dataformat = 'mat';
useAvg = 1;
dump_list = [0];

% saving data
save_flag = 0;
save_format = {'png'};

% plasma properties
plasmaden = 2e14;

% choose property to plot
property_plot = 'wakefields'; % density, wakefields, both

% choose fields to plot
wakefields_direction = 'long'; % trans, long

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 0.02];
xi_ranges = {[19 0]};
lineout_point = '0.005';

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = 1; % true, false

% choose if make pause or not
make_pause = false;

% directory to save the plots
plots_dir = ['DW/','lineoutcompare/',...
    property_plot,'/',wakefields_direction,'/'];
prenames = {'th'};
    
P = Plotty('datadir',datadirs{1},'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list(1),...
    'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'plot_field_lims',[-inf inf],...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'lineout_point',lineout_point,'fig_number',2);

for n = 1:length(dump_list)

    

    P.dump_list = dump_list(n);

    for ixi = 1:length(xi_ranges)
        P.xi_range = xi_ranges{ixi};
        prename = prenames{ixi};

    for d = 1:length(datadirs)
        P.save_flag = save_flag;
        P.datadir = datadirs{d};
        P.lineout_point = lineout_point;

        P.plot_lineout();
        hold on
        P.plot_handle.LineStyle = linestyles{d};
        P.plot_handle.Color = corder_defblack(d,:);
        P.plot_handle.LineWidth = 1;
        if d == 1
        lineout_save = P.lineout;
        elseif d == 2 
            lineout_sum = lineout_save + P.lineout;
        elseif d == 3
            plot(P.dtime + P.simulation_window - P.z,lineout_sum)
        end
        
    end % datadirs
    hold off
    legend(legs{:},'location','best','interpreter','latex','FontSize',P.plot_fontsize-6)
    ylabel('E_z (MV/m)','FontSize', P.plot_fontsize)
    %ylim([-300 300])
    title(['z = ' num2str(0,3),' cm'],'FontSize',P.plot_fontsize);
    %title(['lineout at r = 0.',lineout_point(end-1:end),' mm'])
    %pause(0.5)
    drawnow;

    P.plot_name = [prename,'lineout_comp','',num2str(dump_list(n))];
    P.fig_handle = gcf;
    P.save_flag = 1;
    P.save_plot();
%     clf;
    end % xi ranges
end % dump list

