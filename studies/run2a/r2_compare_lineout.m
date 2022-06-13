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

% datadirs = {'r2l_2','r2l_2_pi','r2l_2_pi2','r2l_2_longp_noe'};
% datadirs = {'r2l_202','r2l_202_pi','r2l_202_pi2','r2l_202_noe'};
% datadirs = {'r2l_302','r2l_302_pi','r2l_302_noe'};
% datadirs = {'r2l_2_2s','r2l_2_pi_2s','r2l_2_noe_2s'};
% datadirs = {'r2l_202_c','r2l_202_c_pi','r2l_202_c_noe'};
datadirs = {'r2l_302_c_e550','r2l_302_c_pi_e550','r2l_302_c_pi_e550'};


% legs ={'e- bunch 6 ps behind den. cut','e- bunch 6 ps + $\pi$ behind den. cut',...
%     'e- bunch 6 ps + $\pi/2$ behind den. cut','no e- bunch'};
% legs ={'e- bunch 1$\sigma$ behind den. cut','e- bunch 1$\sigma$  + $\pi$ behind den. cut',...
%     'e- bunch 1$\sigma$  + $\pi/2$ behind den. cut','no e- bunch'};
% legs ={'e- bunch 1.5$\sigma$ behind den. cut','e- bunch 1.5$\sigma$  + $\pi$ behind den. cut',...
%     'no e- bunch'};
% legs ={'e- bunch 2$\sigma$ behind den. cut','e- bunch 2$\sigma$  + $\pi$ behind den. cut',...
%     'no e- bunch'};
legs ={'e- bunch 0.5$\sigma$ behind den. cut','e- bunch 0.5$\sigma$  + $\pi$ behind den. cut',...
    'no e- bunch'};
linestyles = {'-',':','--','-.','-','-','-','-.','--',':'};


dataformat = 'mat';
useAvg = 0;
% dump_list = [0:1:0];
dump_list = [100:1:100];

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
xi_ranges = {[30 0]};
% xi_ranges = {[19 0]};
lineout_point = '0.003';

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = 1; % true, false

% choose if make pause or not
make_pause = false;

% directory to save the plots
plots_dir = ['r2plot/','summary_lineout/',...
    '302ce550','/'];
prenames = {''};
% prenames = {'z1','z2'};
    
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

        if d > 6
            %             Plineout_save = P.lineout;
            %             Pxi_save = P.z_out;
%             P.dump_list = 2;
            P.useAvg = 1;
            P.dataformat = 'mat';
            if d == 6
                P.useAvg = 1;
                P.dump_list = 0;
            end
        end

        P.plot_lineout();
        hold on
        P.plot_handle.LineStyle = linestyles{d};
        P.plot_handle.Color = corder_defblack(d,:);
        P.plot_handle.LineWidth = 1;
        prop_distance_save = P.dtime;
        
    end % datadirs
    hold off
    legend(legs{:},'location','best','interpreter','latex','FontSize',P.plot_fontsize-4)
    ylabel('$E_z$ (MV/m)','FontSize', P.plot_fontsize,'interpreter','latex')
    %ylim([-300 300])
    title(['z = ' num2str(prop_distance_save,4),' cm'],'FontSize',P.plot_fontsize);
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

