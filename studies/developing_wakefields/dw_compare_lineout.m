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

% datadirs = {'DWdc3_lcode_x1','DWdc3_lcode_x1_pi4','DWdc3_lcode_x1_pi2','DWdc3_lcode_x1_pi'};
% datadirs = {'DWdc_lcode_x666','DWdc_lcode_x666_pi','DWdc_lcode_nofront_nodf'};
datadirs = {'DW_lcode_x5','DW_lcode_x5_r2','DW_lcode_x5_z2'};
datadirs = {'DW_lcode_x5','DW_lcode_x5_t2','DW_lcode_x5_th','DW_lcode_x5_t4','DW_lcode_x5_tq'};
% datadirs = {'DW_lcode_x5','DW_lcode_x5_r2','DW_lcode_x5_rh','DW_lcode_x5_r4'};
% datadirs = {'DW_lcode_x5','DW_lcode_x5_z2','DW_lcode_x5_zh','DW_lcode_x5_z4'};
datadirs = {'DW_lcode_x5','DW_lcode_x5_rz2','DW_lcode_x5_rz2_th'};
% datadirs = {'DW_lcode_x5','DW_lcode_x5_m2','DW_lcode_x5_mh','DW_lcode_x5_m4'};
% datadirs = {'DW_lcode_x5','DW_lcode_x5_th','DW_lcode_x5_tq'};
datadirs = {'DW_lcode_x5_rz2_th','DW_lcode_x5_rz2_th_m2','DW_lcode_x5_rz2_th_mh'};
datadirs = {'DW_lcode_x40','DW_lcode_x40_pi'};
datadirs = {'DW_lcode_x1','DW_lcode_x1_pi2','DW_lcode_x1_pi'};
datadirs = {'DW_lcode_x50','DW_lcode_x50_pi2','DW_lcode_x50_pi'};

datadirs = {'r2l_302_c','r2l_302_c_pi'};


% legs ={'lcode, 15\% density','lcode, 15\% density, $\pi/4$ backwards','lcode, 15\% density, $\pi$ backwards'};
% legs ={'lcode, 100\% density','lcode, 100\% density, $\pi/4$ backwards',...
%     'lcode, 100\% density, $\pi/2$ ba ckwards','lcode, 100\% density, $\pi$ backwards'};
legs ={'baseline','half r','half z'};
legs ={'baseline','timestep x2','timestep x0.5','timestep x4','timestep x0.25'};
legs ={'baseline','t th','rz2 th'};
% legs ={'baseline','r x2','r x0.5','r x4'};
% legs ={'baseline','z x2','z x0.5','z x4'};
% legs = {'x40','x40 pi'};
% legs = {'front 50\%','front 50\% dephased $\pi/2$','front 50\% dephased $\pi$'};
% legs = {'front 1\%','front 1\% dephased $\pi/2$','front 1\% dephased $\pi$'};
legs ={'baseline','pi','noe'};

% legs ={'baseline','completely out-of-phase','completely out-of-phase x2 density'};
% legs ={'baseline','baseline theory','no front of the bunch'};
linestyles = {'-',':','--','-.','-','-','-','-.','--',':'};


dataformat = 'mat';
useAvg = 0;
% dump_list = [0:2:100];
dump_list = [0:1:0];

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
xi_ranges = {[20 0]};
% xi_ranges = {[19 0]};
lineout_point = '0.003';

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = 1; % true, false

% choose if make pause or not
make_pause = false;

% directory to save the plots
plots_dir = ['DW_PEB/','lineoutcompare/',...
    property_plot,'/',wakefields_direction,'/'];
prenames = {'','z1','z2'};
    
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
    legend(legs{:},'location','best','interpreter','latex','FontSize',P.plot_fontsize-6)
    ylabel('E_z (MV/m)','FontSize', P.plot_fontsize)
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

