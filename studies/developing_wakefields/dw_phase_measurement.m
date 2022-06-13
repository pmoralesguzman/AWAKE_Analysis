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
load("colororder_defaultblack.mat"); %corder_defblack
corder_defblack = [corder_defblack;corder_defblack;corder_defblack;corder_defblack;corder_defblack];

% file location variables

% datadirs = {'DWdc3_lcode_x666'};
% datadirs = {'DWdc3_lcode_x1','DWdc3_lcode_x1_pi4',... 100
%     'DWdc3_lcode_x1','DWdc3_lcode_x1_pi2',... 100 pi/2
%     'DWdc3_lcode_x1','DWdc3_lcode_x1_pi',... 100 pi
%     'DWdc3_lcode_x2','DWdc3_lcode_x2_pi4',... 50
%     'DWdc3_lcode_x4','DWdc3_lcode_x4_pi4',... 25
%     'DWdc3_lcode_x5','DWdc3_lcode_x5_pi4',... 20
%     'DWdc3_lcode_x666','DWdc3_lcode_x666_pi4',... 15
%     'DWdc3_lcode_x666','DWdc3_lcode_x666_pi2',... 15 pi/2
%     'DWdc3_lcode_x666','DWdc3_lcode_x666_pi',... 15 pi
%     'DWdc3_lcode_x10','DWdc3_lcode_x10_pi4',... 10
%     'DWdc3_lcode_x100','DWdc3_lcode_x100_pi4'}; % 1

% datadirs = {'DWdc3_lcode_x10','DWdc3_lcode_x10_pi'};

datadirs = {'DW_lcode_x100','DW_lcode_x100_pi',...
    'DW_lcode_x100','DW_lcode_x100_pi2',...
    'DW_lcode_x50','DW_lcode_x50_pi',...
    'DW_lcode_x50','DW_lcode_x50_pi2',...
    'DW_lcode_x40','DW_lcode_x40_pi',...
    'DW_lcode_x40','DW_lcode_x40_pi2',...
    'DW_lcode_x30','DW_lcode_x30_pi',...
    'DW_lcode_x30','DW_lcode_x30_pi2',...
    'DW_lcode_x20','DW_lcode_x20_pi',...
    'DW_lcode_x20','DW_lcode_x20_pi2',...
    'DW_lcode_x15','DW_lcode_x15_pi',...
    'DW_lcode_x15','DW_lcode_x15_pi2',...
    'DW_lcode_x10','DW_lcode_x10_pi',...
    'DW_lcode_x10','DW_lcode_x10_pi2',...
    'DW_lcode_x1','DW_lcode_x1_pi',...
    'DW_lcode_x1','DW_lcode_x1_pi2',...
    };

datadirs = {'r2l_2','r2l_2_pi',...
    'r2l_2','r2l_2_pi2',...
    'r2l_202','r2l_202_pi',...
    'r2l_202','r2l_202_pi2',...
    'r2l_302','r2l_302_pi',...
    'r2l_401','r2l_401_pi',...
    'r2l_401','r2l_401_pi2',...
    };

% datadirs = {'r2l_202','r2l_202_pi',...
%     'r2l_202','r2l_202_pi2',...
%     };

% datadirs = {'DW_lcode_x100','DW_lcode_x100_pi',...
%     'DW_lcode_x100','DW_lcode_x100_pi2',...
%     };


%   NO LEGEND FOR NOW

% legs ={'lcode, 15\% density','lcode, 15\% density, $\pi/4$ backwards'};

linestyles = {'-',':','--','-.','-','-','-','-.','--',':','-',':','--','-.','-','-','-','--','-.','-','-','-','-.','--',':',':','--','-.','-','-','-','-.','--',':'};

dataformat = 'mat';
useAvg = 0;
dump_list = [100];

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
search_points = [2,2,2,2,10,10,10,10,13,13,17,17,17,17];
% search_points = [10,10,13,17,17];
% search_points(1) = 10.0; % xi, cm

% xi_ranges = {[search_points(1) + 0.5, search_points(1) - 0.5],...
%     [search_points(2) + 0.5, search_points(2) - 0.5],...
%     [search_points(3) + 0.5, search_points(3) - 0.5]};




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
prenames = {'','z1','z2'};

P = Plotty('datadir',datadirs{1},'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list(1),...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'plot_field_lims',[-inf inf],...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'lineout_point',lineout_point,'fig_number',2);

ph_i = 1;
phase_from_peak_list = zeros(20,1);

for n = 1:length(dump_list)

    P.dump_list = dump_list(n);

    for d = 1:length(datadirs)
        for ixi = 1:length(1)
            xi_ranges = {[search_points(d) + 0.5, search_points(d) - 0.5]};
            P.xi_range = xi_ranges{ixi};
            search_point = search_points(d);
            %         prename = prenames{ixi};


            P.datadir = datadirs{d};
            P.lineout_point = lineout_point;
%             P.getdata(); P.denorm_distance();

            if d > 10
                %             Plineout_save = P.lineout;
                %             Pxi_save = P.z_out;
                %             P.dump_list = 2;
                P.useAvg = 0;
                P.dataformat = 'mat';
                if d == 9
                    P.useAvg = 1;
                    P.dump_list = 0;
                end
            else
                prop_distance_save = P.propagation_distance;

            end

            P.plot_lineout(); %P.lineout, P.z_out
            hold on
            P.plot_handle.LineStyle = linestyles{d};
            P.plot_handle.Color = corder_defblack(d,:);
            P.plot_handle.LineWidth = 1;

            study_window = min(abs(diff(P.xi_range)),P.simulation_window);
            smooth_profile = smooth(P.lineout,P.plasma_wavelength/study_window/8,'loess');

            peak_dis = 0.85*P.plasma_wavelength;
            [pks,locs] = findpeaks(smooth_profile,-P.z_out,...
                'MinPeakDistance',peak_dis);
            large_pks_ind = pks > 0.1*max(pks);
            locs = -locs(large_pks_ind);
            pks = pks(large_pks_ind);

            [~,loc_closest_to_search_point] = min(abs(search_point - locs));

            locs_to_search = loc_closest_to_search_point + [-1,0,1];

            phase_from_peak = mean(rem(locs(locs_to_search),P.plasma_wavelength))/P.plasma_wavelength*360;

            phase_from_peak_list(ph_i) =  phase_from_peak;
            ph_i = ph_i + 1;

            scatter(locs(locs_to_search),pks(locs_to_search))



        end % datadirs
        hold off
        %     legend(legs{:},'location','best','interpreter','latex','FontSize',P.plot_fontsize-6)
        ylabel('E_z (MV/m)','FontSize', P.plot_fontsize)
        %ylim([-300 300])
        title(['z = ' num2str(prop_distance_save,3),' cm'],'FontSize',P.plot_fontsize);
        %title(['lineout at r = 0.',lineout_point(end-1:end),' mm'])
        %pause(0.5)
        drawnow;

        %     P.plot_name = [prename,'lineout_comp','',num2str(dump_list(n))];
        P.fig_handle = gcf;
        P.save_flag = 0;
        P.save_plot();
        %     clf;
    end % xi ranges
end % dump list

diff_list = abs(phase_from_peak_list(2:2:end) - phase_from_peak_list(1:2:end-1));

densities = [100,100,50,50,40,40,30,30,20,20,15,15,10,10,1,1]';
densities = [0,0,1,1,1.5,2,2]';

save('density_vs_phasediff.mat','diff_list','densities')

