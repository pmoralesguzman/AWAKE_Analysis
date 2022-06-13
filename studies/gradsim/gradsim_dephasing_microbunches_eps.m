%________________________________________________________________________
% Script to calculate, from the waterfall data, the dephasing of the 
% microbunches close to some selected xi.
% Especially developed for the gradient simulations.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 28/07/2020
%________________________________________________________________________

clear;
close all;
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs = {'gm10','gp10'};
grads_sim = [-20,-15,-10,-5,0,5,10,15,20]/10;

plasmaden = 1.81e14;
property = 'density';
dump_list = 10:100;
useAvg = false;
dataformat = 'mat';

trans_range = [0 0.02]; % if density
% dephasing_xis = [1,7,14,20]; % cm
dephasing_xis = [1,7,14]; % cm
% dephasing_xis = 4;
dephasing_search = 'max'; % 0x, max
force_waterfall = false;

% Plotting parameters
plots_dir = 'gradsim_eps/dephasing/';
save_format = {'png'};
save_flag = 1;

% initialize some plotty
P = Plotty('plasmaden',plasmaden,...
    'save_flag',save_flag,...
    'plots_dir',plots_dir,'save_format',save_format);

OPA = OsirisPhaseAnalysis('property',property,'species','proton_beam',...
    'wakefields_direction','long',...
    'trans_range',trans_range,...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,...
    'dephasing_search',dephasing_search,...
    'force_waterfall',force_waterfall);

for ph = 1:length(dephasing_xis)
    dephasing_xi = dephasing_xis(ph);
    OPA.dephasing_xi = dephasing_xi;
    
    close all
    for d = 1:length(datadirs)
        
        datadir = datadirs{d};
        OPA.datadir = datadir;
        OPA.dephasing_first = [];
        
        OPA.dephasing();
        OPA.dump = 100;
        OPA.count_microbunches();
        OPA.microbunch_number;
        
        % waterfall plot datas
        P.waterfall_xi = OPA.waterfall_xi;
        P.waterfall_z = OPA.waterfall_z;
        P.waterfall_mat = OPA.waterfall_mat;
        P.property = OPA.property;
        
        % PLOT WATERFALL
        P.waterfall_plot();
        x_range = dephasing_xi + [-1, 1];
        xlim(x_range);
        ind_x = fliplr((OPA.waterfall_xi < x_range(2)) & (OPA.waterfall_xi > x_range(1)));
        colormap_low = min(OPA.waterfall_mat(:,ind_x),[],'all');
        colormap_upp = max(OPA.waterfall_mat(:,ind_x),[],'all');
        caxis([colormap_low colormap_upp]);
        ax = imgca(P.fig_handle);
        
        % PLOT ADD the dephasing line (zero-crossing or max) to the waterfall plot
       
        
        phase_x_plot = (-OPA.dephasing_line)*OPA.plasma_wavelength + OPA.simulation_window - OPA.dephasing_first;
        phase_y_plot = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line));
        phase_peakvalue = OPA.dephasing_max_amplitude;
        phase_filename = ['save_files/phases/phase_',datadir,'_',num2str(dephasing_xi),'.mat'];
        save(phase_filename,'phase_x_plot','phase_y_plot','phase_peakvalue')
        
        hold on
        plot(phase_x_plot,phase_y_plot,...
            'LineWidth',2,'Parent',ax,'color','r')
        hold off
        
        ax.XDir = 'reverse';
        ax.YDir = 'normal';
        title(['microbunch number ',num2str(OPA.microbunch_number),' (',datadir,')'])
        
        drawnow;
        P.plot_name = ['gradsim_bunch_',datadir,'_',num2str(dephasing_xi)]; % PLOT NAME
        P.save_plot();
        
        lineout_plot(d,:) = fliplr(OPA.lineout);
        dephasing_line_plot(d,:) = OPA.dephasing_line;
        
        microbunch_locations_accumulated(d) = OPA.microbunch_loc;
        microbunch_peak_accumulated(d) = OPA.microbunch_peak;
        
    end % for datadirs
    
    

    % just the dephasing line, holds the figure for the next datadir
%     fig1 = figure(1);
%     for d = 1:length(datadirs)
%         hold on
%         p1 = plot(linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
%             dephasing_line_plot(d,:),'LineWidth',2);
%         hold off
%     end
%     % labels for the dephasing lines
%     ylabel('phase (\lambda_p)')
%     xlabel('propagation distance (m)')
%     legend('gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20','Location','best')
%     xlim([0 10])
%     ylim([-1.6 1]);
%     title(['\xi = ',num2str(dephasing_xi),' cm behind seed. pos., microbunch ',num2str(OPA.microbunch_number)])
%     
%     P.plot_name = ['gradsim_dephasing_bunch_',num2str(dephasing_xi)]; % PLOT NAME
%     P.fig_handle = fig1;
%     P.save_plot();
%     
%     fig2 = figure(2);
%     po = cell(1,length(datadirs));
%     for d = 1:length(datadirs)
%         hold on
%         po{d} = plot(OPA.waterfall_xi,lineout_plot(d,:),'LineWidth',2);
%         hold off
%     end
%     hold on
%     scatter(microbunch_locations_accumulated,microbunch_peak_accumulated,'r','filled')
%     hold off
%     legend(po{:},...
%         'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20','Location','best')
%     ax2 = fig2.CurrentAxes;
%     ax2.XDir = 'reverse';
%     xlabel('\xi (cm)');
%     ylabel('charge (e)')
%     title(['microbunch number ',num2str(OPA.microbunch_number)])
%     xlim(x_range);
%     P.fig_handle = fig2;
%     P.plot_name = ['microbunches_10m_',datadir,'_',num2str(dephasing_xi)];
%     P.save_plot();
%     
%     
%     fig3 = figure(3);
%     for d = 1:length(datadirs)
%         hold on
%         scatter(grads_sim(d),microbunch_locations_accumulated(d)/OPA.plasma_wavelength,100,'filled')
%         hold off
%     end
%     hold on
%     plot(grads_sim,microbunch_locations_accumulated/OPA.plasma_wavelength,'k')
%     hold off
%     xlabel('gradient (%/m)');
%     ylabel('bunch position in \xi (\lambda_p)')
%     title(['microbunch number ',num2str(OPA.microbunch_number)])
%     P.fig_handle = fig3;
%     P.plot_name = ['microbunchposition_10m_',num2str(dephasing_xi)];
%     P.save_plot();
end


