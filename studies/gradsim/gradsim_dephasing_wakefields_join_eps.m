%________________________________________________________________________
% Script to calculate, from the waterfall data, the dephasing of the fields
% close to some selected xi.
% Especially developed for the APS plots.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/07/2020
%________________________________________________________________________

clear;
close all;
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs = {'gm10','g0','gp10'};
g = {'g = -1 %/m','g = 0 %/m','g = +1 %/m'};
letters = {'a) -1 %/m','b)','c)','d) 0 %/m','e)','f)','g) +1 %/m','h)','i)','j)','k)','l)'};
letterbox_x = 0.05;
letterbox_y = 0.15;
% initialize counters
i_letter = 1;
plasmaden = 1.81e14;
property = 'fields';
dump_list = 0:100;
useAvg = 1;
dataformat = 'h5';

trans_range = [0 0.02]; % if density
% dephasing_xis = [1,1.5,2,3,4,5,6,7]; % cm
dephasing_xis = [14,7,1];
dephasing_search = '0x'; % 0x, max
force_waterfall = false;

% initialize some plotty
P1 = Plotty('plots_dir','gradsim_eps/dephasing','save_format',{'png','pdf','eps'});
P2 = Plotty('plots_dir','gradsim_eps/dephasing');

fig_w = figure(1);
P1.fig_handle = fig_w;
tt = tiledlayout(3,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
i = 1;
for d = 1:length(datadirs)
    datadir = datadirs{d};
%     close all;
    for ph = 1:length(dephasing_xis)
        dephasing_xi = dephasing_xis(ph);
        
        
        
        OPA = OsirisPhaseAnalysis('datadir',datadir,...
            'property',property,'species','proton_beam',...
            'direction','z',...
            'wakefields_direction','long',...
            'trans_range',trans_range,...
            'plasmaden',plasmaden,...
            'dump_list',dump_list,'useAvg',useAvg,...
            'dataformat',dataformat,...
            'dephasing_xi',dephasing_xi,...
            'dephasing_search',dephasing_search,...
            'force_waterfall',force_waterfall,...
            'lineout_point',8);
        OPA.dephasing();
        
        % waterfall plot datas
        P1.waterfall_xi = OPA.waterfall_xi;
        P1.waterfall_z = OPA.waterfall_z;
        P1.waterfall_mat = OPA.waterfall_mat;
        P1.property = OPA.property;
        
        ax(i) = nexttile;
       
%         P1.fig_number = 10+d;
        P1.waterfall_plot();
        x_range = dephasing_xi+[-1 1];
        %         x_range = [0 2.5];
        xlim(x_range);
        ind_x = fliplr((OPA.waterfall_xi < x_range(2)) & (OPA.waterfall_xi > x_range(1)));
        colormap_low = min(OPA.waterfall_mat(:,ind_x),[],'all');
        colormap_upp = max(OPA.waterfall_mat(:,ind_x),[],'all');
        caxis(ax(i),[colormap_low colormap_upp]);
        colormap(ax(i),bluewhitered);
        
%         ax = imgca(P1.fig_handle);
%         ax = imgca(gcf);
        
        % load line from the microbunches for double plotting
        phase_filename = ['save_files/phases/phase_',datadir,'_',num2str(dephasing_xi),'.mat'];
        bunches = load(phase_filename);
        bunches_x = bunches.phase_x_plot;
        bunches_y = bunches.phase_y_plot;
        peakvalue = bunches.phase_peakvalue/max(bunches.phase_peakvalue);
        % add the dephasing line (zero-crossing or max) to the waterfall plot
        hold on
        plot((-OPA.dephasing_line)*OPA.plasma_wavelength+OPA.simulation_window-OPA.dephasing_first,...
            linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
            'LineWidth',2,'Parent',ax(i),'color','k')
        scatter(bunches_x,bunches_y,'LineWidth',2,'Parent',ax(i),...
            'markeredgecolor','none','markerfacecolor',[0, 100/256, 0], 'alphadata', peakvalue, 'markerfacealpha', 'flat');
        ax(i).XDir = 'reverse';
        ax(i).YDir = 'normal';
        i = i + 1;
        t = text(letterbox_x,letterbox_y,letters{i_letter},'units','normalized','FontUnits','centimeters','FontSize',0.5);
        i_letter = i_letter + 1;
        if ph == 1
%             ylabel(g{d})
        end
        if ph == 3 && d == 2
        cbar = colorbar;
        cbar.Label.String = 'E_z (MV/m)';
        end
        
        %         switch datadir
        %             case 'g0'
        %               title({'waterfall plot of longitudinal wakefields on axis','constant density'})
        %             case 'gp20'
        %         %         title(['\xi = ',num2str(dephasing_xi),' cm behind seed. pos. (',datadir,')'])
        %         title(['positive gradient = +2 %/m'])
        %             case 'gm20'
        %                title(['negative gradient = -2 %/m'])
        %
        %         end
        hold off
        P1.fig_handle.Units = 'centimeters';
        P1.fig_handle.OuterPosition = [1 1 8.6*1.5 8.6]*2;

        drawnow;
        
        
        % just the dephasing line, holds the figure for the next datadir
        %         fig1 = figure(1);
        %         colororder(cc);
        %         hold on
        %         p1 = plot(linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
        %             OPA.dephasing_line,'LineWidth',2);
        %         hold off
        
    end
    % labels for the dephasing lines
    %     ylabel('phase (\lambda_p)')
    %     xlabel('propagation distance (m)')
    %     legend('gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20','Location','best')
    %     xlim([0 10])
    %     % ylim([-6 1]);
    %     title(['\xi = ',num2str(dephasing_xi),' cm behind seed. pos.'])
    %
    %     P2.plot_name = ['gradsim_dephasing_fields_',num2str(dephasing_xi)]; % PLOT NAME
    %     P2.fig_handle = fig1;
    %     P2.save_plot();
    
end

ylabel(tt,'z (m)');
    xlabel(tt,'\xi (cm)');

P1.plot_name = ['joinwaterfall']; % PLOT NAME
P1.save_plot();



