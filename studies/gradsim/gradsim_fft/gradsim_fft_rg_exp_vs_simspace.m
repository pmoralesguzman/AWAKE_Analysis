%________________________________________________________________________
% Freq vs g for two r's. 
% FFT of the proton or plasma electrons density distribution, or of a
% lineout of the wakefields. Special version to produce the plot for
% the paper, which shows the fft for a narrow and wide window.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% EDA. I. Morales Guzman
% Last update: 05/10/2020
%________________________________________________________________________

close all;
clear;
% set(0,'defaultAxesFontSize',20)
% load experimental data from Fabian
freqs_exp      = load('gradsim_freqs.txt');
grads_exp      = freqs_exp(:,1)/10;
freqs_exp(:,1) = [];
% 1: full green, 2: CTR black, 3: narrow red, 4: wide blue

% data directory
datadirs_exp    = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs    = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
dataformat  = 'mat';

% dataformats  = {'h5','h5','h5','mat','mat','mat','mat','mat','mat'};
dataformats  = {'mat','mat','mat','mat','mat','mat','mat','mat','mat'};

useAvg      = false;
dump        = 100;


% save directory
plots_dir           = ['gradsim/fft/rg/',''];
plot_name_suffix    = [''];
save_format         = {'png','eps','fig'};

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';
grads_sim = [-20,-15,-10,-5,0,5,10,15,20]/10;

% analysis
xi_range        = [14 0.38]; %384
% trans_range is set by trans_lims in this analysis


scan_type       = 'cumulative'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout
freq_range      = [107 128]; % GHz ? maybe norm. to plasma freq
max_dotsize     = 720;
showChargeinDotSize = true; % show dot size to reflect charge


% switches
plot_powerspectra   = false;
save_plot_flag      = false;
normalized_frequency_switch = true;

% calculated parameters
trans_lims      = [0.0604/2 0.0604*3];
nslices         = length(trans_lims);
% measured radius = 0.66 mm,
% first exp limit = 0.868 mm
% last exp limit 2.604 mm
% theo radius at plasma exit = 0.0455, 0.0455*1.3152 = 0.0598 cm (not used
% anymore)

% colors
lincolngreen = [20,82,20]/256;

% initialize variables
peak_freqs = zeros(nslices,length(dump));
dotsize = max_dotsize/10*ones(nslices,length(dump));
charge_in_slice = ones(nslices,length(dump));

% classes
AFFT = AwakeFFT(...
    'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'xi_range',xi_range,'trans_lims',trans_lims,...
    'scan_type',scan_type,'on_axis','sum');

AFFT2 = AwakeFFT('load_data_flag',0,...
    'datadir',datadirs_exp{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'xi_range',xi_range,'trans_lims',trans_lims,...
    'scan_type','cumulative','on_axis','sum');

P = Plotty('plots_dir',plots_dir,'plot_name',['fvsg',''],'plasmaden',plasma_density,'save_flag',save_plot_flag,...
    'datadir',datadirs{1});
EDA = ExperimentalDataAnalyser('plasmaden',plasma_density,'datadir',datadirs_exp{1});
EDA.loadSCdata();

AFFT2.dr = EDA.SCI_r(2) - EDA.SCI_r(1);


for d = 1:length(datadirs)
    
    if d == 2
        EDA.datadir = datadirs_exp{d-1}; EDA.loadSCdata();
    else
        EDA.datadir = datadirs_exp{d}; EDA.loadSCdata();
        
    end
    
    AFFT.datadir = datadirs{d};
    AFFT.dataformat = dataformats{d};
    AFFT.fft_dataload();
    
    AFFT2.r = EDA.SCI_r;
    % trim
    z = EDA.SCI_z*1e-12*EDA.c_cm;
    
    z_ind = z > xi_range(2) & ... %large
        z <= xi_range(1); % small
    AFFT2.z = z(z_ind);
    
    AFFT2.proton_beam = EDA.SCI(:,z_ind);
    
    AFFT2.fft_dataload();
    prop_distance_m = AFFT.propagation_distance/100; % propagation distance in m
    
    % calculates fft and gives AFFT.fft_frequencies
    % and AFFT.fft_powerspectrum (_den or _fld)
    AFFT.get_fft(); AFFT2.get_fft();
    
    switch AFFT.property
        case 'density'
            charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
            powerspectrum_in_slice = AFFT.fft_powerspectrum_den;
            powerspectrum_in_slice2 = AFFT2.fft_powerspectrum_den;
            
        case 'fields'
            powerspectrum_in_slice = AFFT.fft_powerspectrum_fld;
    end % switch property
    
    for r = 1:nslices
        
        AFFT.fft_peaks(powerspectrum_in_slice(r,:));
        AFFT2.fft_peaks(powerspectrum_in_slice2(r,:));

        % max_loc is the location in freq space of the ampltiude peak in
        % the power spectrum
        if isempty(AFFT.max_loc)
            peak_freqs(r,d) = 0;
        else
            peak_freqs(r,d) = AFFT.max_loc;
            peak_freqs2(r,d) = AFFT2.max_loc;
        end
        
        %%% Plots of the power spectrum.
        if plot_powerspectra
            P2 = P;
            trans_lims_plot_title = [0,trans_lims];
            fft_freqs_plot = AFFT.fft_frequencies(AFFT.ind_pksearch)/1e9;
            amplitude_plot = powerspectrum_in_slice(r,AFFT.ind_pksearch);
            P2.plots_dir = [plots_dir,'/power_spectra'];
            P2.plot_name = ['ps',AFFT.datadir,'n',num2str(AFFT.dump),...
                '_s',num2str(r)];
            
            fig_powerspectra = figure('visible','off');
            plot(fft_freqs_plot,amplitude_plot,'LineWidth',2);
            title(['power spec. (proEDA. dist. = ',...
                num2str(AFFT.propagation_distance/100,2),'m; r = [',...
                num2str(trans_lims_plot_title(r)*10,2),', ',num2str(trans_lims_plot_title(r+1)*10,2),']mm)']);
            ylabel('amplitude');
            xlabel('frequency (GHz)');
            xlim([fft_freqs_plot(1),fft_freqs_plot(end)]);
            drawnow;
            P2.fig_handle = fig_powerspectra;
            P2.save_plot();
            close all;
            AFFT.progress_dump('save slice fft',r,nslices);
            
            indz = AFFT.z - AFFT.dtime > 7;
            z_plot = AFFT.z(indz);
            
            density_plot = AFFT.fft_densitymatrix(r,indz);
            P2.plots_dir = [plots_dir,'/density_distributions'];
            P2.plot_name = ['denlong',AFFT.datadir,'n',num2str(AFFT.dump),...
                '_s',num2str(r)];
            
            figy = figure('visible','off');
            plot(z_plot,density_plot,'LineWidth',2);
            title(['long. den. distributon (proEDA. dist. = ',...
                num2str(AFFT.propagation_distance/100,2),'m; r = [',...
                num2str(trans_lims_plot_title(r)*10,2),', ',num2str(trans_lims_plot_title(r+1)*10,2),']mm)']);
            ylabel('density (1/cm)');
            xlabel('z-ct (m)');
            xlim([z_plot(1),z_plot(end)]);
            figy.Units = 'normalized';
            figy.Position = [0.05 0.2 0.9 0.5];
            drawnow;
            P2.fig_handle = figy;
            P2.save_plot();
            close all;
            AFFT.progress_dump('save slice density',r,nslices);
        end
        
    end % for nslices
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize==0) = nan;
    end
end % for datadirs

peak_freqs_plot = peak_freqs;
freqs_exp = peak_freqs2(:,[1,3:9]);

%% plotty
fig_fvsrg = figure();
ax_fvsrg = axes('Parent',fig_fvsrg);
hold on
if normalized_frequency_switch
    marker_size = 10;
    % simulation
    plot_sim(1) = plot(grads_sim,... % narrow trans. range
        peak_freqs_plot(1,:)/AFFT.plasmafreq_GHz,'-o','markerSize',marker_size);
    set(plot_sim(1), 'markerfacecolor', get(plot_sim(1), 'color'),'LineWidth',2);
    plot_sim(2) = plot(grads_sim,... % wide trans. range
        peak_freqs_plot(end,:)/AFFT.plasmafreq_GHz,'-s','markerSize',marker_size);
    set(plot_sim(2), 'markerfacecolor', get(plot_sim(2), 'color'),'LineWidth',2);
    
    % experiment
    plot_exp(1) = plot(grads_exp,... % narrow trans. range
        freqs_exp(1,:)/AFFT.plasmafreq_GHz,'-^','markerSize',marker_size,'LineWidth',2,...
        'color',lincolngreen);
    plot_exp(2) = plot(grads_exp,... % wide trans. range
        freqs_exp(end,:)/AFFT.plasmafreq_GHz,'-d','markerSize',marker_size,'LineWidth',2);
    
    ylim(freq_range/AFFT.plasmafreq_GHz)
    
    ylabel({'normalized modulation frequency (a. u.)'})
    
    legend([plot_sim(1) plot_sim(2) plot_exp(1) plot_exp(2)],...
        'Simulation narrow window','Simulation wide window',...
        'Experiment narrow window','Experiment wide window',...
        'location','southeast','AutoUpdate','off');
    
    ax_handle = gca;
    plot([-2.1 2.1],[1,1],'b','LineWidth',2)
    plot(grads_sim,[0.894394264356599,0.921920505722216,0.948648374071357,0.974643553503204,0.999963186893636,1.02465735438711,1.04877023825847,1.07234105181787,1.09540478827280],'b','LineWidth',2)
    
else
    scatter(trans_lims_m,...
        peak_freqs_plot,dotsize,'filled');
    %         ylim([min(peak_freqs_plot)*0.99,max(peak_freqs_plot)*1.01])
    ylim(freq_range);
    ylabel('frequency (GHz)')
end % if normalized frequency switch

%     ax_handle = gca;
%     plot(3*0.2*sqrt(1+100/4.9^2)*ones(1,2),ax_handle.YLim,'--r','LineWidth',2) % dashed line at 3 sigma_r (bunch travelling in vacuum)
%     plot(ax_handle.XLim,freq_widerange*ones(1,2),'--b','LineWidth',2);

xlim([-2.1 2.1])
hold off
grid on
xlabel('gradient (%/m)')

P.save_plot(fig_fvsrg);




