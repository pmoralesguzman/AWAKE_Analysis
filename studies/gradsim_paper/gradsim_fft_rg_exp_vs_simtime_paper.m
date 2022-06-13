%________________________________________________________________________
% Freq vs g for two r's.
% FFT of the proton or plasma electrons density distribution, or of a
% lineout of the wakefields. Special version to produce the plot for
% the paper, which shows the fft for a narrow and wide window.
% Version which does the plot for publication.
%
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 07/04/2021
%________________________________________________________________________

close all;
clear;

freqs_exp      = load('gradsim_freqs.txt');
grads_exp      = freqs_exp(:,1)/10;

% data directory
datadirs_exp = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs     = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
dataformat   = 'mat';

% dataformats  = {'h5','h5','h5','mat','mat','mat','mat','mat','mat'};
dataformats  = {'mat','mat','mat','mat','mat','mat','mat','mat','mat'};

useAvg      = false;
dump        = 133;
load_dump   = 100;
measurement_distance = 1350;
slit_flag = 1;

% save directory
plots_dir           = ['gradsim_paper/fft/rg/',''];
plot_name_suffix    = [''];
save_format         = {'png','eps','fig'};

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';
grads_sim = [-20,-15,-10,-5,0,5,10,15,20]/10;

% analysis
xi_range          = [14 0.345]; %[14 0.345], new range
% xi_range        = [21 0.742]; % Fabian range, 742

% trans_range is set by trans_lims in this analysis

scan_type       = 'slice_abs'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout
freq_range      = [107 128]; % GHz ? maybe norm. to plasma freq

% switches
plot_powerspectra   = 0;
save_plot_flag      = 1;

% theo radius = 0.02*sqrt(1+13.5^2/4.7382^2) = 0.0604
trans_lims      = [0.0364 0.0364 0.0604*3.06]; %[0.0359 0.0359 0.0604*3.06]
% trans_lims      = [0.018 0.0604*4];
nslices         = length(trans_lims);


% measured radius = 0.66 mm,
% first exp limit = 0.868 mm
% last exp limit 2.604 mm
% theo radius at plasma exit = 0.0455
% theo radius at OTR = 0.02*sqrt(1+13.5^2/4.783^2)


% colors
navyblue = [0,0,0.502];
crimsom = [220,20,60]/256;
% initialize variables
peak_freqs = zeros(nslices,length(dump));
charge_in_slice = ones(nslices,length(dump));

% classes
AFFT = AwakeFFT('load_data_flag',0,'center_around_0_flag',0,...
    'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'xi_range',xi_range,'trans_lims',trans_lims,...
    'scan_type',scan_type,'on_axis',on_axis);

AFFT2 = AwakeFFT('load_data_flag',0,...
    'datadir',datadirs_exp{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'xi_range',xi_range,'trans_lims',trans_lims,...
    'scan_type',scan_type,'on_axis',on_axis);

P = Plotty('plots_dir',plots_dir,'plot_name',['fvsg',''],'plasmaden',plasma_density,'save_flag',save_plot_flag,...
    'datadir',datadirs{1});

EDA = ExperimentalDataAnalyser('plasmaden',plasma_density,'datadir',datadirs_exp{1});
EDA.loadSCdata();

AFFT2.dr = EDA.SCI_r(2) - EDA.SCI_r(1);

if slit_flag
    slit_text = '_slit';
else
    slit_text = '';
end

for d = 1:length(datadirs)
    
    if d == 2
        EDA.datadir = datadirs_exp{d-1}; EDA.loadSCdata();
    else
        EDA.datadir = datadirs_exp{d}; EDA.loadSCdata();
    end
    
    AFFT.datadir = datadirs{d};
    AFFT.dataformat = dataformats{d};
    
    %     timetemp = load([datadirs{d},'densitytimeprofile_slit.mat']);
    timetemp = load([datadirs{d},'/n',num2str(load_dump),'_m',num2str(measurement_distance),slit_text,'.mat']);
    simtimeprofile = timetemp.densitymatrix;
    z_sim = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2));
    r_sim = linspace(0,timetemp.ylims(2),size(simtimeprofile,1));
    AFFT.r = r_sim;
    zz_sim = z_sim*1e-12*EDA.c_cm;
    
    z_ind = zz_sim > AFFT.xi_range(2) & ... %large
        zz_sim <= AFFT.xi_range(1); % small
    AFFT.z = zz_sim(z_ind);
    AFFT.proton_beam = simtimeprofile(:,z_ind);
    AFFT.dr = r_sim(2)-r_sim(1);
    
    AFFT.fft_dataload();
    xx = AFFT.fft_densitymatrix;
    
    AFFT2.r = EDA.SCI_r;
    % trim
    z = EDA.SCI_z*1e-12*EDA.c_cm;
    
    z_ind = z > AFFT2.xi_range(2) & ... %large
        z <= AFFT2.xi_range(1); % small
    AFFT2.z = z(z_ind);
    
    AFFT2.proton_beam = EDA.SCI(:,z_ind);
    AFFT2.fft_dataload();
    prop_distance_m = AFFT.propagation_distance/100; % propagation distance in m
    
    % calculates fft and gives AFFT.fft_frequencies
    % and AFFT.fft_powerspectrum (_den or _fld)
    AFFT.get_fft(); AFFT2.get_fft();
    yy = AFFT2.fft_densitymatrix;
    switch AFFT.property
        case 'density'
            %             charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
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
            peak_freqs2(r,d) = 0;
        else
            peak_freqs(r,d) = AFFT.max_loc;
            peak_freqs2(r,d) = AFFT2.max_loc;
        end
        
        %%% Plots of the power spectrum.
        if plot_powerspectra && d < 4 %&& r == 2
            P2 = P;
            trans_lims_plot_title = [trans_lims];
            fft_freqs_plot = AFFT.fft_frequencies(AFFT.ind_pksearch)/1e9;
            amplitude_plot = powerspectrum_in_slice(r,AFFT.ind_pksearch);
            P2.plots_dir = [plots_dir,'/power_spectra'];
            P2.plot_name = ['ps',AFFT.datadir,'n',num2str(AFFT.dump),...
                '_s',num2str(r)];
            
            fig_powerspectra = figure('visible','on');
            plot(fft_freqs_plot,amplitude_plot,'LineWidth',2);
            title(['power spec. ',datadirs{d},...
                'r = [',...
                num2str(0),', ',num2str(trans_lims_plot_title(r)*10,2),']mm)']);
            ylabel('amplitude');
            xlabel('frequency (GHz)');
            xlim([fft_freqs_plot(1),fft_freqs_plot(end)]);
            drawnow;
            P2.fig_handle = fig_powerspectra;
            P2.save_plot();
%             close all;
            AFFT.progress_dump('save slice fft',r,nslices);
            
%             AFFT.dtime = 0;
%             indz = AFFT.z - AFFT.dtime > 7;
%             z_plot = AFFT.z(indz);
%             
%             density_plot = AFFT.fft_densitymatrix(r,indz);
%             P2.plots_dir = [plots_dir,'/density_distributions'];
%             P2.plot_name = ['denlong',AFFT.datadir,'n',num2str(AFFT.dump),...
%                 '_s',num2str(r)];
%             
%             figy = figure('visible','off');
%             plot(z_plot,density_plot,'LineWidth',2);
%             title(['long. den. distributon (proEDA. dist. = ',...
%                 num2str(AFFT.propagation_distance/100,2),'m; r = [',...
%                 num2str(trans_lims_plot_title(r)*10,2),', ',num2str(trans_lims_plot_title(r+1)*10,2),']mm)']);
%             ylabel('density (1/cm)');
%             xlabel('z-ct (m)');
% %             xlim([z_plot(1),z_plot(end)]);
%             figy.Units = 'normalized';
%             figy.Position = [0.05 0.2 0.9 0.5];
%             drawnow;
%             P2.fig_handle = figy;
%             P2.save_plot();
% %             close all;
%             AFFT.progress_dump('save slice density',r,nslices);
        end
        
    end % for nslices
    
    if strcmp(datadirs{d},'gm10')
        a = 1;
    end
    
end % for datadirs

peak_freqs_plot = peak_freqs;
freqs_exp = peak_freqs2(:,[1,3:9]);

%% plotty
fig_fvsrg = figure();
fig_fvsrg.Units = 'centimeters';
fig_fvsrg.Position = [1,1,8.6,8.6*3/4]*1.5;

ax_fvsrg = axes('Parent',fig_fvsrg);
ax_fvsrg.FontUnits = 'centimeters';
ax_fvsrg.FontSize = 0.4; % cm
hold on
norm_freq = 1; % AFFT.plasmafreq_GHz
marker_size = 10;
% simulation
plot_sim(1) = plot(grads_sim,... % narrow trans. range
    peak_freqs_plot(1,:)/norm_freq,'--o','markerSize',marker_size,...
    'color',crimsom,'markerfacecolor',crimsom,'Linewidth',1);
plot_sim(2) = plot(grads_sim,... % wide trans. range
    peak_freqs_plot(end,:)/norm_freq,'--o','markerSize',marker_size,...
    'color',crimsom,'Linewidth',1);

% experiment
plot_exp(1) = plot(grads_exp,... % narrow trans. range
    freqs_exp(1,:)/norm_freq,'--s','markerSize',marker_size,'LineWidth',1,...
    'markerfacecolor',navyblue,'color',navyblue);
plot_exp(2) = plot(grads_exp,... % wide trans. range
    freqs_exp(end,:)/norm_freq,'--s','markerSize',marker_size,'LineWidth',1,...
    'color',navyblue);

ax_handle = gca;
plot_fpe(1) = plot([-2.1 2.1],AFFT.plasmafreq_GHz*[1,1],'--k','LineWidth',2);
plot_fpe(2) = plot(grads_sim,AFFT.plasmafreq_GHz*[0.894394264356599,0.921920505722216,0.948648374071357,0.974643553503204,0.999963186893636,1.02465735438711,1.04877023825847,1.07234105181787,1.09540478827280],...
    'k','LineWidth',2);

frequencies_save = [grads_sim',peak_freqs_plot'];
save('simulation_frequencies.txt','frequencies_save','-ascii')

ylim(freq_range/norm_freq)
ylabel('$f_{\mathrm{mod}}$ (GHz)','interpreter','latex')

legend([plot_sim(1) plot_sim(2) plot_exp(1) plot_exp(2) plot_fpe(1) plot_fpe(2)],...
    'Simulation narrow window','Simulation wide window',...
    'Experiment narrow window','Experiment wide window',...
    'Plasma frequency entrance','Plasma frequency exit',...
    'location','southeast','AutoUpdate','off','FontSize',9);
legend('boxoff')

xlim([-2.1 2.1])
hold off
xlabel('$g$ (\%/m)','interpreter','latex')


P.save_plot(fig_fvsrg);




