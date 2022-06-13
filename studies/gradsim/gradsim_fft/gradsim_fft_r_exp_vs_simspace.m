%________________________________________________________________________
% Freq vs r for 1 gradient.
% Comparison of simulation and experiment. 
% Loads experimental data using ExperimentalDataLoader.
% 
% FFT of the proton density distribution.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 11/01/2021
%________________________________________________________________________
clear;

% data directory
datadirs = {'g0'};
datadir_exp = {'g0'};
% datadirs = {'gm20','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
dataformat  = 'mat';
useAvg      = false;
dump        = 133;

% save directory
plots_dir = ['gradsim/fft/r/',datadirs{1},'n',num2str(dump),''];

save_format = {'png','eps','fig'};

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';

% analysis
xi_range        = [14 0.656]; % cm % 2nd microbunch end = 0.656 cm, 1st 384
trans_lims_sim     = 0.018*(1:7);
trans_lims_exp     = 0.018*(-7:7);
nslices_sim = length(trans_lims_sim); % does not include 0
nslices_exp = length(trans_lims_exp) - 1; % includes a 0

trans_upperlimit   = 0.2; % cm


scan_type       = 'slice'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout
max_dotsize     = 150;

% switches
use_raw             = false;
plot_powerspectra   = false;
save_all_plots      = false;
save_plot_flag      = false;
showChargeinDotSize = false; % show dot size to reflect charge
normalized_frequency_switch = true;


for d = 1:length(datadirs)
    
    % classes
    AFFT = AwakeFFT(...
        'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
        'plasmaden',plasma_density,'property',property,'species',species,...
        'xi_range',xi_range,'trans_lims',trans_lims_sim,...
        'scan_type',scan_type,'on_axis',on_axis);
    
    AFFT2 = AwakeFFT('load_data_flag',0,...
        'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
        'plasmaden',plasma_density,'property',property,'species',species,...
        'xi_range',xi_range,'trans_lims',trans_lims_exp,...
        'scan_type',scan_type,'on_axis',on_axis);
    
    
    P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'save_format',save_format,...
        'datadir',datadirs{1});
    EDA = ExperimentalDataAnalyser('datadir',datadir_exp{1});
    EDA.loadSCdata();
    
    AFFT2.dr = EDA.SCI_r(2) - EDA.SCI_r(1);
    
    datadir = datadirs{d};
    plot_name = ['fvsr',datadir,'r',num2str(trans_upperlimit*100)];
    
    AFFT.datadir = datadir;
    AFFT.scan_type = scan_type;
    
    if use_raw
        AFFT.fft_rawdataload();
    else
        AFFT.fft_dataload();
        AFFT2.r = EDA.SCI_r;
        
        % trim
        z = EDA.SCI_z*1e-12*P.c_cm;
        
        z_ind = z > xi_range(2) & ... %large
            z <= xi_range(1); % small
        AFFT2.z = z(z_ind);
        
        AFFT2.proton_beam = EDA.SCI(:,z_ind);
        AFFT2.fft_dataload();
    end
    prop_distance_m = AFFT.propagation_distance/100; % propagation distance in m
    
    % calculates fft and gives AFFT.fft_frequencies
    % and AFFT.fft_powerspectrum (_den or _fld)
    AFFT.get_fft(); AFFT2.get_fft();
    
    switch AFFT.property
        case 'density'
            charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
            data_in_slice = AFFT.fft_powerspectrum_den;
            data_in_slice2 = AFFT2.fft_powerspectrum_den;
        case 'fields'
            data_in_slice = AFFT.fft_powerspectrum_fld;
    end % switch property
    
    % initialize variables
    peak_freqs = zeros(nslices_sim,1);
    peak_freqs2 = zeros(nslices_exp,1);
    for r = 1:nslices_sim
        
        AFFT.fft_peaks(data_in_slice(r,:));
        % max_loc is the location in freq space of the ampltiude peak in
        % the power spectrum
        if isempty(AFFT.max_loc)
            peak_freqs(r) = nan;
        else
            peak_freqs(r) = AFFT.max_loc;
        end
        
    end % for nslices
    
    for r = 1:nslices_exp
        AFFT2.fft_peaks(data_in_slice2(r,:));
        peak_freqs2(r) = AFFT2.max_loc;
    end % for nslices
    
    % Calculate DFT for whole range
    AFFT.scan_type = 'cumulative';
    AFFT.on_axis = 'int';
    AFFT.trans_lims = trans_upperlimit;
    if use_raw
        AFFT.fft_rawdataload();
    else
        AFFT.fft_dataload();
    end
    AFFT.get_fft();
    AFFT.fft_peaks(AFFT.fft_powerspectrum_den);
    freq_widerange = AFFT.max_loc;
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize == 0) = nan;
    else
        dotsize = max_dotsize*ones(length(charge_in_slice),1);
        dotsize(dotsize == 0) = nan;
    end % end if charge dot size

    peak_plot = peak_freqs2;
    
    %     hold on
    
    %     if normalized_frequency_switch
    %         scatter(trans_lims_mm,...
    %             peak_freqs/AFFT.plasmafreq_GHz,dotsize,'filled');
    %         scatter(trans_lims_mm,peak_freqs2/AFFT.plasmafreq_GHz,dotsize,'x','LineWidth',1); % experiment
    %         plot(xlimits,freq_widerange*ones(1,2)/AFFT.plasmafreq_GHz,'--b','LineWidth',2);
    %         freq_lims = [0.98*min([peak_freqs2';peak_freqs]/AFFT.plasmafreq_GHz),...
    %             1.02*max([peak_freqs2';peak_freqs]/AFFT.plasmafreq_GHz)];
    %         ylim(freq_lims)
    %         ylabel({'transverse slice frequency /','plasma frequency at z = 0 m'})
    %         set(gca,'xdir','reverse')
    %     else
    %         scatter(trans_lims_mm,peak_freqs,dotsize,'filled'); % simulation
    %         scatter(exp_translims,exp_freqs,dotsize,'x','LineWidth',1); % experiment
    %         plot(xlimits,freq_widerange*ones(1,2),'--b','LineWidth',2);
    %         freq_lims = [0.98*min([peak_freqs2';peak_freqs]),1.02*max([peak_freqs2';peak_freqs])];
    %         ylim(peak_freqs);
    %         ylabel('freq. (GHz)')
    %     end % if normalized frequency switch
    %
    %     hold off
    %     grid on
    %     xlim(xlimits)
    %     xlabel('r (mm)')
    %     legend('simulation','experiment','location','best');
    %     drawnow;
    %
    %     P.plot_name = plot_name;
    %     P.fig_handle = fig_fvsr;
    %     P.save_plot();
    
    
    
end % for datadirs

fig_fvsr = figure;

% plot limits
trans_lims_mm = (10*trans_lims_sim)'; % trans lims in mm
xlimits = [-1.01*trans_lims_mm(end) 1.01*trans_lims_mm(end)];

hold on
scatter([flipud(-trans_lims_mm);trans_lims_mm],...
    [flipud(peak_freqs);peak_freqs]/AFFT.plasmafreq_GHz,[dotsize;dotsize],'filled');

scatter([flipud(-trans_lims_mm);trans_lims_mm],peak_plot/AFFT.plasmafreq_GHz,...
    [dotsize;dotsize],'x','LineWidth',1); % experiment
plot(xlimits,freq_widerange*ones(1,2)/AFFT.plasmafreq_GHz,'--b','LineWidth',2);

hold off

freq_lims = [0.98*min(peak_freqs/AFFT.plasmafreq_GHz),...
    1.02*max(peak_freqs/AFFT.plasmafreq_GHz)];
ylim(freq_lims); xlim(xlimits);
ylabel({'normalized transverse slice frequency (a.u.)'})
xlabel('r (mm)')
legend('simulation','experiment','location','best');
P.plot_name = plot_name;
P.fig_handle = fig_fvsr;
P.save_plot();
