%________________________________________________________________________
% freq in one slice in r, for all z, for all gradients
% also plot and save the powerspectra and a waterfall plot with all ps vs z
% FFT of the proton or plasma electrons density distribution, or of a
% lineout of the wakefields.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 02/07/2021
%________________________________________________________________________

% load files
load('color_red_to_blue.mat');
cc = ccrb;
ccn = 1; % 1 red, 9 blue

% data directory
datadirs = {'gm20'};
leg = {'-2 %/m'};
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
dataformat  = 'mat';
useAvg      = 0;
dump_list   = 0:1:134;

save_format = {'png','eps','fig'}; % {'png','eps','fig'}

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';
field           = 'e';
direction       = 'z';
wakefields_direction = 'long';

% analysis
xi_range        = [14 0.38];
% trans_range is set by trans_lims in this analysis
trans_lims       = [0 0.018]; % cm

scan_type       = 'slice'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout
nfreq_range     = [0.7 1.3]; % maybe norm. to plasma freq
plot_nfreq_range= [0.85 1.04]; %[0.85 1.04];%[0.85 1.04]; ; % maybe norm. to plasma freq [0.95 1.09]

% switches
plot_powerspectra = 1;
plot_ps_waterfall_flag = 1;
normalize_waterfall = 1;

save_plot_powerspectra = 0;
% save_powerspectra = 0;

save_plot_flag = 0; % general

AFFT = AwakeFFT('datadir',datadirs{1},...
    'dataformat',dataformat,'useAvg',useAvg,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'field',field,'direction',direction,'wakefields_direction',wakefields_direction,...
    'xi_range',xi_range,'trans_lims',trans_lims,...
    'scan_type',scan_type,'on_axis',on_axis,...
    'n_zeropadding',10); 

P = Plotty('plasmaden',plasma_density,'save_format',save_format,'save_flag',save_plot_flag);

% load data of first dump to calculate fft array size
AFFT.dump = dump_list(1); AFFT.fft_dataload(); AFFT.get_fft();
% calculated parameters
low_freqrange   = 1e9*AFFT.plasmafreq_GHz*nfreq_range(1); % min frequency to plot (norm. to plasma freq)
upp_freqrange   = 1e9*AFFT.plasmafreq_GHz*nfreq_range(2); % max frequency to plot (norm. to plasma freq)
freq_ind = (AFFT.fft_frequencies > low_freqrange) & (AFFT.fft_frequencies < upp_freqrange);
plot_freq = AFFT.fft_frequencies(freq_ind)/1e9;

% initialize variables
peak_freqs = zeros(length(datadirs),length(dump_list));
prop_distance_m = zeros(length(datadirs),length(dump_list));
plot_ps = zeros(sum(freq_ind),length(dump_list));

for d = 1:length(datadirs)
    AFFT.datadir = datadirs{d};
    
    switch AFFT.datadir
        case {'gp20','gp20d2','gp20_20m'}
            title_g = 'g = +2 %/m';
            grad_sim = 0.02;
        case 'gp15'
            title_g = 'g = +1.5 %/m';
            grad_sim = 0.015;
        case 'gp10'
            title_g = 'g = +1 %/m';
            grad_sim = 0.01;
        case 'gp5'
            title_g = 'g = +0.5 %/m';
            grad_sim = 0.005;
        case {'g0','g0d2'}
            title_g = 'g = 0 %/m';
            grad_sim = 0.0;
        case 'gm5'
            title_g = 'g = -0.5 %/m';
            grad_sim = -0.005;
        case 'gm10'
            title_g = 'g = -1 %/m';
            grad_sim = -0.01;
        case 'gm15'
            title_g = 'g = -1.5 %/m';
            grad_sim = -0.015;
        case 'gm20'
            title_g = 'g = -2 %/m';
            grad_sim = -0.02;
        case 'R2gap0_2e14'
            title_g = 'density step';
    end % switch datadir
    
    for n = 1:length(dump_list)
        
        AFFT.dump = dump_list(n);
        
        AFFT.fft_dataload();
        prop_distance_m(d,n) = AFFT.propagation_distance/100; % propagation distance in m
        
        % calculates fft and gives AFFT.fft_frequencies
        % and AFFT.fft_powerspectrum (_den or _fld)
        AFFT.get_fft();
        
        switch AFFT.property
            case 'density'
                powerspectrum_in_slice = AFFT.fft_powerspectrum_den(1,:);
            case 'fields'
                powerspectrum_in_slice = AFFT.fft_powerspectrum_fld;
        end % switch property
        
        plot_ps(:,n) = powerspectrum_in_slice(freq_ind);
        
        if plot_powerspectra
            figure(1000); %#ok<*UNRCH>
            plot(plot_freq,plot_ps(:,n),'k');
            
            xlabel('frequencies (GHz)');
            switch AFFT.property
                case 'fields'
                    ylabel('amplitude');
                case 'density'
                    ylabel('amplitude');
            end % switch property
            title(['power spectrum ',title_g,'; z = ',num2str(prop_distance_m(n),3),...
                ' m; r = [',num2str(trans_lims(1)*10),' ',num2str(trans_lims(end)*10),'] mm']);
            xlim(nfreq_range*AFFT.plasmafreq_GHz)
            ylim([0 1.05*max(plot_ps(:,n))])
            
            drawnow;
            P.save_flag = save_plot_powerspectra;
            P.plots_dir = ['gradsim/fft/z/',AFFT.datadir,'/ps/'];
            P.plot_name = ['r',num2str(trans_lims(end)*100),'n',num2str(n)];
            P.save_plot(gcf);
        end % if plot_powerspectra
        
        AFFT.fft_peaks(powerspectrum_in_slice);
        % maxloc is the location in freq space of the ampltiude peak in
        % the power spectrum
        if isempty(AFFT.max_loc)
            peak_freqs(d,n) = 0;
        else
            peak_freqs(d,n) = AFFT.max_loc;
        end
        AFFT.progress_dump(['fft z ',AFFT.datadir],n,length(dump_list));
    end %  for dump list
    
    
end % for datadirs

peak_freqs_plot = peak_freqs;
peak_freqs_plot(peak_freqs < low_freqrange/1e9) = nan;
peak_freqs_plot(peak_freqs > upp_freqrange/1e9) = nan;

%% Plotty
fig_gz = figure(4);
colororder(cc);
for d = 1:length(datadirs)
    hold on
    scatter(prop_distance_m(d,:),peak_freqs_plot(d,:),'filled')
    hold off
end % for datadirs

plasmafreq_ini = AFFT.plasmafreq_GHz;
AFFT.plasmaden = plasma_density*(1+grad_sim*prop_distance_m); % !!!!!!!!!!!!!!!!!!!!
AFFT.den2freq();
fpe_plot = AFFT.plasmafreq_GHz;
hold on
plot(prop_distance_m(prop_distance_m<10),fpe_plot(prop_distance_m<10)) % !!!!!!!!!!!!!!!!!!!!
hold off
AFFT.plasmaden = plasma_density;
AFFT.den2freq();

grid on
ylim(plot_nfreq_range*AFFT.plasmafreq_GHz)
xlim([0 13.5])
xlim([0 20])
xlabel('z (m)')
ylabel({'modultion freq. (GHz)'});
% title(['trans. lim. = 0.1 mm'])
title(['',title_g,' r = [',num2str(trans_lims(1)*10),' ',num2str(trans_lims(end)*10),'] mm']);
legend(leg,'Location','southwest')

P.plots_dir = ['gradsim/fft/z/',AFFT.datadir];
P.plot_name = ['fvsz_','smoothn',num2str(trans_lims(end)*100)];
P.save_plot(gcf);

if plot_ps_waterfall_flag
    figure % waterfall fields
    if normalize_waterfall
        max_ps = max(plot_ps);
    else
        max_ps = 1;
    end
    imagesc(prop_distance_m,plot_freq,plot_ps./max_ps);
%     ylim(plot_nfreq_range*AFFT.plasmafreq_GHz);
    colormap(flipud(gray)); % ------------------------ colormap
%     title([waterfall_plot_title,' ',title_g,' r = [',num2str(trans_lims(1)*10),' ',num2str(trans_lims(end)*10),'] mm']);
    ylabel('frequency (GHz)')
    xlabel('z (m)')
    set(gca,'YDir','normal')
%     colorbar;
    hold on
    plot(prop_distance_m(d,:),peak_freqs_plot(d,:),'linestyle','--','color',cc(ccn,:))
    plot(prop_distance_m(prop_distance_m<10),fpe_plot(prop_distance_m<10),'color',cc(ccn,:)) % !!!!!!!!!!!!!!!!!!!!
    hold off
    legend('maximum amplitude','$f_{pe}(z)$','interpreter','latex','location','southwest','box','off')

    ylim([106 132])
    P.plots_dir = ['gradsim_paper/fft/z/',AFFT.datadir];
    P.plot_name = ['fvszwaterfallps','n','18'];
    P.save_plot(gcf);
end








