%________________________________________________________________________
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
% Last update: 08/06/2020
%________________________________________________________________________

% load color order for 9 gradients
load('color_red_to_blue.mat');
cc = ccrb; %

run_it = 0;

% data directory
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};

% plot properties
leg = {'$-2$\,\%/m','$-1.5$\,\%/m','$-1$\,\%/m','$-0.5$\,\%/m',...
    '\ $0$\,\%/m','$+0.5$\,\%/m','$+1$\,\%/m','$+1.5$\,\%/m','$+2$\,\%/m'};
markers = {'o','+','x','s','d','p','h','.','*'};
grads = [-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]/100;

% parameters
plasma_density = 1.81e14;

% properties
property = 'density';
species = 'proton_beam';
field = 'e';
direction = 'z';

% simulation parameters
dump_list = 0:1:133;
dataformat = 'mat';
useAvg = 0;

% limits
xi_range = [14 0.345];

% analysis parameters
scan_type = 'cumulative'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout
low_freqrange = 0.89; %0.88; % min frequency to plot (norm. to plasma freq)
upp_freqrange = 1.11; %1.01 % max frequency to plot (norm. to plasma freq)

% calculated parameters
trans_lims = 0.018;

% initialize variables
peak_freqs = zeros(1,length(dump_list));
prop_distance_m = zeros(1,length(dump_list));
peak_freqs_den = zeros(1,length(dump_list));
peak_freqs_fld = zeros(1,length(dump_list));

AFFT = AwakeFFT('datadir',datadirs{1},...
    'plasmaden',plasma_density,'property','density','species',species,...
    'direction',direction,'wakefields_direction','long',...
    'dump',dump_list(1),'dataformat',dataformat,'useAvg',useAvg,...
    'trans_lims',trans_lims,'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis,...
    'lineout_point',3,'n_zeropadding',20);

P = Plotty('plots_dir','gradsim_paper/fft/gz/','plot_name','fvsgz');



if run_it
    for d = 1:length(datadirs)
        
        AFFT.datadir = datadirs{d};
        
        AFFT.fft_dataload();
        AFFT.get_fft();
        
        freq_ind = (AFFT.fft_frequencies*1e-9/AFFT.plasmafreq_GHz > low_freqrange) & ...
            (AFFT.fft_frequencies*1e-9/AFFT.plasmafreq_GHz < upp_freqrange);
        powerspectra_den = zeros(sum(freq_ind),length(dump_list));
        
        for n = 1:length(dump_list)
            
            AFFT.dump = dump_list(n);            
            AFFT.fft_dataload();
            prop_distance_m(n) = AFFT.propagation_distance/100; % propagation distance in m
            
            % calculates fft and gives AFFT.fft_frequencies
            % and AFFT.fft_powerspectrum (_den or _fld)
            AFFT.get_fft();
            powerspectra_den(:,n) = AFFT.fft_powerspectrum_den(freq_ind);

            AFFT.fft_peaks(AFFT.fft_powerspectrum_den);
            % max_loc is the location in freq space of the ampltiude peak in
            % the power spectrum
            if isempty(AFFT.max_loc)
                peak_freqs_den(d,n) = 0;
            else
                peak_freqs_den(d,n) = AFFT.max_loc;
            end
            
            AFFT.progress_dump('frequency',n,length(dump_list))
        end % for dump list

        peak_freqs_den_plot = peak_freqs_den;
        
    end % for datadirs
    
else
    load('gradsim_fft_fvsgz18.mat')
end

peak_freqs_den_plot(peak_freqs_den_plot/(AFFT.plasmafreq_GHz) < low_freqrange) = nan;
peak_freqs_den_plot(peak_freqs_den_plot/(AFFT.plasmafreq_GHz) > upp_freqrange) = nan;

%%
fig_fvsz = figure(2);
fig_fvsz.Units = 'centimeters';
fig_fvsz.Position = [1,1,8.6,8.6*3/4]*1.5;

ax_fvsz = gca;
ax_fvsz.FontUnits = 'centimeters';
ax_fvsz.FontSize = 0.4; % cm

hold on

for d = 1:length(datadirs)
    AFFT.plasmaden = plasma_density*(1+grads(d)*prop_distance_m(prop_distance_m < 10.2)); % !!!
    AFFT.den2freq();
    fpe_plot = AFFT.plasmafreq_GHz;
    plot(prop_distance_m(prop_distance_m < 10.2),fpe_plot,'LineWidth',2,'Color',[cc(d,:),1]) % !!!!
end % for datadirs

AFFT.plasmaden = plasma_density;
AFFT.den2freq();
plasmafreq_ini = AFFT.plasmafreq_GHz;

for d = 1:length(datadirs)
    p(d) = plot(prop_distance_m,peak_freqs_den_plot(d,:),...
        'Color',cc(d,:),'MarkerEdgeColor',cc(d,:),'MarkerFaceColor',cc(d,:),...
        'Marker',markers{d});
end % for datadirs

hold off

xlim([0 prop_distance_m(end)])
ylim([0.89,1.1]*plasmafreq_ini)

xlabel('$z$ (m)','interpreter','latex')
ylabel('$f_{\mathrm{mod}}$ (GHz)','interpreter','latex')

legend(p(5:-1:1),leg(5:-1:1),'NumColumns',1,'Box','off','location','southwest','interpreter','latex')
ah1 = axes('position',get(gca,'position'),'visible','off');
ah1.FontUnits = 'centimeters';
ah1.FontSize = 0.4; % cm
legend(ah1,p(9:-1:6),leg(9:-1:6),'NumColumns',1,'Box','off','location','northwest','interpreter','latex')

drawnow;
P.fig_handle = fig_fvsz;
P.save_plot();

if run_it
    save('gradsim_fft_fvsgz18.mat','prop_distance_m','peak_freqs_den_plot')
end


