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
% Last update: 05/11/2021
%________________________________________________________________________

% load color order for 9 gradients
load('color_red_to_blue.mat');
cc = ccrb; %

datadirs = {'DWdc3_lcode_x666'};
% datadirs = {'DW_lcode_nofront'};

run_it = 1;

% data directory
line_style = {':','--','-.','-','-','-','-.','--',':'};

% plot properties
leg = {'baseline','density feature pos. 1','density feature pos. 2'};

markers = {'o','+','x','s','d','p','h','.','*'};
% grads = [-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]/100;

% parameters
plasma_density = 2e14;

% properties
property = 'fields';
species = 'proton_beam';
field = 'e';
direction = 'z';

% simulation parameters
dump_list = 200:1:200;
dataformat = 'mat';
useAvg = 0;

% limits
xi_range = [7.5 0.0];
xi_range = [15 7.5];

% analysis parameters
scan_type = 'cumulative'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout
low_freqrange = 0.89; %0.88; % min frequency to plot (norm. to plasma freq)
upp_freqrange = 1.11; %1.01 % max frequency to plot (norm. to plasma freq)

% calculated parameters
trans_lims = 0.02;

% initialize variables
peak_freqs = zeros(1,length(dump_list));
prop_distance_m = zeros(1,length(dump_list));
peak_freqs_den = zeros(1,length(dump_list));
peak_freqs_fld = zeros(1,length(dump_list));

AFFT = AwakeFFT('datadir',datadirs{1},...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'direction',direction,'wakefields_direction','long',...
    'dump',dump_list(1),'dataformat',dataformat,'useAvg',useAvg,...
    'trans_lims',trans_lims,'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis,...
    'lineout_point',3,'n_zeropadding',20);

P = Plotty('plasmaden',plasma_density,'plots_dir','DW/fft/','plot_name','fvsgz');



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
            powerspectra_den(:,n) = AFFT.fft_powerspectrum_fld(freq_ind);

            AFFT.fft_peaks(AFFT.fft_powerspectrum_fld);
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
    %load('gradsim_fft_fvsgz18.mat')
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

legend(p,leg,'NumColumns',1,'Box','off','location','southwest','interpreter','latex')
ah1 = axes('position',get(gca,'position'),'visible','off');
ah1.FontUnits = 'centimeters';
ah1.FontSize = 0.4; % cm

drawnow;
P.fig_handle = fig_fvsz;
P.save_plot();

if run_it
    save('dw_fft_fvsz1.mat','prop_distance_m','peak_freqs_den_plot')
end


