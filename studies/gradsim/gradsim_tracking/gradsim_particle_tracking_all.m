%__________________________________________________________________________
% Script that takes tracks according to some conditions and plots them.
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 08/07/2020
%__________________________________________________________________________
% comment: particle tracks are saved each 5 dumps, they start at dump 0,
% the second save is dump 4, third save dump 9 and so on. So the last point
% for the gradients which has 200 points corresponds to dump 99.5. To have
% corresponding values, one should use dump 99 with point 199.

% clear;
% close all;
% data location
datadirs = {'gm20'};
dataformat = 'mat';
useAvg = false;

dump_fft = 129; % must be at 99

% simulation parameters
plasmaden = 1.81e14; %cm^-3
species = 'proton_beam';

trans_ranges = [0,0.018*(1:8)]' + [0.0 0.01]; %(0:0.01:0.14)' + [0.0 0.01]; % cm

% FFT parameters
scan_type = 'slice';
on_axis = 'sum';
xi_range = [14 0.38];
nslices = size(trans_ranges,1);

% plotting
P = Plotty('plasmaden',plasmaden);
color_meanline = [
     0.1900    0.0718    0.2322
    0.1900    0.0718    0.2322
    0.1900    0.0718    0.2322
    0.1567    0.7378    0.9208
    0.6541    0.0793    0.0041
    0.5588    0.0404    0.0058
    0.4796    0.0158    0.0106
    0.4796    0.0158    0.0106
    0.4796    0.0158    0.0106];

datadir = datadirs{1};

OPT = OsirisParticleTracking('datadir',datadir,'plasmaden',plasmaden,...
    'dataformat','mat',...
    'property','tracks','trackfile_suffix','1',...
    'track_dataset','q');
OPT.getdata();

if strcmp(OPT.dataformat,'mat')
    OPT.track_dataset = 'x'; OPT.direction = 'z'; OPT.getdata();
    OPT.track_dataset = 'x'; OPT.direction = 'r'; OPT.getdata();
end

fprintf('loading done\n');

P.plots_dir = ['gradsim_paper/tracking/',datadir];
peak_freqs_t = zeros(nslices,1);
peak_amplitude = zeros(nslices,1);
peak_phase = zeros(nslices,1);

for rr = 1:length(trans_ranges)
    trans_range = trans_ranges(rr,:);
    trans_lim = trans_range;
    OPT.trans_range = trans_range;
    
    on_axis = 'int';
    AFFT = AwakeFFT('datadir',datadir,...
        'plasmaden',plasmaden,'property','density','species',species,...
        'dump',dump_fft,'dataformat',dataformat,'useAvg',useAvg,...
        'trans_lims',trans_lim,'xi_range',xi_range,...
        'scan_type',scan_type,'on_axis',on_axis);
    
    
    AFFT.fft_dataload();
    
    % calculates fft and gives AFFT.fft_frequencies
    % and AFFT.fft_powerspectrum (_den or _fld)
    AFFT.get_fft();
    
    charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
    data_in_slice = AFFT.fft_powerspectrum_den;
    
    
    
    tracks_r = OPT.denorm_distance(OPT.tracks_r(:,[1:250,301:end]));
    tracks_z = OPT.denorm_distance(OPT.tracks_z(:,[1:250,301:end]));
    
    AFFT.fft_peaks(data_in_slice(end,:),AFFT.fft_phase_den(end,:),AFFT.fft_densitymatrix(end,:));
    % max_loc is the location in freq space of the ampltiude peak in
    % the power spectrum
    if ~isempty(AFFT.max_loc)
        peak_freqs_t(rr,1) = AFFT.max_loc;
        peak_amplitude(rr,1) = AFFT.max_peak;
        peak_phase = AFFT.max_phase;
    end
    
    wavenumber = peak_freqs_t(rr,1)*1e9/AFFT.c_cm;
    
    cos_signal = cos(2*pi*AFFT.z*wavenumber + peak_phase);
    z_cos_positive = AFFT.z(cos_signal > 0);
    
    ind_r = ((tracks_r(:,end-1) < trans_range(2)) ...
        & (tracks_r(:,end-1) > trans_range(1)));
    
    ind_morethanzero = OPT.tracks_z(:,end-1) > 0;
    ind_z = [(tracks_z(:,end-1)-min(tracks_z(ind_morethanzero,end-1)) > xi_range(2)) ...
        & (tracks_z(:,end-1)-min(tracks_z(ind_morethanzero,end-1) < xi_range(1)))];
    
    ind_rz = find(ind_r & ind_z);
    
    
    par_last_z = tracks_z(ind_rz,end-1);
    %         par_last_r = OPT.tracks_r(ind_rz,end-1);
    diff_matrix = par_last_z - z_cos_positive;
    
    ind_diff = logical(sum(abs(diff_matrix) <= AFFT.dz/2,2)); %indices for particles which are on the positive zones of cosine
    ind_final = ind_rz(ind_diff);
    
    
    par_z = tracks_z(ind_final,1:end-1); % z position for particles in the positive region of the cosine from the fft
    par_r = tracks_r(ind_final,1:end-1); % r position for particles in the positive region of the cosine from the fft
    
    par_q = OPT.tracks_q(ind_final,1);
    
    on_axis = 'int';
    switch on_axis
        case 'int'
            par_q2 = par_q;
            par_rr = par_r;
        case 'sum'
            par_rr = par_r;
            if rr > 3
                ind_rlow = par_r < 10*4e-4;
            else
                ind_rlow = par_r < 2*4e-4;
            end
            par_rr(ind_rlow) = nan;
            
            par_q2 = par_q./par_rr; %par_q
            
    end
    
    
    
    rave(rr,:) = 10*sum((par_q2.*par_rr),'omitnan')./sum(par_q2,'omitnan');
    zave(rr,:) = sum((par_q2.*par_z),'omitnan')./sum(par_q2,'omitnan')/100;
    
    
    OPT.progress_dump('slice',rr,length(trans_ranges))
    
end %% rr trans_ranges

%% new
powerspectra = zeros(dump_fft,9,313);
for p = 1:dump_fft+1
    
    fft_filename  = ['save_files/fft/rz/freqpowerspec_',AFFT.datadir,'n',num2str(p),...
        'xi',num2str(xi_range(1),2),...
        'xi',num2str(xi_range(2),2),'nslice',num2str(9),'sum','.mat'];
    fft_temp = load(fft_filename);
    if ~isempty(fft_temp.frequency)
        fft_frequencies = fft_temp.frequency;
        powerspectra(p,:,:) = fft_temp.powerspectrum;
    else
        powerspectra(p,:,:) = 0.001;
    end
    
    
    z_filename  = ['save_files/fft/rz/zlongprofile_',AFFT.datadir,'n',num2str(p),...
        'xi',num2str(xi_range(1),2),...
        'xi',num2str(xi_range(2),2),'nslice',num2str(9),'sum','.mat'];
    charge = load(z_filename,'-mat','charge');
    chargematrix(:,p) = charge.charge;
    OPT.progress_dump('load data',p,dump_fft+1)
end


%% powerspectra
% powerspectra ( dump, radius, power spectrum  )
for rr = 1:length(trans_ranges)
    
    ind_freq = (fft_frequencies > 0.998*peak_freqs_t(rr)) & (fft_frequencies < 1.002*peak_freqs_t(rr));
    
    r_position_tracks = ceil(rave(rr,:)/0.2); % look here
    z_position_tracks = repelem(1:130,2); z_position_tracks(end) = [];
    
    AFFT.fft_frequencies = fft_frequencies*1e9;
    powerspectramatrix = nan(size(chargematrix));
    for p = 1:dump_fft+1
        for rrr = min(r_position_tracks):max(r_position_tracks)
            AFFT.fft_peaks(squeeze(powerspectra(p,rrr,:)));
            ind_pks = (AFFT.locs > 0.995*peak_freqs_t(rr)) & (AFFT.locs < 1.005*peak_freqs_t(rr));
%             [maxpk,maxpk_loc] = max(AFFT.pks);
%             ind_pks = (AFFT.locs(maxpk_loc) > 0.99*peak_freqs_t(rr)) & (AFFT.locs(maxpk_loc) < 1.01*peak_freqs_t(rr));

            if any(ind_pks)
                [maxamplitude,max_loc] = max(AFFT.pks(ind_pks));
%                 maxamplitude = maxpk;
                powerspectramatrix(rrr,p) = maxamplitude;
                
                
                if sum(ind_pks) > 1
                    %                 pause
                end
            end
        end
    end
    
    %     powerspectramatrix = sum(powerspectra(:,:,ind_freq),3)';
    
    ind_charge_r = sub2ind(size(chargematrix), r_position_tracks, z_position_tracks);
    
    charge_position_r = chargematrix(ind_charge_r);
    integral_around_peak = powerspectramatrix(ind_charge_r);
    
    linewidths_temp = integral_around_peak./charge_position_r;
    %     linewidths(rr,:) = 5*linewidths_temp/max(linewidths_temp,[],'all');
    linewidths(rr,:) = linewidths_temp;
end

linewidths = 6*linewidths/max(linewidths,[],'all');
linewidths_isnan = isnan(linewidths);
linewidths(linewidths_isnan) = 1e-10;
% linewidths(~linewidths_isnan) = 4; %1e-10;

ind_nan = linewidths < 1.5;
rave_plot = rave;
rave_plot(ind_nan) = nan;
zave_plot = zave;
zave_plot(ind_nan) = nan;

disp('done')


%% plot
fig1 = figure;

% hold on
% plot(zave',rave','color',[0.7 0.7 0.7],'LineWidth',3);
% hold off


for rr = length(trans_ranges):-1:1
    
    for p = 0:dump_fft-1
        hold on
        pave = plot(zave_plot(rr,2*(p+0.5):2*(p+0.5)+2), rave_plot(rr,2*(p+0.5):2*(p+0.5)+2),...
            'color',color_meanline(rr,:),'LineWidth',linewidths(rr,2*p+1)); %(rr,2*p+1)
        hold off
        OPT.progress_dump('plotillo',p,100)
    end
end

fig1.CurrentAxes.FontSize = 15;
xlabel('z (m)');
ylabel('r (mm)');

xlim([0 13.3])
ylim([0 1.5])
fig1.Units = 'normalized';
fig1.Position = [0.2 0.2 0.7 0.5];


yticks(0:0.18:1.5)
grid on

P.plot_name = 'tracksblank';
P.fig_handle = fig1;
P.save_plot();





