%________________________________________________________________________
% Freq vs r for 1 gradient.
% Comparison of simulation and experiment.
% FFT of the proton density distribution.
% Version for paper.
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
datadir     = 'gm20'; % gm20, g0, gp20
dataformat  = 'mat';
useAvg      = 0;
dump        = 133;
load_dump   = 130;
measurement_distance = 1350;
slit_flag = 1;

% save directory
plots_dir = ['gradsim_paper/fft/r/','n',num2str(dump),''];

save_format = {'png','eps','fig'};

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';

% analysis
xi_range        = [14 0.374]; % cm % 2nd microbunch end = 0.656 cm, 1st 384
r_step          = 0.018;
n_r_step        = 9;
trans_lims_sim  = r_step*(1:n_r_step); % cm
trans_lims_exp  = r_step*(-n_r_step:n_r_step); % cm
nslices_sim     = length(trans_lims_sim); % does not include 0
nslices_exp     = length(trans_lims_exp) - 1; % includes a 0

scan_type       = 'slice'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout
max_dotsize     = 150;

% switches
save_plot_flag      = 0;
showChargeinDotSize = false; % show dot size to reflect charge

% classes
% AwakeFFT for the simulation
AFFT = AwakeFFT('load_data_flag',0,'center_around_0_flag',0,...
    'datadir',datadir,'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'xi_range',xi_range,'trans_lims',trans_lims_sim,...
    'scan_type',scan_type,'on_axis',on_axis);

% AwakeFFT for the experimental data
AFFT2 = AwakeFFT('load_data_flag',0,'center_around_0_flag',0,...
    'datadir',datadir,'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'xi_range',xi_range,'trans_lims',trans_lims_exp,...
    'scan_type',scan_type,'on_axis',on_axis,'center_around_0_points',151);

P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'save_format',save_format,...
    'datadir',datadir,'save_flag',save_plot_flag);
plot_name = ['fvsr',datadir,'m',num2str(measurement_distance)];

EDA = ExperimentalDataAnalyser('datadir',datadir,'plasmaden',plasma_density);

if slit_flag
    slit_text = '_slit';
else
    slit_text = '';
end

timetemp = load([datadir,'/n',num2str(load_dump),'_m',num2str(measurement_distance),slit_text,'.mat']);
simtimeprofile = timetemp.densitymatrix;
z_sim = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2))-2.6;
r_sim = linspace(0,timetemp.ylims(2),size(simtimeprofile,1));
AFFT.r = r_sim;
zz_sim = z_sim*1e-12*EDA.c_cm;

z_ind = zz_sim > xi_range(2) & ... %large
    zz_sim <= xi_range(1); % small
AFFT.z = zz_sim(z_ind);
AFFT.proton_beam = simtimeprofile(:,z_ind);
AFFT.dr = r_sim(2)-r_sim(1);
AFFT.dz = zz_sim(2)-zz_sim(1);

% load simulation data and build densitymatrix
AFFT.fft_dataload();
xx = AFFT.fft_densitymatrix;

% load experimental data
EDA.loadSCdata();
% set dr in AFFT2 to build densitymatrix
AFFT2.dr = EDA.SCI_r(2) - EDA.SCI_r(1);
AFFT2.r = EDA.SCI_r;
% trim
z_exp = EDA.SCI_z*1e-12*P.c_cm;
z_ind = z_exp > xi_range(2) & ... %large
    z_exp <= xi_range(1); % small
AFFT2.z = z_exp(z_ind);

% set experimental data to proton_beam for analysis
AFFT2.proton_beam = EDA.SCI(:,z_ind);

% build densitymatrix for the experimental data
AFFT2.fft_dataload();

% calculates fft and gives AFFT.fft_frequencies and AFFT.fft_powerspectrum
AFFT.get_fft(); AFFT2.get_fft();

charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
exp_charge_in_slice = sum(AFFT2.fft_densitymatrix,2);
sim_densitymatrix = AFFT.fft_densitymatrix; % for checking
exp_densitymatrix = AFFT2.fft_densitymatrix; % for checking

simdata_in_slice = AFFT.fft_powerspectrum_den;
expdata_in_slice = AFFT2.fft_powerspectrum_den;

% initialize variables
simpeak_freqs = zeros(nslices_sim,1);
exppeak_freqs = zeros(nslices_exp,1);
for r = 1:nslices_sim
    AFFT.fft_peaks(simdata_in_slice(r,:));
    if isempty(AFFT.max_loc)
        simpeak_freqs(r) = 0;
    else
        simpeak_freqs(r) = AFFT.max_loc;
    end
end % for nslices sim

for r = 1:nslices_exp
    AFFT2.fft_peaks(expdata_in_slice(r,:));
    exppeak_freqs(r) = AFFT2.max_loc;
end % for nslices exp

dotsize = max_dotsize*ones(length(charge_in_slice),1);

fig_fvsr = figure();
fig_fvsr.OuterPosition = [336.1111  444.5556  862.2222  409.3333];
ax_fvsr = axes('Parent',fig_fvsr);
ax_fvsr.FontSize = 12;
% plot limits
trans_lims_mm = (10*trans_lims_sim)'; % trans lims in mm
xlimits = [-1.01*trans_lims_mm(end) 1.01*trans_lims_mm(end)];


hold on
s1 = scatter([flipud(-trans_lims_mm);trans_lims_mm],...
    [flipud(simpeak_freqs);simpeak_freqs],[dotsize;dotsize],'filled');

s2 = scatter([flipud(-trans_lims_mm);trans_lims_mm],exppeak_freqs,...
    [dotsize;dotsize],'x','LineWidth',1); % experiment
y1 = plot([-1.85,1.85],[AFFT.plasmafreq_GHz,AFFT.plasmafreq_GHz],'--r','LineWidth',2);
switch datadir
    case 'gp20'
        y2 = plot([-1.85,1.85],[AFFT.plasmafreq_GHz,AFFT.plasmafreq_GHz]*sqrt(1.2),'--b','LineWidth',2);
    case 'gm20'
        y2 = plot([-1.85,1.85],[AFFT.plasmafreq_GHz,AFFT.plasmafreq_GHz]*sqrt(0.8),'--b','LineWidth',2);  
end

hold off

freq_lims = [0.98*min(simpeak_freqs),1.02*max(simpeak_freqs)];

ylim(freq_lims); xlim(xlimits);
ylim([0.88 1.0762]*AFFT.plasmafreq_GHz);
ylabel({'transverse slice','frequency (GHz)'})
xlabel('x (mm)')
legend([s1 s2],'simulation','experiment','location','northeast');
% legend([y1 y2],'f_{pe} entrance','f_{pe} exit','experiment','location','northwest');


P.plot_name = plot_name;
P.fig_handle = fig_fvsr;
P.save_plot();
