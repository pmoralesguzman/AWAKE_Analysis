%________________________________________________________________________
% freq in one slice in r, for all z, for all gradients, in one plot
% The frequency is obtained from the distance between maximum points of the
% microbunch density profile on axis or the zero-crossings of the fields.

% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 26/03/2021
%________________________________________________________________________

clear;
close all;

plots_dir = ['gradsim_test/freq/'];

% load color order for 9 gradients
load('color_red_to_blue.mat'); % ccrb

% cell plotting parameters
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp1','gp2','gp5','gp10'};
% leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
% line_style = {':','--','-.','-','-','-','-.','--',':'};
% colors = {'r','k',[0 0.5 0],'b',[0.5,0,0.5]};

datadirs = {'gm20'};
leg = {'0.5 %/m'};
line_style = {'-'};
colors = {'k'};

% study parameters
plasmaden = 1.81e14;
dump_list = 10:1:134;
useAvg = false;
dataformat = 'mat';
trans_range = [0 0.018]; %cm
xi_range = [4 0.38];

%
property        = 'density';
scan_type       = 'slice'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout

%


% initialize classes
P = Plotty('plasmaden',plasmaden,'plots_dir',plots_dir);
AF = AwakeFFT('datadir',datadirs{1},...
    'property',property,'wakefields_direction','long',...
    'plasmaden',plasmaden,...
    'useAvg',useAvg,...
    'dataformat',dataformat,...
    'trans_range',trans_range,'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis);

% dephasing_xi = [0.9:OPA.plasma_wavelength:15]; % cm
% dephasing_xi0_plot = dephasing_xi(1:end-1) + OPA.plasma_wavelength/2;
% dephasing_lines_number = 1:(length(dephasing_xi) - 1);

% OPA.dephasing();

% initialize variables
% dephasing_z = zeros(length(OPA.dephasing_line),1);
% dephasing_lines = zeros(length(dephasing_xi),length(datadirs),length(OPA.dephasing_line));
% dephasing_lines = zeros(length(dephasing_xi),length(OPA.dephasing_line));

%% start the script
% fig_dephase = figure(1);
% fig_dephase.OuterPosition = [100 100 1200 400];
% colororder(ccrb);
% tt = tiledlayout(1,3);
% tt.TileSpacing = 'compact';
% tt.Padding = 'compact';

for d = 1:length(datadirs)
    
    AF.datadir = datadirs{d};
    
    for n = 1:length(dump_list)
        AF.dump = dump_list(n);
        
        AF.fft_dataload(); %AF.assign_density(); AF.denorm_density();
        study_window = min(abs(diff(xi_range)),AF.simulation_window);
        den_profile = smooth(AF.fft_densitymatrix,AF.plasma_wavelength/(study_window*2),'loess');
        %         den_profile = AF.fft_densitymatrix;
        peak_dis = 0.85*AF.plasma_wavelength;
        [pks,locs] = findpeaks(den_profile,AF.z,...
            'MinPeakDistance',peak_dis);
        large_pks_ind = pks > 0.05*max(pks);
        %         locs = locs(large_pks_ind);
        locs_save(n,:) = -locs(end-10:end) + AF.dtime + AF.simulation_window;
        %         pks = pks(large_pks_ind);
        figure(1)
        plot(AF.z,den_profile)
        hold on
        scatter(locs,pks)
        hold off
        xlim([min(AF.z),max(AF.z)])
        ylim([0 max(pks)*1.1]);
        AF.progress_dump('frequency through peaks',n,length(dump_list))
        pause(0.25)
        
    end % for dump list
    %     nfreq = freq/AF.plasmafreq*2*pi;
    
    AF.progress_dump('DATADIR',d,length(datadirs))
end % datadirs
%% plotty
locs_diff = diff(locs_save);

plot(linspace(0,13.5,length(locs_save)),locs_save(:,[5:11]) - locs_save(1,[5:11]));

legend('microbunch 8','microbunch 7','microbunch 6','microbunch 5','microbunch 4','microbunch 3','microbunch 2',...
    'location','southeast')
xlim([0 13.5])

xlabel('z (m)');
ylabel('microbunch position in \xi - initial position in \xi (cm)');

P = Plotty();
P.plots_dir = 'gradsim/microbunch_location';
P.plot_name = 'gm20';
P.save_plot(gcf);

