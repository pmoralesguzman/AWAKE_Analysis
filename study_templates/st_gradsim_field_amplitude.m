%________________________________________________________________________
% Calculate and plot the mean or maximum fields 
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 15/03/2021
%________________________________________________________________________

clear;
% close all;

load('color_red_to_blue.mat');
cc = ccrb;


datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};

% save file name
plots_dir = ['gradsim/amplitude/'];
plot_name = 'standard';
saveplot = true;
colornumber = 5;
y_label = ['max. E_z (MV/m) around \xi_0 = 7 cm'];
y_label_position = ['max. E_z position ([z - z_c]/\sigma_z)'];
y_lim = [0 350];

% simulation parameters
plasmaden = 1.81e14;
dump_list = 0:1:100;
sigma_z = 6.98516;% 6.98516; % cm
bunch_center = 3.7;
wakefields_direction = 'long';

search_type = 'max';
search_xi = 6.98516;
trans_range = [0.003 0.005];
useAvg = false;
dataformat = 'mat';

% Initialization
amplitudes = cell(length(datadirs),1);
positions = cell(length(datadirs),1);
propagations = cell(length(datadirs),1);

OWA = OsirisWakefieldAnalysis('dump_list',dump_list,...
    'wakefields_direction',wakefields_direction,...
    'useAvg',useAvg,...
    'search_type',search_type,...
    'plasmaden',plasmaden,'trans_range',trans_range,...
    'sigma_z',sigma_z,'bunch_center',bunch_center,...
    'dataformat',dataformat);

P = Plotty('plots_dir',plots_dir,'save_flag',saveplot,'plasmaden',plasmaden);


for d = 1:length(datadirs)
    OWA.datadir = datadirs{d};
    
    % OWA.xi_range = search_xi + 0.6*[OWA.plasma_wavelength,-OWA.plasma_wavelength];
    OWA.xi_range = [18 0];
    
    % find max field
    OWA.amplitude_vs_z();
    amplitudes{d} = OWA.denorm_Efield(OWA.amplitude_z);
    propagations{d} = OWA.propagation_z/100;
    positions{d} = OWA.pos_amplitude_z;
    
end

%% plot results
fig_amplitude = figure;
fig_amplitude.Units = 'normalized';
fig_amplitude.OuterPosition = [0 0.25 0.4 0.4];

for d = 1:length(datadirs)
    
    switch datadirs{d}
        case 'gp20'
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
    
    
    if isempty(colornumber)
        colornumber = d;
    end
    plot(propagations{d},amplitudes{d},'LineWidth',2,'color',cc(colornumber,:));
    
end

axis('tight');
xlim([0 10]);
ylim(y_lim);
xlabel('z (m)');
ylabel(y_label);
% title(title_g);
% legend(leg,'Location','best')
drawnow;
P.plot_name = ['st',plot_name];
P.fig_handle = fig_amplitude;
P.save_plot();


fig_position = figure;
fig_position.Units = 'normalized';
fig_position.OuterPosition = [0 0.25 0.6 0.55];

for d = 1:length(datadirs)
    plot(propagations{d},positions{d},'LineWidth',2);
    
    axis('tight');
    xlim([0 10]);
    
    % title(['',meanormax,' E field position along bunch']);
    xlabel('z (cm)');
    ylabel(y_label_position);
    %     legend(leg,'Location','best')
    drawnow;
    
    P.plot_name = ['pos_','st',plot_name];
    P.fig_handle = fig_position;
    P.save_plot();
end

