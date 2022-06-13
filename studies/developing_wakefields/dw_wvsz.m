%________________________________________________________________________
% Calculate and plot the mean fields,
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 25/10/2021
%________________________________________________________________________

answer = questdlg('Are you sure you want to run? It might take long time.');
switch answer
    case 'Yes'
        % continue
    otherwise
        return;
end

clear;
% close all;

load('color_red_to_blue.mat'); % load plot colors
% ccrb = ccrb(5:end,:);

run_it = 1; % run the code without loading data from cache

datadirs = {'DWb','DWdc','DWdc2','DWdc3'};
datadirs = {'DWb2','DWdc4'};
leg = {'baseline proton bunch','density feature','somewhat out-of-phase','very out-of-phase'};

% save file name
plots_dir = ['DW/wzvsz'];
plot_name = 'longmean';
cache_dir = ['longmean'];
saveplot = 1;

% simulation parameters
plasmaden = 1.81e14;
dump_list = 0:3:97;
wakefields_direction = 'long';

search_type = 'mean';
trans_range = [0.002 100];%[0.002 0.01];
useAvg = 1;
dataformat = 'h5';

% Initialization
amplitudes = cell(length(datadirs),1);
positions = cell(length(datadirs),1);
propagations = cell(length(datadirs),1);

P = Plotty('plots_dir',plots_dir,'save_flag',saveplot,'plasmaden',plasmaden);

if run_it == 1
    
    OWA = OsirisWakefieldAnalysis('dump_list',dump_list,...
        'wakefields_direction',wakefields_direction,...
        'useAvg',useAvg,...
        'search_type',search_type,...
        'plasmaden',plasmaden,'trans_range',trans_range,...
        'dataformat',dataformat);
    
    for d = 1:length(datadirs)
        datadir = datadirs{d};
        OWA.datadir = datadirs{d};

        % find max field
        OWA.amplitude_vs_z();
        amplitudes{d} = OWA.denorm_Efield(OWA.amplitude_z);
        propagations{d} = OWA.propagation_z/100; % propagation in meters
        positions{d} = OWA.pos_amplitude_z; 
        
        OWA.progress_dump('directory',d,length(datadirs));
        
    end % for datadirs
    
else
    load('DW_fieldamplitude_transmean.mat');
end % if run it

%% plot results

fontsize_annotation = 9; % points (1 point = 1/72 inches = 0.0353 cm; 9 point = 0.388)
fontsize_label = 0.4; % cm

fig_amplitude = figure(1);
colororder(ccrb);
line_style = {':','--','-.','-','-','-','-.','--',':'};

hold on
for d = 1:length(datadirs)
    plot(propagations{d},amplitudes{d},line_style{d},'LineWidth',2);
end
hold off
ax = gca;
ax.FontUnits = 'centimeters';
ax.FontSize = fontsize_label;
% fig_amplitude.Units = 'normalized';
% fig_amplitude.OuterPosition = [0 0.25 0.6 0.55];
fig_amplitude.Units = 'centimeters';
fig_amplitude.Position = [1,1,8.6,8.6*3/4]*1.5;
% axis('tight');
xlim([0 10]);
xlabel('$z$ (m)','interpreter','latex');
ylabel(['mean long. wakefields (MV/m)'],'interpreter','latex');

legend(leg,'Location','best')
P.plot_name = plot_name;
P.fig_handle = fig_amplitude;
P.save_plot();

if ~isfolder('loading_files/DW_cache/')
    mkdir('loading_files/DW_cache/')
end
save(['loading_files/DW_cache/DW_fieldamplitude_',cache_dir,'.mat'],...
    'propagations','amplitudes');



