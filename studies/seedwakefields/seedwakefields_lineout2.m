%________________________________________________________________________
% Calculate and plot the mean transverse fields, together with manually
% positioned tags.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 03/06/2021
%________________________________________________________________________

% answer = questdlg('Are you sure you want to run? It might take long time.');
% switch answer
%     case 'Yes'
%         % continue
%     otherwise
%         return;
% end

clear;
close all;

load('color_red_to_blue.mat'); % load plot colors
% ccrb = ccrb(5:end,:);

run_it = 0; % run the code without loading data from cache

datadirs = {'g0'};

% save file name
plots_dir = ['gradsim_paper/ch/wrvsz'];
plot_name = 'meandefocusing';
cache_dir = ['transmean'];
saveplot = 1;

% simulation parameters
plasmaden = 1.81e14;
dump_list = 1:1:2;
wakefields_direction = 'trans';

search_type = 'mean';
trans_range = [0 100];%[0.002 0.01];
useAvg = 0;
dataformat = 'mat';

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
    load('gradsim_fieldamplitude_transmean.mat');
end % if run it

%% plot results

fontsize_annotation = 9; % points (1 point = 1/72 inches = 0.0353 cm; 9 point = 0.388)
fontsize_label = 0.4; % cm

fig_amplitude = figure(5);
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
ylabel(['mean defocusing wakefields (MV/m)'],'interpreter','latex');
addtextlines = 0;
% if addtextlines
    % Create textarrow
   % Create textarrow
annotation(fig_amplitude,'textarrow',[0.490553602811951 0.510544815465729],...
    [0.32751744765703 0.383848454636092],'Color',[0.03125 0.31640625 0.609375],...
    'String','$g = +2$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.516322262235806 0.486005917771834],...
    [0.823395752068133 0.801461554660356],...
    'Color',[0.0390625 0.0390625 0.0390625],...
    'String','$g = 0$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.550727214837622 0.559544438036216],...
    [0.618605300950436 0.52778775359251],...
    'Color',[0.12890625 0.44140625 0.70703125],...
    'String','$g = +1.5$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.660477743516462 0.63016139905249],...
    [0.723633186954593 0.701698989546816],...
    'Color',[0.2578125 0.5703125 0.7734375],...
    'String','$g = +1$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.826573421898078 0.843269379718816],...
    [0.666675177835237 0.716525626489275],...
    'Color',[0.41796875 0.6796875 0.8359375],...
    'String','$g = +0.5$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.296901293604336 0.358852084465496],...
    [0.623158861056795 0.531434035533365],...
    'Color',[0.64453125 0.05859375 0.08203125],...
    'String','$g = -2$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.310115148799567 0.365475429994646],...
    [0.694288252096996 0.597578381708162],...
    'Color',[0.79296875 0.09375 0.11328125],...
    'String','$g = -1.5$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.336234835331237 0.386762076104523],...
    [0.75027234308177 0.654559481666018],...
    'Color',[0.93359375 0.23046875 0.171875],...
    'String','$g = -1$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.37387657953999 0.421328249135772],...
    [0.826386778987375 0.717712800921572],...
    'Color',[0.98046875 0.4140625 0.2890625],...
    'String','$g = -0.5$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);
% end % if add text lines
if addtextlines
% Create textarrow
annotation(fig_amplitude,'textarrow',[0.422459086993971 0.466214470284238],...
    [0.878074456885826 0.797856254187004],...
    'Color',[0.0390625 0.0390625 0.0390625],...
    'String','$g = 0$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.336234835331237 0.379615274188918],...
    [0.75027234308177 0.650789549239162],...
    'Color',[0.93359375 0.23046875 0.171875],...
    'String','$g = -1$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.296901293604336 0.353179730117715],...
    [0.623158861056795 0.521954253995598],...
    'Color',[0.64453125 0.05859375 0.08203125],...
    'String','$g = -2$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.310115148799567 0.358649153028998],...
    [0.694288252096996 0.588802756244617],...
    'Color',[0.79296875 0.09375 0.11328125],...
    'String','$g = -1.5$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.34041774332472 0.396023542922768],...
    [0.823380227772992 0.716422624174562],...
    'Color',[0.98046875 0.4140625 0.2890625],...
    'String','$g = -0.5$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.773413723801321 0.808964972724663],...
    [0.868351038376878 0.792994544932529],...
    'Color',[0.41796875 0.6796875 0.8359375],...
    'String','$g = +0.5$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.632120298593167 0.592922767728969],...
    [0.766255144032922 0.773547707914633],...
    'Color',[0.2578125 0.5703125 0.7734375],...
    'String','$g = +1$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.372322710307207 0.350445018662073],...
    [0.355440712029859 0.438089769355919],'Color',[0.03125 0.31640625 0.609375],...
    'String','$g = +2$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);

% Create textarrow
annotation(fig_amplitude,'textarrow',[0.770679012345679 0.736950904392765],...
    [0.681175232079625 0.538970236386257],...
    'Color',[0.12890625 0.44140625 0.70703125],...
    'String','$g = +1.5$\,\%/m',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',9);
end
% legend(leg,'Location','best')
P.plot_name = plot_name;
P.fig_handle = fig_amplitude;
P.save_plot();

save(['loading_files/gradsim_cache/gradsim_fieldamplitude_',cache_dir,'.mat'],...
    'propagations','amplitudes');



