%________________________________________________________________________
% Script to compare the max long fields of the simulations with gap and no
% gap with density step
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 24/02/2020
%________________________________________________________________________

clear;
close all;

% run switch
run_it = 1;

plots_dir = ['gap/max_long_fields/',''];
plot_name = ['max_long_field'];

% load files
load('color_red_to_blue.mat'); % ccrb
% color selection
i_color = [1,9];

% cell plotting parameters
datadirs = {'R2gap0_2e14','R2gap100_2e14'};
leg = {'gap = 0 m','gap = 1 m'};
line_style = {':','--',};

% plotting parameters
fontsize_annotation = 12;
fontsize_label  = 14;

% study parameters
dataformat  = 'h5';
useAvg      = true;
dump_list   = 0:1:200;
plasmaden   = 2e14;
property    = 'fields';
trans_limit = 0.0536;
xi_range    = [7,5];


% Load the analysis class and initial charge
O = OsirisDenormalizer(...
    'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg',...
    'dump',0,'plasmaden',plasmaden,...
    'property',property,'wakefields_direction','long',...
    'trans_range',[0,trans_limit],'xi_range',xi_range);
P = Plotty('plots_dir',plots_dir,'plasmaden',plasmaden,...
    'plot_name',plot_name,'save_flag',1);

fieldvsz = zeros(length(datadirs),length(dump_list));
plot_z = zeros(length(datadirs),length(dump_list));

% begin loop
if run_it
    for d = 1:length(datadirs)
        O.datadir = datadirs{d};
        
        for n = 1:length(dump_list)
            O.dump = dump_list(n); O.getdata(); O.assign_fields();
            
            O.denorm_Efield(); O.denorm_distance();
            
            z_ind = O.z > O.dtime+O.simulation_window - O.xi_range(1) & ... %large
                O.z <= O.dtime+O.simulation_window - O.xi_range(2); % small
            r_ind = O.r >= O.trans_range(1) & ...
                O.r < O.trans_range(2);
            long_fieldinxi = O.longfield(r_ind,z_ind);
            fieldvsz(d,n) = max(long_fieldinxi,[],'all');
            %             field_temp_gauss = imgaussfilt(fields_for_gauss,67,'FilterDomain','spatial');
            %             fields_temp_mean = mean(fields_for_gauss,'all');
            %             fieldvsz(x,d,n) = mean(field_temp_gauss(field_temp_gauss > 0));
            plot_z(d,n) = O.propagation_distance/100;
            
            O.progress_dump('dump',n,length(dump_list))
            
        end % for dump
        
        O.progress_dump('directory',d,length(datadirs))
    end % for datadirs
else
    load('loading_files/gap_fieldamplitude_longmax.mat')
end

%% plotting

fig_cvsz = figure(1);
fig_cvsz.OuterPosition = [100 100 600 400];

hold on
for d = 1:length(datadirs)
    plot(plot_z(d,:),fieldvsz(d,:),...
        line_style{d},'LineWidth',2,'color',ccrb(i_color(d),:))
end % datadir
hold off

ylim([0,max(fieldvsz,[],'all')])
xlim([0,20])

legend(leg,'location','best','FontSize',fontsize_annotation)
xlabel('z (m)')
ylabel('max long field (MV/m)');

P.fig_handle = fig_cvsz;
P.save_plot();

save('loading_files/gap_fieldamplitude_longmax.mat','plot_z','fieldvsz');







