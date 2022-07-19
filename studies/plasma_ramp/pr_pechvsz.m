%________________________________________________________________________
% Evolution of the total charge along propagation
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 19/07/2022
%#ok<*NBRAK2> 
%________________________________________________________________________

% answer = questdlg('Are you sure?');
% 
% switch answer
%     case 'Yes'
%         % continue
%     otherwise
%         return;
% end

% clear;
close all;

% load files
load('color_red_to_blue.mat');

run_it = 1;

% data directory
datadirs    = {'fdr_26'};

species     = 'electrons';

leg         = {'plasma ramp'};
dataformat  = 'mat';
useAvg      = false;
initialdump = 0;
dump_list   = 0:1:200;

% save directory
plots_dir           = ['plasma_ramp/pechvsgz/']; 
plot_name_suffix    = [''];
save_format         = {'png','eps'};
plot_name           = 'chargeevo';

% properties
plasma_density  = 2e14;
property        = 'density';
trans_range     = inf; 

% switches
save_plot_flag      = true;

% initialize charges
chargevsz = zeros(length(datadirs),length(dump_list));
plot_z = zeros(length(datadirs),length(dump_list));
charge_fraction = zeros(length(datadirs),length(dump_list));
exp_normalization = zeros(length(datadirs),length(dump_list));
sigma_beam = zeros(length(datadirs),length(dump_list));
bunch_trans_fraction = zeros(length(datadirs),length(dump_list));

% Load the analysis class and initial charge
O = OsirisDenormalizer(...
    'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg',...
    'dump',initialdump,'plasmaden',plasma_density,...
    'property',property,'raw_dataset','q',...
    'species',species,'direction','r',...
    'trans_range',[0,trans_range]);

P = Plotty('plots_dir',plots_dir,'plasmaden',plasma_density,...
    'plot_name',plot_name,'save_flag',save_plot_flag);

if run_it == 1
    
    % begin loop
    for d = 1:length(datadirs)
        
        % select study directory
        
        O.datadir = datadirs{d};
        
        for n = 1:length(dump_list)
            O.dump = dump_list(n); O.getdata();
            plot_z(d,n) = O.denorm_distance(O.ntime)/100;
            %sigma_beam(d,n) = 0.02*sqrt(1 + plot_z(d,n)^2/4.7382^2);
            %O.trans_range = [0,sigma_beam(d,n)];
            switch O.property
                case 'raw'
                    
                    O.raw_dataset = 'x'; O.getdata(); O.assign_raw();
                    O.r_raw = O.denorm_distance(O.nr_raw);
                    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
                    O.denorm_distance();
                    
                    ind_translim = O.r_raw < trans_limit;
                    if n == 1
                        if dump_list(n) < 2
                            initial_charge = sum(O.q_raw);
                        else
                            error('for raw the first dump must be 0 or 1 to calculate initial charge')
                        end % if dump list 2
                    end % if n 1
                    
                    chargevsz(d,n) = sum(O.q_raw(ind_translim));
                    O.propagation_distance = O.dtime - 21/2;
                    
                case 'density'

                    O.trim_data(); O.denorm_density();
                    chargevsz(d,n) = O.cylindrical_integration(O.r,O.z,O.(O.species),'sum_density');
                    if n == length(dump_list)
                        if dump_list(n) == length(dump_list)
                            initial_charge = chargevsz(d,n);
                        else
                            error('dump must be the last one to calculate end density')
                        end % if dump list 2
                    end % if n 1

            end %switch property
            
            O.progress_dump(['charge evo. ',O.datadir],n,length(dump_list));
        end % for dump list
        
    end % for datadirs
    
    charge_fraction = (chargevsz/initial_charge);
    
else
    
    load(['loading_files/plasma_ramp_charge_evolution',obj.datadir,'.mat']);
    
end % if run it

%% plotting

fig_cvsz = figure(1);
fig_cvsz.Units = 'centimeters';
fig_cvsz.Position = [1,1,8.6,8.6*3/4]*1.5;
colororder(ccrb);
line_style = {':','--','-.','-','-','-','-.','--',':'};
fontsize_annotation = 9; % points (1 point = 1/72 inches = 0.0353 cm; 9 point = 0.388)
fontsize_label = 0.4; % cm

for d = 1:length(datadirs)
    hold on
    plot(plot_z(d,:),charge_fraction(d,:),line_style{d},'LineWidth',2)
    hold off
end
% hold on
% xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
% xline(10,'--','LineWidth',1,'color',[0 0.0 0.0]);
% hold off
% xlim([0,12])
ax = gca;
ax.FontUnits = 'centimeters';
ax.FontSize = fontsize_label;
xlabel('$z$ (m)','interpreter','latex')
ylabel('charge fraction','interpreter','latex');
legend(leg,'Location','best','Autoupdate','off',...
    'FontSize',fontsize_annotation,'Interpreter','latex')
legend('boxoff')

if false
% Create textarrow
annotation(fig_cvsz,'textarrow',[0.62440476190476 0.573809523809524],...
    [0.808994708994716 0.666402116402117],...
    'Color',[0.12890625 0.44140625 0.70703125],...
    'String','$g > 0$\,\%/m',...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',fontsize_annotation);

% Create textarrow
annotation(fig_cvsz,'textarrow',[0.60477318403675 0.566487223657766],...
    [0.515877117427505 0.543831945640731],'String','$g = 0$\,\%/m',...
    'Interpreter','latex',... 
    'HeadStyle','none',...
    'FontSize',fontsize_annotation);

% Create textarrow
annotation(fig_cvsz,'textarrow',[0.473691604117961 0.481628112054468],...
    [0.310599100392382 0.450281640074922],...
    'Color',[0.79296875 0.09375 0.11328125],...
    'String','$g < 0$\,\%/m',...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',fontsize_annotation);
end


drawnow;

P.fig_handle = fig_cvsz;
P.save_plot();

if run_it
    save(['loading_files/plasma_ramp_charge_evolution',O.datadir,'.mat'],'plot_z','charge_fraction');
end
