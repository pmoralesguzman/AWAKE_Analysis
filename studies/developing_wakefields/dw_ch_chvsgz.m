%________________________________________________________________________
% Evolution of the total charge along propagation for each gradient
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 04/06/2030
%________________________________________________________________________

answer = questdlg('Are you sure?');

switch answer
    case 'Yes'
        % continue
    otherwise
        return;
end

% clear;
close all;

% load files
load('color_red_to_blue.mat');

run_it = 1;

% data directory
datadirs    = {'DWb','DWdc','DWdc2','DWdc3'};
% datadirs    = {'gm20'};


leg         = {'baseline proton bunch','in-phase','somewhat out-of-phase','very out-of-phase'};
dataformat  = 'h5';
useAvg      = 1;
initialdump = 0;
dump_list   = 0:2:100;

% save directory
plots_dir           = ['DW/ch/chvsgz/',''];
plot_name_suffix    = [''];
save_format         = {'png'};
plot_name           = 'chargeevolution';

% properties
plasma_density  = 2e14;
property        = 'density';


% switches
save_plot_flag      = true;


% now that we have the analytical bunch fractions, calculate the initial charge
initial_charge = 1;

% initialize charges
chargevsz = zeros(length(datadirs),length(dump_list));
plot_z = zeros(length(datadirs),length(dump_list));

% Load the analysis class and initial charge
O = OsirisDenormalizer(...
    'datadir','g0','dataformat',dataformat,'useAvg',useAvg',...
    'dump',initialdump,'plasmaden',plasma_density,...
    'property',property,'raw_dataset','q',...
    'species','density_feature','direction','r',...
    'trans_range',[0,0.02]);

P = Plotty('plots_dir',plots_dir,'plasmaden',plasma_density,...
    'plot_name',plot_name,'save_flag',save_plot_flag);

if run_it == 1
    
    % begin loop
    for d = 1:length(datadirs)

        % select study directory

        O.datadir = datadirs{d};
        if d == 1
            O.species = 'proton_beam';
        else
            O.species = 'density_feature';
        end
        
        for n = 1:length(dump_list)

            O.dump = dump_list(n); O.getdata();
            plot_z(d,n) = (O.denorm_distance(O.ntime)-O.denorm_distance(O.n_simulation_window)/2)/100;
            switch O.property
                case 'raw'
                    
                    O.raw_dataset = 'x'; O.getdata(); O.assign_raw();
                    O.r_raw = O.denorm_distance(O.nr_raw);
                    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
                    O.denorm_distance();
                    
                    ind_translim = O.r_raw < trans_limit;
                    if n == 1
                        if dump_list(n) < 2
                            initial_charge = sum(O.q_raw)/(bunch_long_fraction - bunch_back_fraction);
                        else
                            error('for raw the first dump must be 0 or 1 to calculate initial charge')
                        end % if dump list 2
                    end % if n 1
                    
                    chargevsz(d,n) = sum(O.q_raw(ind_translim));
                    O.propagation_distance = O.dtime - 21/2;

                case 'density'

                    O.trim_data(); O.denorm_density();
                    chargevsz(d,n) = O.cylindrical_integration(O.r,O.z,O.(O.species),'sum');

                    if n == 1
                        if dump_list(n) < 2
                            initial_charge(d) = chargevsz(d,n);
                        else
                            error('for raw the first dump must be 0 or 1 to calculate initial charge')
                        end % if dump list 2
                    end % if n 1
                    
            end %switch property
            
            O.progress_dump(['charge evo. ',O.datadir],n,length(dump_list));
        end % for dump list
        
        
    end % for datadirs
    charge_fraction = chargevsz./initial_charge';
    
    
else
    
    load('loading_files/DW_cache/DW_charge_evolution.mat');
    
end % if run it

%% plotting

fig_cvsz = figure(1);
fig_cvsz.Units = 'centimeters';
fig_cvsz.Position = [1,1,8.6,8.6*3/4]*1.5;
colororder(ccrb);
% line_style = {'-.','-.','-.','-.','-','--','--','--','--'};
line_style = {':','--','-.','-','-','-','-.','--',':'};
fontsize_annotation = 9; % points (1 point = 1/72 inches = 0.0353 cm; 9 point = 0.388)
fontsize_label = 0.4; % cm

for d = 1:length(datadirs)
    hold on
    plot(plot_z(d,:),charge_fraction(d,:),line_style{d},'LineWidth',2)
    hold off
end
hold on
xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
xline(10,'--','LineWidth',1,'color',[0 0.0 0.0]);
hold off
xlim([0,12])
ax = gca;
ax.FontUnits = 'centimeters';
ax.FontSize = fontsize_label;
xlabel('$z$ (m)','interpreter','latex')
ylabel('charge fraction','interpreter','latex');
legend(leg,'Location','southwest','Autoupdate','off',...
    'FontSize',fontsize_annotation,'Interpreter','latex')
legend('boxoff')

drawnow;

P.fig_handle = fig_cvsz;
P.save_plot();

if run_it
    save('loading_files/DW_cache/DW_charge_evolution.mat','plot_z','charge_fraction');
end
