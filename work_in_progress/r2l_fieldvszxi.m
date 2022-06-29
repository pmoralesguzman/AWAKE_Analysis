%________________________________________________________________________
% 
% Script to calculate the fields within a region +/- lambda/2
% around a three initial xi0, following the zero-crossing of the fields 
% starting at that xi0%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 18/02/2020
%________________________________________________________________________

clear;
close all;

% run switch
run_it = 1;

plots_dir = ['eps2022'];
plot_name = ['fieldvszxi'];

% load files
load('color_red_to_blue.mat'); % ccrb
% color selection
% i_color = [3,5,7];
i_color = [1:9];

% cell plotting parameters
datadirs = {'r2l_302_c_e550_l','r2l_302_c_pi_e550_l'};

leg = {'e- bunch no shift','e- bunch shifted by $\lambda_{pe}/2$'};
line_style = {':','--','-.','-','-','-','-.','--',':'};
% line_style = {'--','-','-.'};

% plotting parameters
fontsize_annotation = 12;
fontsize_label  = 14;

% study parameters
dataformat  = 'mat';
useAvg      = false;
dump_list   = 6:1:100;
plasmaden   = 2e14;
property    = 'fields';
trans_limit = 0.02;
dephasing_xi= [7,3];


% Load the analysis class and initial charge
O = OsirisDenormalizer(...
    'datadir','r2l_302_c_e550_l','dataformat',dataformat,'useAvg',useAvg',...
    'dump',0,'plasmaden',plasmaden,...
    'property',property,'wakefields_direction','trans',...
    'trans_range',[0,trans_limit]);
P = Plotty('plots_dir',plots_dir,'plasmaden',plasmaden,...
    'plot_name',plot_name,'save_flag',1);

OWA = OsirisWakefieldAnalysis('datadir','r2l_302_c_e550_l','dataformat',dataformat,'useAvg',useAvg',...
    'dump',0,'plasmaden',plasmaden,...
    'property',property,'wakefields_direction','long',...
    'trans_range',[0,trans_limit]);

xi_ranges = dephasing_xi + [O.plasma_wavelength/2,-O.plasma_wavelength/2];
fieldvsz = zeros(length(dephasing_xi),length(datadirs),length(dump_list));
plot_z = zeros(length(dephasing_xi),length(datadirs),length(dump_list));
field_ini = ones(length(dephasing_xi),length(datadirs),1);
% begin loop

    for d = 1:length(datadirs)
        O.datadir = datadirs{d};
        
        for n = 1:length(dump_list)
            O.plasmaden = plasmaden; 
            O.dump = dump_list(n); O.getdata(); O.assign_fields();
            O.denorm_distance(); plasma_wavelength0 = O.plasma_wavelength;
            
            for xi = 1:length(dephasing_xi)

                O.xi_range = dephasing_xi(xi) ...
                    + [O.plasma_wavelength/2,-O.plasma_wavelength/2];

                O.denorm_Efield();

                z_ind = O.z > O.dtime+O.simulation_window - O.xi_range(1) & ... %large
                    O.z <= O.dtime+O.simulation_window - O.xi_range(2); % small
                
                r_ind = O.r >= O.trans_range(1) & ...
                    O.r < O.trans_range(2);
                trans_fieldinxi = O.transfield(r_ind,z_ind);

                %             fields_for_gauss = trans_fieldinxi;
                %             fields_for_gauss(~i_field) = 0;
                fields_for_mean = trans_fieldinxi;
                if n == 1
                    field_ini(xi,d,n) = max(abs(fields_for_mean),[],'all');
                end

                fieldvsz(xi,d,n) = max(abs(fields_for_mean),[],'all')/field_ini(xi,d,1);

                %             field_temp_gauss = imgaussfilt(fields_for_gauss,67,'FilterDomain','spatial');
                %             fields_temp_mean = mean(fields_for_gauss,'all');
                %             fieldvsz(x,d,n) = mean(field_temp_gauss(field_temp_gauss > 0));
                plot_z(xi,d,n) = O.propagation_distance/100;
                
            end % for xi
            O.progress_dump('dump',n,length(dump_list))
            
        end % for dump
        
        O.progress_dump('directory',d,length(datadirs))
    end % for datadirs


%% plotting

fig_cvsz = figure(2);
fig_cvsz.OuterPosition = [100 100 900 400];

tt = tiledlayout(1,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
for xi = 1:length(dephasing_xi)
    ax_fld(xi) = nexttile;
    ax_fld(xi).FontSize = fontsize_label;
    
    hold on
    for d = 1:length(datadirs)

        plot(squeeze(plot_z(xi,d,:)),squeeze(fieldvsz(xi,d,:)),...
            line_style{d},'LineWidth',2,'color',ccrb(i_color(d),:))
    end % datadir
    xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
    hold off
    
    switch dephasing_xi(xi)
        case 3
            position_word = '(front of thebunch)';
        case 7
            position_word = '(back of the bunch)';
        case 1
            position_word = '(front of the bunch)';
    end
    title(['$\xi_0 = ',num2str(dephasing_xi(xi)),'$ cm ',position_word],'Interpreter','latex');
    ylim([0,max(fieldvsz,[],'all')])
    xlim([0,10])
    
end % xi

legend(ax_fld(2),leg,'location','northeast','FontSize',fontsize_annotation,'Interpreter','latex')
xlabel(tt,'z (m)','Interpreter','latex')
ylabel(tt,{'$W_{r,\mathrm{max}}(z) / W_{r,\mathrm{max}}(z=50\,\mathrm{cm})$'},'Interpreter','latex');

P.fig_handle = fig_cvsz;
P.save_plot();








