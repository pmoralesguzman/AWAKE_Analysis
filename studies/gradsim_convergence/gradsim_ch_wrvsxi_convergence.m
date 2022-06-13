%________________________________________________________________________
% gradsim paper
% Script to calculate the mean defocusing fields within a region +/- lambda/2
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
% close all;

% run switch
run_it = 0;

plots_dir = ['gradsim_convergence/ch/wrvszxi'];
plot_name = ['wrvszxi'];

% load files
load('color_red_to_blue.mat'); % ccrb
load('loading_files/gradsim_cache/gradsim_dephasing_convergence.mat');
% color selection
% i_color = [3,5,7];
i_color = [1:9];

% cell plotting parameters
datadirs = {'gm20','gm15','gm10d2','gm5','g0','gp5','gp10d2','gp15','gp20'};
grads = [-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]/100;
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
line_style = {':','--','-.','-','-','-','-.','--',':'};
% line_style = {'--','-','-.'};

% plotting parameters
fontsize_annotation = 12;
fontsize_label  = 14;

% study parameters
dataformat  = 'mat';
useAvg      = false;
dump_list   = 0:1:100;
plasmaden   = 1.81e14;
property    = 'fields';
trans_limit = 0.0536;
dephasing_xi= [14;7;1];


% Load the analysis class and initial charge
O = OsirisDenormalizer(...
    'datadir','g0','dataformat',dataformat,'useAvg',useAvg',...
    'dump',0,'plasmaden',plasmaden,...
    'property',property,'wakefields_direction','trans',...
    'trans_range',[0,trans_limit]);
P = Plotty('plots_dir',plots_dir,'plasmaden',plasmaden,...
    'plot_name',plot_name,'save_flag',1);

xi_ranges = dephasing_xi + [O.plasma_wavelength/2,-O.plasma_wavelength/2];
fieldvsz = zeros(length(dephasing_xi),length(datadirs),length(dump_list));
plot_z = zeros(length(dephasing_xi),length(datadirs),length(dump_list));

% begin loop
if run_it
    for d = 1:length(datadirs)
        O.datadir = datadirs{d};
        switch O.datadir
            case {'gm10d2','gp10d2','gm20r17','gp20r17','gp20r15','gm20r15'}
                O.dataformat = 'h5';
                O.useAvg = 1;
%                 O.trans_range = [0,0.0477];
            case {'gp20r15','gm20r15'}
                O.dataformat = 'h5';
                O.useAvg = 0;
%                 O.trans_range = [0,0.0408];
            otherwise
                O.dataformat = 'mat';
                O.useAvg = 0;
%                 O.trans_range = [0,trans_limit];
        end
        
        for n = 1:length(dump_list)
            O.plasmaden = plasmaden; 
            O.dump = dump_list(n); O.getdata(); O.assign_fields();
            O.denorm_distance(); plasma_wavelength0 = O.plasma_wavelength;
            for xi = 1:length(dephasing_xi)

                if n < 4
                    O.xi_range = dephasing_xi(xi) ...
                        + [O.plasma_wavelength/2,-O.plasma_wavelength/2];
                else
                    O.plasmaden = plasmaden.*(1 + grads(d)*O.propagation_distance_m);
                    O.den2freq();
                    O.xi_range = dephasing_xi(xi) + dephasing_lines(xi,d,n-3)*plasma_wavelength0/2 ...
                        + [O.plasma_wavelength/2,-O.plasma_wavelength/2];
                    O.plasmaden = plasmaden;
                    O.den2freq();
                end
                
                O.denorm_Efield(); 
                
                z_ind = O.z > O.dtime+O.simulation_window - O.xi_range(1) & ... %large
                    O.z <= O.dtime+O.simulation_window - O.xi_range(2); % small
                
                r_ind = O.r >= O.trans_range(1) & ...
                    O.r < O.trans_range(2);
                trans_fieldinxi = O.transfield(r_ind,z_ind);
                i_field = trans_fieldinxi > 0;
                %             fields_for_gauss = trans_fieldinxi;
                %             fields_for_gauss(~i_field) = 0;
                fields_for_mean = trans_fieldinxi;
                fields_for_mean(~i_field) = 0;
                fieldvsz(xi,d,n) = mean(fields_for_mean,'all');
                %             field_temp_gauss = imgaussfilt(fields_for_gauss,67,'FilterDomain','spatial');
                %             fields_temp_mean = mean(fields_for_gauss,'all');
                %             fieldvsz(x,d,n) = mean(field_temp_gauss(field_temp_gauss > 0));
                plot_z(xi,d,n) = O.propagation_distance/100;
                
            end % for xi
            O.progress_dump('dump',n,length(dump_list))
            
        end % for dump
        
        O.progress_dump('directory',d,length(datadirs))
    end % for datadirs
else
    load('loading_files/gradsim_cache/gradsim_fieldamplitude_transmean_xi_convergence.mat')
end

%% plotting

fig_cvsz = figure(1);
fig_cvsz.OuterPosition = [100 100 1200 400];

tt = tiledlayout(1,3);
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
        case 14
            position_word = '(back of the bunch)';
        case 7
            position_word = '(middle of the bunch)';
        case 1
            position_word = '(front of the bunch)';
    end
    title(['\xi_0 = ',num2str(dephasing_xi(xi)),' cm ',position_word]);
    ylim([0,max(fieldvsz,[],'all')])
    xlim([0,10])
    
end % xi

legend(ax_fld(3),leg,'location','northeast','FontSize',fontsize_annotation)
xlabel(tt,'z (m)')
ylabel(tt,'mean defocusing field (MV/m)');

P.fig_handle = fig_cvsz;
P.save_plot();

save('loading_files/gradsim_cache/gradsim_fieldamplitude_transmean_xi_convergence.mat','plot_z','fieldvsz');







