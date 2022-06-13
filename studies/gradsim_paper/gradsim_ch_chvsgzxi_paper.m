%________________________________________________________________________
% gradsim paper
% Script to calculate the charge fraction within a region +/- lambda/2
% around a three initial xi0, following the zero-crossing of the fields 
% starting at that xi0
% 
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 18/02/2021
%________________________________________________________________________

answer = questdlg('Are you sure?');

switch answer
    case 'Yes'
        % continue
    otherwise
        return;
end

clear;
close all;

% run switch
run_it = 1;

plots_dir = ['gradsim_paper/ch/chvsgzxi/',''];
plot_name = ['chargeinxi'];

% load files
load('color_red_to_blue.mat'); % ccrb
load('loading_files/gradsim_cache/gradsim_dephasing.mat');
i_color = [1:9];

% cell plotting parameters
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
grads = [-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]/100;

leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
line_style = {':','--','-.','-','-','-','-.','--',':'};

% plotting parameters
fontsize_annotation = 12;
fontsize_label = 14;

% study parameters
dataformat  = 'mat';
useAvg      = 0;
dump_list   = 0:1:100;
plasmaden   = 1.81e14;
property    = 'density';
trans_limit = 0.0545;
dephasing_xi= [14;7;1];

% Load the analysis class and initial charge
O = OsirisDenormalizer(...
    'datadir','g0','dataformat',dataformat,'useAvg',useAvg',...
    'dump',0,'plasmaden',plasmaden,...
    'property',property,'raw_dataset','q',...
    'species','proton_beam','direction','r',...
    'trans_range',[0,trans_limit]);
P = Plotty('plots_dir',plots_dir,'plasmaden',plasmaden,...
    'plot_name',plot_name,'save_flag',true);

xi_ranges = dephasing_xi + [O.plasma_wavelength/2,-O.plasma_wavelength/2];
chargevsz = zeros(length(dephasing_xi),length(datadirs),length(dump_list));
plot_z = zeros(length(dephasing_xi),length(datadirs),length(dump_list));
plot_z_sigma = zeros(length(dump_list),1);
sigma_beam = zeros(length(dump_list),1);
 
if run_it
    % begin loop
    for d = 1:length(datadirs)  
        O.datadir = datadirs{d};
        
        
        for n = 1:length(dump_list)
            O.plasmaden = plasmaden; plasma_wavelength0 = O.plasma_wavelength;
            O.dump = dump_list(n); O.getdata(); O.assign_density(); O.denorm_distance(); O.denorm_density();
            plot_z_sigma(n) = (O.denorm_distance(O.ntime)-O.denorm_distance(O.n_simulation_window)/2)/100;
            sigma_beam(n) = 0.02*sqrt(1 + (plot_z_sigma(n))^2/4.7382^2);
            O.trans_range = [0,sigma_beam(n)];
            for xi = 1:length(dephasing_xi)

                if n < 4
                    O.xi_range = dephasing_xi(xi) ...
                    + [O.plasma_wavelength/2,-O.plasma_wavelength/2];
                else
                    O.plasmaden = plasmaden.*(1 + grads(d)*O.propagation_distance_m);
                    O.den2freq();
                    O.xi_range = dephasing_xi(xi) + dephasing_lines(xi,d,n-3)*plasma_wavelength0/2 ...
                    + [O.plasma_wavelength/2, -O.plasma_wavelength/2];
                    O.plasmaden = plasmaden;
                    O.den2freq();
                end

                z_ind = O.z > O.dtime + O.simulation_window - O.xi_range(1) & ... %large
                    O.z <= O.dtime + O.simulation_window - O.xi_range(2); % small
                
                r_ind = O.r >= O.trans_range(1) & O.r < O.trans_range(2);
                
                chargevsz(xi,d,n) = O.cylindrical_integration(O.r(r_ind),O.z(z_ind),O.(O.species)(r_ind,z_ind),'trapz');
                plot_z(xi,d,n) = O.propagation_distance/100;
                if n == 1
                    proton_n1 = O.(O.species); r_n1 = r_ind;
                end
                charge_norm_xi(xi,d,n) = O.cylindrical_integration(O.r(r_n1),O.z(z_ind),proton_n1(r_n1,z_ind),'trapz');
                
            end % for xi
            O.progress_dump('dump',n,length(dump_list))
            
        end % for dump
        
        O.progress_dump('xi',d,length(datadirs))
    end % for datadirs
else
    load('loading_files/gradsim_chargeinxi.mat');
end % if run it


%% plotting

fig_cvsz = figure(1);
fig_cvsz.OuterPosition = [100 100 1200 400];

tt = tiledlayout(1,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

for xi = 1:length(dephasing_xi)
    ax_ch(xi) = nexttile;
    ax_ch(xi).FontSize = fontsize_label;
%     charge_norm_xi = mean(squeeze(chargevsz(xi,:,1)));
    
    hold on
    for d = 1:length(datadirs)        
        plot(squeeze(plot_z(xi,d,:)),squeeze(chargevsz(xi,d,:))./squeeze(charge_norm_xi(xi,d,:)),...
            line_style{d},'LineWidth',2,'color',ccrb(i_color(d),:))
    end % datadir
    xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
    hold off
    
    switch dephasing_xi(xi)
        case 14
            position_word = '(back)';
        case 7
            position_word = '(middle)';
        case 1
            position_word = '(front)';
    end
    title(['\xi_0 = ',num2str(dephasing_xi(xi)),' cm ',position_word]);
%     ylim([0,max(chargevsz,[],'all')])
    ylim([0,1.5])
    xlim([0,10])
    
end % xi


legend(ax_ch(3),leg,'location','southwest','FontSize',11)
xlabel(tt,'z (m)')
ylabel(tt,'total charge (arb. units)');

P.fig_handle = fig_cvsz;
P.save_plot();

save('loading_files/gradsim_chargeinxi.mat','plot_z','chargevsz');








