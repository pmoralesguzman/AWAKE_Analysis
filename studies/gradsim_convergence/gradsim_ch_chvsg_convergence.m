%________________________________________________________________________
% gradsim paper
% Comparison between total charge measured experimentally and in
% simulations
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

clear;
% close all;
% load experimental data from Tatiana

navyblue = [30,144,255]/256;
crimsom = [220,20,60]/256;

charge_exp    = load('totalcharge_1sigma.txt');
errorbars_exp = load('errorbars_totalcharge_1sigma.txt');

% data directory
% datadirs    = {'gm20d2','gm15d2','gm10d2','gm5d2','g0d2','gp5d2','gp10d2','gp15d2','gp20d2'};
% datadirs    = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs = {'gm20r15','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20r15'};
grads_exp   = [-19.4,-9.3,-5.16,0.3,4.3,8.7,13,20]/1;
grads_sim   = [-20,-15,-10,-5,0,5,10,15,20];
dataformat  = 'h5';
useAvg      = 0;
initialdump = 0;
dump        = 119;

% save directory
plots_dir           = ['gradsim_convergence/ch/chvsg/',''];
plot_name_suffix    = ['15'];
save_format         = {'png','eps','fig'};


% properties
plasma_density      = 1.81e14;

% plasma parameters
seeding_position    = 3.8074; % (127) ps ps*1e-12*O.c_cm
sigma_z             = 6.98516; % cm
sigma_exp           = 0.0545; % cm 0.0536, beta = 0.02*srqt(1+12^2/4.7382^2) = 0.0545
exp_upperlimit      = sigma_exp;

% analysis
limitr              = 1; %

% switches
save_plot_flag      = 1;

% plot format
fontsize = 0.4; %cm

% calculated variables
trans_limit = sigma_exp*limitr; % trans. limit in cm
exp_steps   = size(charge_exp,1);
exp_r       = linspace(sigma_exp/(exp_steps+1),3*exp_upperlimit,exp_steps);

% analytical charge fraction
bunch_trans_fraction= 1-exp(-(trans_limit/sigma_exp).^2/2);
bunch_long_fraction = normcdf(seeding_position,0,sigma_z);
bunch_back_fraction = normcdf(seeding_position - 3*sigma_z,0,sigma_z);
bunch_fraction_outside_simulation_window ...
    = bunch_trans_fraction*(1 - bunch_long_fraction + bunch_back_fraction);

% Load the analysis class and initial charge
O = OsirisDenormalizer('datadir','g0','dataformat',dataformat,'useAvg',useAvg',...
    'dump',initialdump,'plasmaden',plasma_density,...
    'property','density','raw_dataset','q',...
    'species','proton_beam','direction','r',...
    'trans_range',[0,trans_limit]);

plot_name = ['bunchfrac',num2str(dump),'s',num2str(limitr),'den',plot_name_suffix];
P = Plotty('plots_dir',plots_dir,'plasmaden',plasma_density,...
    'plot_name',plot_name,'save_flag',save_plot_flag);

% now that we have the analytical bunch fractions, calculate the initial charge
initial_charge      = 3e11;

% the charge within 1 sigma of the unmodulated proton bunch is taken as 1.0
% in the experiment
sigma_normalization = 1.0;
exp_normalization   = 1 - exp(-(sigma_normalization)^2/2);

% initialize charges
charge_translim = zeros(length(trans_limit),length(grads_sim));
charge_fraction = zeros(length(trans_limit),length(grads_sim));


%% begin loop

O.dump = dump;

for d = 1:length(datadirs)
    
    % select study directory
    O.datadir = datadirs{d};
    switch O.datadir
        case {'gm20r17','gp20r17','gp20r15','gm20r15'}
            O.dataformat = 'h5';
            O.useAvg = 1;
            O.trans_range = [0,0.0477];
        case {'gp20r15','gm20r15'}
            O.dataformat = 'h5';
            O.useAvg = 1;
            O.trans_range = [0,0.0408];
        otherwise
            O.dataformat = 'mat';
            O.useAvg = 0;
            O.trans_range = [0,trans_limit];
    end
    
    O.getdata(); O.assign_density(); O.denorm_density(); O.denorm_distance();
    for r = 1:length(trans_limit)
        ind_r = O.r < trans_limit(r);
        charge_translim(r,d) = O.cylindrical_integration(O.r(ind_r),O.z,O.proton_beam(ind_r,:),'sum');
        charge_fraction(r,d) = charge_translim(r,d)/initial_charge ...
            + bunch_fraction_outside_simulation_window(r);
        O.progress_dump('charge fraction r',r + (d-1)*length(trans_limit),length(trans_limit)*length(datadirs));
    end
    
end
sim_charge = charge_fraction/exp_normalization;

%% plot bunch population of selected bunch for the gradients
marker_size = 10;

fig_cvsg = figure;
fig_cvsg.Units = 'centimeters';
fig_cvsg.Position = [1,1,8.6,8.6*3/4]*1.5;
[~,indr_exp] = min(abs(trans_limit-exp_r));
plot_exp = errorbar(grads_exp/10,charge_exp(indr_exp,:),errorbars_exp(indr_exp,:),...
    's','MarkerSize',marker_size,'Color',navyblue,'MarkerFaceColor',navyblue); %    DATA
hold on
plot_sim = plot(grads_sim/10,sim_charge,'o','MarkerSize',marker_size,...
    'Color',crimsom,'MarkerFaceColor',crimsom); % SIM
hold off
legend([plot_sim(1) plot_exp(1)],'Simulation','Experiment','location','northwest',...
    'FontSize',12);
legend('boxoff')
% 'Experiment \pm std. dev.'
ax = gca;
ax.FontUnits = 'centimeters';
ax.FontSize = fontsize;

% grid on
xlabel('$g$ (\%/m)','interpreter','latex')
ylabel('charge fraction','interpreter','latex')
xlim([-2.2,2.2])
%     ylim([0 0.85])

drawnow;
P.fig_handle = fig_cvsg;
P.save_plot();


