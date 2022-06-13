
clear;
close all;

load('color_red_to_blue.mat');
letters = {'a) $g = -2$\,\%/m','b) $g = 0$\,\%/m','c) $g = +2$\,\%/m'};
hletter = 0.15;
wletter = 0.015;
xlims_fft = [-1.6362,1.6262];
ylims_fft = [107,133]; % GHz
rdir = 'r15';

fig1 = openfig(['gradsim_convergence/fft/r/n133',rdir,'/fig/fvsrgm20',rdir,'m1350.fig']);
oax1 = fig1.Children(2);

fig2 = openfig(['gradsim_paper/fft/r/n133/fig/fvsrg0m1350.fig']);
oax2 = fig2.Children(2);

fig3 = openfig(['gradsim_convergence/fft/r/n133',rdir,'/fig/fvsrgp20',rdir,'m1350.fig']);
oax3 = fig3.Children(2);

P = Plotty('plasmaden',1.81e14);


figx = figure(10);
figx.Units = 'centimeters';
figx.OuterPosition = [1 1 8.6 8.6]*2;
% colororder(cc);

tt = tiledlayout(3,1);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

ax(1) = nexttile;
oax1.Children(4).MarkerFaceColor = ccrb(1,:);
copyobj(oax1.Children,ax(1));
ylim(ylims_fft);
xlim(xlims_fft);
ax(1).FontUnits = 'centimeters';
ax(1).FontSize = 0.6; %cm
ax(1).XTickLabel = [];
text(wletter,hletter,letters{1},'Units','normalized','FontUnits','centimeters',...
    'FontSize',0.5,'interpreter','latex')

ax(2) = nexttile;
oax2.Children(3).MarkerFaceColor = [0 0 0];
copyobj(oax2.Children,ax(2));
ylim(ylims_fft);
xlim(xlims_fft);
ax(2).FontUnits = 'centimeters';
ax(2).FontSize = 0.6; %cm
ax(2).XTickLabel = [];
text(wletter,hletter,letters{2},'Units','normalized','FontUnits','centimeters',...
    'FontSize',0.5,'interpreter','latex')

ax(3) = nexttile;
oax3.Children(4).MarkerFaceColor = ccrb(9,:);
copyobj(oax3.Children,ax(3));
ylim(ylims_fft);
xlim(xlims_fft);
ax(3).FontUnits = 'centimeters';
ax(3).FontSize = 0.6; %cm
text(wletter,hletter,letters{3},'Units','normalized','FontUnits','centimeters',...
    'FontSize',0.5,'interpreter','latex')

% ylabel(tt,['$f_{\mathrm{mod}}$ in ','transverse slice ','(GHz)'],'FontSize',15,'interpreter','latex');
ylabel(tt,['$f_{\mathrm{mod}}$ (GHz)'],'FontSize',15,'interpreter','latex');

xlabel(tt,'$x$ (mm)','FontSize',15,'interpreter','latex');
legend(ax(3),'Simulation','Experiment','location','southeast',...
    'FontSize',12,'box','off');
% legend('boxoff')
% handle


P.plots_dir = 'gradsim_convergence/fft/join';
P.plot_name = ['fvsrg',rdir];
P.fig_handle = figx;
P.save_plot();

