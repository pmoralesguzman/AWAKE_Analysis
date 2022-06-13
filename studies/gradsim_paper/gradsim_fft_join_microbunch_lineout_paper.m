
clear;
close all;

letters = {'a) $g = -2$\,\%/m','b) $g = +2$\,\%/m'};
legs = {'$z = 3.8$\,m','$z = 10$\,m','$z = 10$\,m (envelope)'};
hletter = 0.9;
% wletter = 0.85;
wletter = 0.01;

fig1 = openfig(['gradsim_paper/fft/lineout/fig/microbunchgm20.fig']);
oax1 = fig1.Children(1);

fig2 = openfig(['gradsim_paper/fft/lineout/fig/microbunchgp20.fig']);
oax2 = fig2.Children(2);

P = Plotty('plasmaden',1.81e14);

figx = figure(10);

figx.Units = 'centimeters';
figx.OuterPosition = [1,1,8.6*1.5,8.6*1.5*600/1400]*1.5*1.5;

tt = tiledlayout(2,1);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

ax(1) = nexttile;
ax(1).XDir = 'reverse';
copyobj(oax1.Children,ax(1));

ylim(1e11*[0 7]);
xlim([0.073,2]);
ax(1).FontUnits = 'centimeters';
ax(1).FontSize = 0.4; %cm
ax(1).YTickLabel = [];
text(wletter,hletter,letters{1},'Units','normalized',...
    'FontUnits','centimeters','FontSize',0.4,'interpreter','latex')

ax(2) = nexttile;
ax(2).XDir = 'reverse';
copyobj(oax2.Children,ax(2));
ylim(1e11*[0 8]);
xlim([0.073,14]);
ax(2).YTickLabel = [];
ax(2).FontUnits = 'centimeters';
ax(2).FontSize = 0.4;
text(wletter,hletter,letters{2},'Units','normalized',...
        'FontUnits','centimeters','FontSize',0.4,'interpreter','latex')

ylabel(tt,['density (arb. units)'],'FontSize',15,'interpreter','latex')
xlabel(tt,'$\xi$ (cm)','FontSize',15,'interpreter','latex')
% legend(ax(2),legs{:},'location','west','FontSize',10,'interpreter','latex');
legend(ax(2),legs{:},'position',[0.234192953707497 0.35756552020026 0.153634227235412 0.13659793814433],...
    'FontSize',10,'interpreter','latex');

legend('boxoff')


P.plots_dir = 'gradsim_paper/fft/join';
P.plot_name = 'microbunch';
P.fig_handle = figx;
P.save_plot();

