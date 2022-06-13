
clear;
close all;

load('color_red_to_blue.mat');

fig1 = openfig(['gradsim_paper/fft/z/gm20/fig/fvszwaterfallps','n','18','.fig']);
oax1 = fig1.Children(2);

fig2 = openfig(['gradsim_paper/fft/z/gp20/fig/fvszwaterfallps','n','18','.fig']);
oax2 = fig2.Children(2);



figx = figure(10);
figx.Units = 'centimeters';
figx.OuterPosition = [1,1,8.6,8.6*2*3/4]*1.5;
% colororder(cc);

tt = tiledlayout(2,1);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

ax(1) = nexttile;
ax(1).FontUnits = 'centimeters';
ax(1).FontSize = 0.4;
oax1.Children(1).Color = ccrb(1,:);
oax1.Children(2).Color = [255,185,185]/256; %cc(3,:)/0.9336;
copyobj(oax1.Children,ax(1));


legend('maximum amplitude                                                  \','$f_{\mathrm{pe}}(z)$','interpreter','latex',...
    'Fontsize',12,'location','southwest','box','off')

colormap(flipud(gray));
xlim([0 13.5])
ylim([106 132])
xticks([]);

t1 = text(0.03,0.951,'a)','Units','normalized','FontUnits','centimeters','FontSize',0.4);

ax(2) = nexttile;
ax(2).FontUnits = 'centimeters';
ax(2).FontSize = 0.4;
oax2.Children(1).Color = ccrb(9,:);
oax2.Children(2).Color = [172, 217, 230]/256; %cc(6,:)/0.9;
copyobj(oax2.Children,ax(2));

legend('maximum amplitude                                                  \ ','$f_{\mathrm{pe}}(z)$','interpreter','latex',...
    'Fontsize',12,'location','southwest','box','off')

colormap(flipud(gray));
xlim([0 13.5])
ylim([106 132])
xlabel(tt,'z (m)','fontsize',12,'interpreter','latex')
ylabel(tt,'$f_{\mathrm{mod}}$ (GHz)','fontsize',12,'interpreter','latex')

% yticks([]);
t2 = text(0.03,0.951,'b)','Units','normalized','FontUnits','centimeters','FontSize',0.4);

P = Plotty('plasmaden',1.81e14);
P.plots_dir = 'gradsim_paper/fft/join/';
P.plot_name = 'waterfall';
P.fig_handle = figx;
P.save_plot();

