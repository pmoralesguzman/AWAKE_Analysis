
clear;
close all;

load('color_red_to_blue.mat');
cc = ccrb;
leg = {'$-2$\,\%/m','$-1.5$\,\%/m','$-1$\,\%/m','$-0.5$\,\%/m',...
    '\ $0$\,\%/m','$+0.5$\,\%/m','$+1$\,\%/m','$+1.5$\,\%/m','$+2$\,\%/m'};

fig1 = openfig(['gradsim_paper/ene/wz_ene/fig/wz_int.fig']);
oax1 = fig1.Children(2);

fig2 = openfig(['gradsim_paper/ene/wz_ene/fig/field_positive.fig']);
oax2 = fig2.Children(1);



figx = figure(10);
% figx.Units = 'pixels';
% figx.OuterPosition = [100 100 800 400];
figx.Units = 'centimeters';
figx.OuterPosition = [1,1,8.6,8.6*2*3/4]*1.5;
colororder(cc);

tt = tiledlayout(2,1);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

ax(1) = nexttile;
ax(1).FontUnits = 'centimeters';
ax(1).FontSize = 0.4; % cm

copyobj(oax1.Children,ax(1));
xlim([0 10])
ylim([-50 2000])
ylabel('energy gain (MeV)')
xlabel('distance from plasma end (m)')
ll = legend(ax(1),leg,'Location','northwest','AutoUpdate','off','interpreter','latex');
legend('boxoff')
text(0.9265,0.951,'a)','FontSize',14,'Units','normalized')

ax(2) = nexttile;
ax(2).FontUnits = 'centimeters';
ax(2).FontSize = 0.4; % cm
colororder(ax(2),cc(5:end,:));

copyobj(oax2.Children,ax(2));
% xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
xlim([0 10.0])
ylabel('maximum E_z (MV/m)')
xlabel('z (m)')
text(0.9265,0.951,'b)','FontSize',14,'Units','normalized')

P = Plotty();
P.plots_dir = 'gradsim_paper/ene/join/';
P.plot_name = 'energygainamplitude';
P.fig_handle = figx;
P.save_plot();

