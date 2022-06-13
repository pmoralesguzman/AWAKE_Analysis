
clear;
close all;

% file location variables for insets
datadirs = {'gm20','g0','gp20'};
datadirs_sim = {'gm20','g0','gp20'};

% saving data
save_flag = true;
save_format = {'pdf','eps','png'};

% choose limits (in cm for transverse, in ps for longitudinal)
xlims_fft = [-1.6362,3.1];
ylims_fft = [107,133]; % GHz
trans_range = [0 0.2];
xi_range1 = [50,90]; %ps
xi_range2 = [50,90]; %ps
xi_range3 = [90,130]; %ps

% figure number
fig_number = 10;

EDA = ExperimentalDataAnalyser('plasmaden',1.81e14);


%--------- LOAD FFT VS R FIGS ---------------------------------------------

load('color_red_to_blue.mat');
letters = {'a) $g = -2$\,\%/m','b) $g = 0$\,\%/m','c) $g = +2$\,\%/m'};
hletter = 0.85;
hletter2 = 0.15;
wletter = 0.015;

fig1 = openfig(['gradsim_paper/fft/r/n133/fig/fvsrgm20m1350.fig']);
oax1 = fig1.Children(2);

fig2 = openfig(['gradsim_paper/fft/r/n133/fig/fvsrg0m1350.fig']);
oax2 = fig2.Children(2);

fig3 = openfig(['gradsim_paper/fft/r/n133/fig/fvsrgp20m1350.fig']);
oax3 = fig3.Children(2);

%--------- PLOT FFT VS R FIGS ---------------------------------------------

P = Plotty('plasmaden',1.81e14,'save_format',save_format);

figx = figure(fig_number);
figx.Units = 'centimeters';
figx.OuterPosition = [1 1 8.6 8.6]*2;

tt = tiledlayout(6,3);
% tt = tiledlayout(3,1);
tt.TileSpacing = 'loose';
tt.Padding = 'compact';

ax(1) = nexttile(1,[2 3]);
% ax(1) = nexttile;
oax1.Children(4).MarkerFaceColor = ccrb(1,:);
copyobj(oax1.Children,ax(1));
ylim(ylims_fft);
xlim(xlims_fft);
ax(1).FontUnits = 'centimeters';
ax(1).FontSize = 0.6; %cm
ax(1).XTickLabel = [];
text(wletter,hletter,letters{1},'Units','normalized','FontUnits','centimeters',...
    'FontSize',0.5,'interpreter','latex')

% -- INSET g = -2 --
d = 1; % d = 1, g = -2

EDA.datadir = datadirs{d};
EDA.loadSCdata(); % P.SC
z_plot = EDA.SCI_z; %ps
r_plot = EDA.SCI_r*10; %mm
% density_plot = EDA.SCI;

%first manual trimming
ind_r = (r_plot > -trans_range(2)*10) & (r_plot < trans_range(2)*10);
ind_xi = (z_plot > xi_range1(1)) & (z_plot < xi_range1(2));
r_plot = r_plot(ind_r);
density_plot = EDA.SCI(ind_r,:);

ind_opaque = density_plot;
ind_opaque(density_plot > 4*std(density_plot(:,ind_xi),[],'all')) = 4*std(density_plot(:,ind_xi),[],'all');
ind_opaque = ind_opaque/max(ind_opaque,[],'all');

ax_density_exp(d) = axes(tt);
ax_density_exp(d).Layout.Tile = 6;

% set limits to the axis EXP
ax_density_exp(d).XLim = xi_range1;  % d = 1, g = -2
ax_density_exp(d).YLim = [min(r_plot),max(r_plot)]; % d = 1, g = -2
ax_density_exp(d).XDir = 'reverse'; % d = 1, g = -2
ax_density_exp(d).Box = 'on'; % d = 1, g = -2
ax_density_exp(d).LineWidth = 1; % d = 1, g = -2
ax_density_exp(d).TickLength = [0.03,0.025]; % d = 1, g = -2
ax_density_exp(d).XAxisLocation = 'top'; % d = 1, g = -2
ax_density_exp(d).XTickLabel = []; % d = 1, g = -2
ax_density_exp(d).YTickLabel = []; % d = 1, g = -2

imagesc(ax_density_exp(d),'XData',z_plot,'YData',r_plot,'CData',double(density_plot>0),'alphadata',ind_opaque);
grad = [1 1 1; 0 0 0];
colormap(ax_density_exp(d),grad);

% -- SIM g = -2

timetemp = load([datadirs_sim{d},'/n100_m1350_slit.mat']);
simtimeprofile = timetemp.densitymatrix;
z_plot = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2))-2.6352;
r_plot = linspace(0,timetemp.ylims(2),size(simtimeprofile,1))*10;

%first manual trimming
ind_r_sim = (r_plot > trans_range(1)*10) & (r_plot < trans_range(2)*10);
ind_xi = (z_plot > xi_range1(1)) & (z_plot < xi_range1(2));
density_plot_half = simtimeprofile(ind_r_sim,ind_xi);
r_plot = r_plot(ind_r_sim);
% z_plot = z_plot(ind_xi);

density_plot = [flip(simtimeprofile);simtimeprofile];
meanstd_density = 3*std(density_plot_half,[],'all');
max_opaqueness = 1;
ind_opaque = max_opaqueness*density_plot;
ind_opaque(density_plot > meanstd_density) = max_opaqueness*meanstd_density;
if d == 1
    max_indopaque_exp = max(ind_opaque,[],'all');
end

ind_opaque = ind_opaque/max_indopaque_exp;

ax_density_sim(d) = axes(tt);
ax_density_sim(d).Layout.Tile = 3;


% set limits to the axis
ax_density_sim(d).XLim = xi_range1; 
ax_density_sim(d).YLim = [-max(r_plot),max(r_plot)];
ax_density_sim(d).XDir = 'reverse';
ax_density_sim(d).Box = 'on';
ax_density_sim(d).XTickLabel = [];
ax_density_sim(d).YTickLabel = [];
ax_density_sim(d).LineWidth = 1;

imagesc(ax_density_sim(d),'XData',z_plot,'YData',[-max(r_plot),max(r_plot)],'CData',double(density_plot>0),'alphadata',ind_opaque);
colormap(ax_density_sim(d),grad);

% -------------------------- next g 0 ----


ax(2) = nexttile(7,[2 3]);
% ax(2) = nexttile;
oax2.Children(3).MarkerFaceColor = [0 0 0];
copyobj(oax2.Children,ax(2));
ylim(ylims_fft);
xlim(xlims_fft);
ax(2).FontUnits = 'centimeters';
ax(2).FontSize = 0.6; %cm
ax(2).XTickLabel = [];
text(wletter,hletter,letters{2},'Units','normalized','FontUnits','centimeters',...
    'FontSize',0.5,'interpreter','latex')



% -- INSET g = 0 --
d = 2;

EDA.datadir = datadirs{d};
EDA.loadSCdata(); % P.SC
z_plot = EDA.SCI_z; %ps
r_plot = EDA.SCI_r*10; %mm
ind_r = (r_plot > -trans_range(2)*10) & (r_plot < trans_range(2)*10);

ind_xi = (z_plot > xi_range2(1)) & (z_plot < xi_range2(2));
r_plot = r_plot(ind_r);
density_plot = EDA.SCI(ind_r,:);

ind_opaque = density_plot;
ind_opaque(density_plot > 4*std(density_plot(:,ind_xi),[],'all')) = 4*std(density_plot(:,ind_xi),[],'all');
ind_opaque = ind_opaque/max(ind_opaque,[],'all');

ax_density_exp(d) = axes(tt);
ax_density_exp(d).Layout.Tile = 12;

% set limits to the axis
ax_density_exp(d).XLim = xi_range2;
ax_density_exp(d).YLim = [min(r_plot),max(r_plot)];
ax_density_exp(d).XDir = 'reverse';
ax_density_exp(d).Box = 'on';
ax_density_exp(d).LineWidth = 1;
ax_density_exp(d).TickLength = [0.03,0.025];
ax_density_exp(d).XAxisLocation = 'top';
ax_density_exp(d).XTickLabel = [];
ax_density_exp(d).YTickLabel = [];

imagesc(ax_density_exp(d),'XData',z_plot,'YData',r_plot,'CData',double(density_plot>0),'alphadata',ind_opaque);
colormap(ax_density_exp(d),grad);

% --

timetemp = load([datadirs_sim{d},'/n100_m1350_slit.mat']);
simtimeprofile = timetemp.densitymatrix;
z_plot = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2))-2.6352*0-1;
r_plot = linspace(0,timetemp.ylims(2),size(simtimeprofile,1))*10;

%first manual trimming
ind_r_sim = (r_plot > trans_range(1)*10) & (r_plot < trans_range(2)*10);
ind_xi = (z_plot > xi_range2(1)) & (z_plot < xi_range2(2));
density_plot_half = simtimeprofile(ind_r_sim,ind_xi);
r_plot = r_plot(ind_r_sim);
% z_plot = z_plot(ind_xi);

density_plot = [flip(simtimeprofile);simtimeprofile];
meanstd_density = 3*std(density_plot_half,[],'all');
max_opaqueness = 1;
ind_opaque = max_opaqueness*density_plot;
ind_opaque(density_plot > meanstd_density) = max_opaqueness*meanstd_density;
if d == 1
    max_indopaque_exp = max(ind_opaque,[],'all');
end

ind_opaque = ind_opaque/max_indopaque_exp;

ax_density_sim(d) = axes(tt);
ax_density_sim(d).Layout.Tile = 9;

% set limits to the axis
ax_density_sim(d).XLim = xi_range2; 
ax_density_sim(d).YLim = [-max(r_plot),max(r_plot)];
ax_density_sim(d).XDir = 'reverse';
ax_density_sim(d).Box = 'on';
ax_density_sim(d).XTickLabel = [];
ax_density_sim(d).YTickLabel = [];
ax_density_sim(d).LineWidth = 1;

imagesc(ax_density_sim(d),'XData',z_plot,'YData',[-max(r_plot),max(r_plot)],'CData',double(density_plot>0),'alphadata',ind_opaque);
colormap(ax_density_sim(d),grad);


% -------------------------- next g ----


ax(3) = nexttile(13,[2 3]);
% ax(3) = nexttile;
oax3.Children(4).MarkerFaceColor = ccrb(9,:);
copyobj(oax3.Children,ax(3));
ylim(ylims_fft);
xlim(xlims_fft);
ax(3).FontUnits = 'centimeters';
ax(3).FontSize = 0.6; %cm
ax(3).XTickLabel = {'-1','0','1'}; %cm
ax(3).XTick = [-1,0,1]; %cm
text(wletter,hletter2,letters{3},'Units','normalized','FontUnits','centimeters',...
    'FontSize',0.5,'interpreter','latex')

% -- INSET g = +2 --
d = 3;

EDA.datadir = datadirs{d};
EDA.loadSCdata(); % P.SC
z_plot = EDA.SCI_z; %ps
r_plot = EDA.SCI_r*10; %mm

%first manual trimming
ind_r = (r_plot > -trans_range(2)*10) & (r_plot < trans_range(2)*10);
ind_xi = (z_plot > xi_range3(1)) & (z_plot < xi_range3(2));
r_plot = r_plot(ind_r);
density_plot = EDA.SCI(ind_r,:);

ind_opaque = density_plot;
gp20_opq = 2;
ind_opaque(density_plot > gp20_opq*std(density_plot(:,ind_xi),[],'all')) = gp20_opq*std(density_plot(:,ind_xi),[],'all');
ind_opaque = ind_opaque/max(ind_opaque,[],'all');

ax_density_exp(d) = axes(tt);
ax_density_exp(d).Layout.Tile = 18;

% set limits to the axis
ax_density_exp(d).XLim = xi_range3;
ax_density_exp(d).YLim = [min(r_plot),max(r_plot)];
ax_density_exp(d).XDir = 'reverse';
ax_density_exp(d).Box = 'on';
ax_density_exp(d).LineWidth = 1;
ax_density_exp(d).TickLength = [0.03,0.025];
ax_density_exp(d).XAxisLocation = 'top';
% xlabel('$t$ (ps)','FontSize',9,'Interpreter','latex');
% ylabel('$x$ (mm)','FontSize',9,'Interpreter','latex')
ax_density_exp(d).XTickLabel = [];
ax_density_exp(d).YTickLabel = [];

imagesc(ax_density_exp(d),'XData',z_plot,'YData',r_plot,'CData',double(density_plot>0),'alphadata',ind_opaque);
colormap(ax_density_exp(d),grad);

% --

timetemp = load([datadirs_sim{d},'/n100_m1350_slit.mat']);
simtimeprofile = timetemp.densitymatrix;
z_plot = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2))-2.6352;
r_plot = linspace(0,timetemp.ylims(2),size(simtimeprofile,1))*10;

%first manual trimming
ind_r_sim = (r_plot > trans_range(1)*10) & (r_plot < trans_range(2)*10);
ind_xi = (z_plot > xi_range3(1)) & (z_plot < xi_range3(2));
density_plot_half = simtimeprofile(ind_r_sim,ind_xi);
r_plot = r_plot(ind_r_sim);
% z_plot = z_plot(ind_xi);

density_plot = [flip(simtimeprofile);simtimeprofile];
meanstd_density = 3*std(density_plot_half,[],'all');
max_opaqueness = 1;
ind_opaque = max_opaqueness*density_plot;
ind_opaque(density_plot > meanstd_density) = max_opaqueness*meanstd_density;
if d == 1
    max_indopaque_exp = max(ind_opaque,[],'all');
end

ind_opaque = ind_opaque/max_indopaque_exp;

ax_density_sim(d) = axes(tt);
ax_density_sim(d).Layout.Tile = 15;


% set limits to the axis
ax_density_sim(d).XLim = xi_range3;
ax_density_sim(d).YLim = [-max(r_plot),max(r_plot)];
ax_density_sim(d).XDir = 'reverse';
ax_density_sim(d).Box = 'on';
ax_density_sim(d).LineWidth = 1;

imagesc(ax_density_sim(d),'XData',z_plot,'YData',[-max(r_plot),max(r_plot)],'CData',double(density_plot>0),'alphadata',ind_opaque);
colormap(ax_density_sim(d),grad);


ax_density_sim(d).XTickLabel = [];
ax_density_sim(d).YTickLabel = [];

ylabel(tt,['$f_{\mathrm{mod}}$ (GHz)'],'FontSize',15,'interpreter','latex');

xlabel(tt,'$x$ (mm)','FontSize',15,'interpreter','latex');

legend(ax(3),'Simulation','Experiment','Position',[0.378760355573917 0.133767405586064 0.207491289198608 0.0783450704225352],...
    'FontSize',12,'box','off');
% legend('boxoff')

P.plots_dir = 'gradsim_paper/fft/join';
P.plot_name = 'fvsrginset';
P.fig_handle = figx;
P.save_plot();

