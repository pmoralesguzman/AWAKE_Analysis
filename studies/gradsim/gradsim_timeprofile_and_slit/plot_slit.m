%________________________________________________________________________
% Plot a random sample of particles as they would fill 3D spaces in a
% ring-like manner, in order to simulate the same effect as adding a slit
% through which the proton beam goes, as in the experiment. 
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 30/11/2020
%________________________________________________________________________

clear;
close all;

% maximum delay in the experiment: 618.2097 ps --> 18.5335 cm
% longitudinal pixels after the seeding: 1500
% to get 1500 pixel res., binsize = 0.41214 ps

% exp. total transverse distance = 672*0.00217 = 1.4582 cm
% exp. half distance ~= 0.7291 cm in ~336 pixels
% --> exp pixels in 0.16 cm ->
% transverse pixels total = 672 --> 73.73 ~ 74

% parameters

datadir = 'gm10';
measurement_point = 1350; % in cm, from plasma beginning
dump = 100;
binsize = 0.41214;
trans_lims = linspace(0,0.3,138); %cm (74 for 0.16 cm (see above))
xi_range = [18.5335,0];
dataformat = 'h5';

% get macroparticle info

O = OsirisDenormalizer('plasmaden',1.81e14,'datadir',datadir,...
    'property','raw','xi_range',xi_range,...
    'dump',dump,'dataformat',dataformat,'species','proton_beam');
P = Plotty('plasmaden',1.81e14,'plots_dir','slit','plot_name','3D','save_format','png');
O.raw_dataset = 'x'; O.direction = 'z'; O.getdata(); O.assign_raw();
O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();


O.raw_dataset = 'p'; O.direction = 'z'; O.getdata(); O.assign_raw();
% pz = O.denorm_distance(O.npz_raw);
O.raw_dataset = 'p'; O.direction = 'r'; O.getdata(); O.assign_raw();
% pr = O.denorm_distance(O.npr_raw);

O.raw_dataset = 'ene'; O.getdata(); O.assign_raw();
O.raw_dataset = 'q'; O.getdata(); O.assign_raw();


% to get simulation window size, needed for the distance
O.dataformat = 'mat'; % done for gm20
O.property = 'density'; O.getdata(); O.trim_data(); O.denorm_distance();
O.dataformat = 'h5';

% manual trim of raw data

ixi = (O.dtime + O.simulation_window - O.denorm_distance(O.nz_raw)) < xi_range(1);

z = O.denorm_distance(O.nz_raw(ixi));
r = O.denorm_distance(O.nr_raw(ixi));
O.npr_raw = O.npr_raw(ixi);
O.nr_raw = O.nr_raw(ixi);
O.q_raw = O.q_raw(ixi);
O.nE_raw = O.nE_raw(ixi);

% front and back propagate

prop_distances = measurement_point + O.simulation_window - z;

new_r = O.denorm_distance(O.charge_pusher(O,prop_distances));

dy = (200/2)*1e-4; % 74 um width
dy_slit = (74/2)*1e-4*10; % 74 um width (just for the plot, plotted in mm)
n_points = 10;

eq_spaced_v = linspace(0,1,n_points);
th_end_points = asin(dy./new_r);
i_complex = (imag(th_end_points) ~= 0);
th_end_points(i_complex) = pi/2;
q_temp = O.q_raw.*th_end_points;

th_mat = eq_spaced_v.*(th_end_points);
new_x = 10*new_r.*cos(th_mat);
% new_r = new_x(:);

new_y = 10*new_r.*sin(th_mat);

new_q = q_temp.*(ones(1,n_points)./n_points);
new_q = new_q(:);

new_z = z.*(ones(1,n_points));
% new_z = new_z(:);


delay_time = (O.dtime + O.simulation_window - new_z)/O.c_cm*1e12; % ps
t_simulation_window = (O.simulation_window)/O.c_cm*1e12; % ps

i_in_slit = find(new_r < 2*dy);

%%
N = 15;
% random_picker = randperm(N,N); %length(O.q_raw)
% random_first_step = randperm(length(O.q_raw),N);

% random_first_step = [i_in_slit(randi(length(i_in_slit))),i_in_slit(randi(length(i_in_slit))),randperm(length(O.q_raw),N-2)];
random_first_step = i_in_slit(randi(length(i_in_slit),[N 1]))';
random_picker = repmat(random_first_step,n_points,1); random_picker = random_picker(:);
random_picker30 = repmat((1:n_points)',N,1); %randi(2,1,N);

p_x = new_x(sub2ind(size(new_x),random_picker,random_picker30));
p_y = new_y(sub2ind(size(new_y),random_picker,random_picker30));
p_z = delay_time(sub2ind(size(new_z),random_picker,random_picker30));

cmp2D = colormap(parula(N));
cmp3D = repelem(cmp2D,n_points*ones(N,1),1);

close;

plot_x = [p_x;p_x;-p_x;-p_x];
plot_y = [p_y;-p_y;p_y;-p_y];
plot_z = [p_z;p_z;p_z;p_z];

plot_cmp = [cmp3D;cmp3D;cmp3D;cmp3D];
S = 50;
% 

fig_2D = figure(2);
scatter(delay_time(random_first_step),new_x(random_first_step),S,cmp2D,'filled');
hold on
yline(dy_slit,'k')
hold off
ax2D = gca; ax2D.XDir = 'reverse';
xlabel('z-ct (ps)'); ylabel('r (mm)');
ylim([0,max(new_x(random_first_step))])
P.fig_handle = fig_2D;
P.plot_name = 'slit2D_100';
P.save_plot();


fig_3D = figure(3);
hold on
scatter3(plot_z,plot_x,plot_y,S,plot_cmp,'filled')
X = [min(plot_x);max(plot_x);max(plot_x);min(plot_x);min(plot_x)];
Z = [min(p_z);min(p_z);max(p_z);max(p_z);min(p_z)];
Y = [0;0;0;0;0];

patch(Z,X,Y+dy_slit,'k','FaceAlpha',0.3);   % draw a square in the xy plane with z = 0
patch(Z,X,Y-dy_slit,'k','FaceAlpha',0.3); % draw a square in the xy plane with z = 1
for k = 1:length(X)-1
    patch([Z(k);Z(k)],[X(k);X(k)],[-dy_slit;dy_slit],'k');
end
hold off
ax3D = gca;
ax3D.XDir = 'reverse';
ylabel('x (mm)'); zlabel('y (mm)'); xlabel('z-ct (ps)'); 

view([61.7650,7.4127])
P.fig_handle = fig_3D;
P.plot_name = 'slit3D_100_1';
drawnow;
% a = 1;
P.save_plot();

view([90,0])
P.plot_name = 'slit3D_100_2';
drawnow;
% a = 1;
P.save_plot();

