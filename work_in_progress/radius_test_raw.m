%________________________________________________________________________
% Script to measure the radius of the proton bunch from the raw
% macroparticles.
%
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 23/03/2021
%________________________________________________________________________

close all;

datadir = 'gradsim_shortbunch';
useavg = 1;
dataformat = 'h5';

dump_list = 0:1:100;
save_flag = 1;

plasmaden = 1.81e14;

plots_dir = ['gradsim/radius_check'];
plot_name = ['g0sb'];

% gradsim sim window = 21cm
% plasma wavelength = 0.25 cm
% 1st microbunch = 0.25
% nbins = 21/0.25 = 84;
nbins = 1;




O = OsirisDenormalizer('datadir',datadir,'dataformat',dataformat,'useAvg',useavg,...
    'property','density',...
    'dump',0,'plasmaden',1.81e14);

Or = OsirisDenormalizer('datadir',datadir,'dataformat',dataformat,'useAvg',useavg,...
    'property','density',...
    'dump',0,'plasmaden',1.81e14);

P = Plotty('plasmaden',plasmaden,'plots_dir',plots_dir,'plot_name',plot_name,...
    'save_flag',save_flag);

% raw
Or.property = 'raw'; Or.dataformat = 'h5';
Or.raw_dataset = 'x'; Or.direction = 'r'; Or.getdata(); Or.assign_raw();
Or.direction = 'z'; Or.getdata(); Or.assign_raw();
Or.raw_dataset = 'q'; Or.getdata(); Or.assign_raw();
Or.raw_dataset = 'ene'; Or.getdata(); Or.assign_raw();
Or.raw_dataset = 'p'; Or.direction = 'r'; Or.getdata(); Or.assign_raw();
Or.direction = '\theta'; Or.getdata(); Or.assign_raw();

% initialize vaiables
r_rms = zeros(length(dump_list),nbins);
avg_z = zeros(length(dump_list),1);
r_den = zeros(length(dump_list),1);
r_push = zeros(length(dump_list),1);
% loop over dumps
for n = 1:length(dump_list)
    
    O.dump = dump_list(n);
    % density
    O.property = 'density';
    O.getdata(); O.assign_density(); O.denorm_density(); O.denorm_distance();
    
    % raw
    O.property = 'raw'; O.dataformat = 'h5';
    O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();
    r = O.denorm_distance(O.nr_raw);
    O.direction = 'z'; O.getdata(); O.assign_raw();
    z = O.denorm_distance(O.nz_raw);
    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
    q = O.q_raw;
    O.raw_dataset = 'ene'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'p'; O.direction = 'r'; O.getdata(); O.assign_raw();
    
    
    % rms
    % binning
%     [z_bin,z_edges] = discretize(z,nbins);
    [z_bin,z_edges] = discretize(z,O.dtime + O.simulation_window - [0.1 0]);
    avg_z(n) = O.dtime;
    %avg_z(n) = sum(z.*q)/sum(q);
    
    for b = 1:1
        z_ind = (z_bin == b);
        r_bin = r(z_ind);
        q_bin = q(z_ind);
        r_rms(n,b) = sqrt(sum((r_bin.^2).*q_bin)/sum(q_bin))/sqrt(2);
    end
    
    % trimming
    z_ind = O.z > z_edges(end-1);
    
    
    % density
    dz = O.z(2) - O.z(1);
    trans_density_profile = dz*sum(O.proton_beam(:,z_ind),2);
    trans_density_profile_mirrored = [flipud(trans_density_profile);trans_density_profile];
    r_fit = [-fliplr(O.r),O.r];
    
    f = fit(r_fit',trans_density_profile_mirrored,'gauss1');
    r_den(n) = f.c1/sqrt(2);
%     plot(r_fit,trans_density_profile_mirrored);
    
    
    % push
    
    if n == 1
        zraw_ind = Or.denorm_distance(Or.nz_raw) > z_edges(end-1);
        Or.nr_raw = Or.nr_raw(zraw_ind);
        Or.nE_raw = Or.nE_raw(zraw_ind);
        Or.npr_raw = Or.npr_raw(zraw_ind);
        Or.npth_raw = Or.npth_raw(zraw_ind);
        Or.q_raw = Or.q_raw(zraw_ind);
    end
    r_new = Or.denorm_distance(Or.charge_pusher(Or,avg_z(n) - avg_z(1)));
    r_push(n) = sqrt(sum((r_new.^2).*Or.q_raw)/sum(Or.q_raw))/sqrt(2);
    
    O.progress_dump('calc radius ',n,length(dump_list));
end

z_plot = z_edges + (z_edges(2) - z_edges(1))/2;
z_plot(end) = [];
% beta = 0.0002^2/(3.6e-6/426.439)
r_beta = 0.02*sqrt(1 + ((avg_z - avg_z(1))/100).^2/(4.7382*sqrt(1))^2);

figure
hold on
plot(avg_z,r_rms(:,end))
plot(avg_z,r_beta)
plot(avg_z,r_den)
plot(avg_z,r_push)
hold off

xlim([0 1000])

xlabel('z (cm)');
ylabel('radius (cm)');
% legend('beta','push','location','best');
legend('rms','beta','den','push','location','best');

P.fig_handle = gcf;
P.save_plot();

%
% p = O.proton_beam;
%
% pm = [flipud(p);p];
%
% figure(1);
% imagesc(pm);
%
% % lo = pm(:,35300);
% lo = sum(pm,2);
% x = [-fliplr(O.r),O.r];
%
% figure(2);
% plot(x,lo)
%
% f = fit(x(:),lo(:),'gauss1');
%
% cr = O.radial_integration(O.r,O.z,p);