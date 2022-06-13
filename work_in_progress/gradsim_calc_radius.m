%________________________________________________________________________
%
% Script to calculate the radius of the proton bunch.
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

% save directories
clear;
close all;

% paramaters
plasmaden = 1.81e14;
datadir = 'g0';
dump_list = 0:1:100;
binsize = OsirisDenormalizer.ps2cm(0.41214);
trans_lims = linspace(0,0.16,74); %cm (74 for 0.16 cm (see above))
prop_distances = 9.958450708031729*dump_list;

r_in = trans_lims(1:end-1);
r_ex = trans_lims(2:end);
ringvolume = binsize*pi*(r_ex.^2 - r_in.^2);

xi_range = [0.005,0];
dataformat = 'h5';


% get macroparticle info

O = OsirisDenormalizer('plasmaden',plasmaden,'datadir',datadir,...
    'property','raw','xi_range',xi_range,...
    'dump',100,'dataformat',dataformat,'species','proton_beam');
P = Plotty('plasmaden',plasmaden,'datadir',datadir,...
    'species','proton_beam',...
    'property_plot','density',...
    'plots_dir','radius_test');


O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();
O.raw_dataset = 'x'; O.direction = 'z'; O.getdata(); O.assign_raw();
O.raw_dataset = 'p'; O.direction = 'r'; O.getdata(); O.assign_raw();
O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
O.raw_dataset = 'ene'; O.getdata(); O.assign_raw();

% to get simulation window size, needed for the distance
O.dataformat = 'mat'; % done for gm20
O.property = 'density'; O.getdata(); O.trim_data(); O.denorm_distance();
O.dataformat = 'h5';


% manual trim of raw data
ixi = (O.denorm_distance(O.nz_raw) - O.dtime) < xi_range(1);
z = O.denorm_distance(O.nz_raw(ixi));
r = O.denorm_distance(O.nr_raw(ixi));
q = O.q_raw(ixi);

O.npr_raw = O.npr_raw(ixi);
O.nr_raw = O.nr_raw(ixi);
O.nE_raw = O.nE_raw(ixi);

    % delay_time = (O.dtime + O.simulation_window - z)/O.c_cm*1e12; % ps
    % t_simulation_window = (O.simulation_window)/O.c_cm*1e12; % ps
%     z_pos = (O.dtime + O.simulation_window - z); %z with respect to time
    z_pos = z - O.dtime;
    trim_window_size = max(O.z) - min(O.z);
    rms_raw = zeros(1,length(dump_list));

for n = 1:length(dump_list)
        
%     new_r = O.denorm_distance(O.charge_pusher(O,prop_distances(n)));
    new_r = r;
    rms_raw(n) = sqrt(sum((new_r.^2).*q)/sum(q))/sqrt(2);
    
    O.progress_dump('radius calculation',n,length(dump_list));
    
end


rms_theorie = 0.02*sqrt(1+(prop_distances/100).^2/4.7381^2);

hold on
plot(prop_distances,rms_raw)
plot(prop_distances,rms_theorie)
hold off

legend({'rms raw','rms theorie \beta = m'})

xlabel('z (m)');
ylabel('rms radius (cm)')

P.plot_name = 'beta67';
P.fig_handle = gcf;
P.save_plot();


