%________________________________________________________________________
% Create a time density profile of the proton bunch, similar to what is obtained in
% the experiment (density profile at one point in space, while time
% advances, instead of a profile frozen in time while z varies).
%
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
% close all;

% EXPERIMENTAL INFO:
% maximum delay in the experiment: 618.2097 ps --> 18.5335 cm
% longitudinal pixels after the seeding: 1500
% to get 1500 pixel res., binsize = 0.41214 ps

% exp. total transverse distance = 672*0.00217 = 1.4582 cm
% exp. half distance ~= 0.7291 cm in ~336 pixels
% --> exp pixels in 0.16 cm ->
% transverse pixels total = 672 --> 73.73 ~ 74
% exp pixels in 0.2/0.00217 = 92.17 ~ 92

% parameters

datadir = 'gp20';
measurement_point = 1350; % in cm, from plasma beginning
dump_list = 95;
binsize = 0.41214;
trans_lims = linspace(0,0.2,92); %cm (74 for 0.16 cm (see above))
% trans_lims = linspace(0,0.01,40); %cm (74 for 0.16 cm (see above))
xi_range = [18.5335,0];
dataformat = 'h5';
initial_protons = 3e11;
seeding_position = 3.8074;
sigma_z = 6.98516; % cm

% slit parameters
slit_flag = 0;
dy = (74/2)*1e-4; % 74 um width
n_theta_divisions = 30; % number of divisions in theta to create particles in 3D, that will be considered in the slit

% save parameters
save_format = {'png','eps'};
save_flag = 1;
save_plot_flag = 0;

if slit_flag
    slit_text = '_slit';
else
    slit_text = '';
end

plots_dir = ['gradsim/densitytimeprofile/',datadir];
save_dir = ['loading_files/gradsim/densitytimeprofile',slit_text,'/',datadir];

% load classes
O = OsirisDenormalizer('plasmaden',1.81e14,'datadir',datadir,...
    'property','raw','xi_range',xi_range,...
    'dump',0,'dataformat',dataformat,'species','proton_beam');
P = Plotty('plasmaden',1.81e14,'property_plot','density',...
    'save_format',save_format,'save_flag',save_plot_flag,'plots_dir',plots_dir,...
    'dump_list',0);

O.raw_dataset = 'q'; O.getdata(); O.assign_raw();

% normalize to charge
% analytical charge fraction
bunch_long_fraction = normcdf(seeding_position,0,sigma_z);
bunch_back_fraction = normcdf(seeding_position - xi_range(1),0,sigma_z);
bunch_fraction_inside_simulation_window = bunch_long_fraction - bunch_back_fraction;
bunch_protons_inside_simulation_window = initial_protons*bunch_fraction_inside_simulation_window;
charge_norm_factor = bunch_protons_inside_simulation_window./sum(O.q_raw);

for n = 1:length(dump_list)
    O.dump = dump_list(n);
    
    % get data
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
%     O.useAvg = 1;
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
    
    
    if slit_flag
        % calculate the value of theta (angle from x plane to slit, with origin at the axis)
        % for each macroparticle.
        th_end_points = asin(dy./new_r);
        % imaginary values would mean particle ring is fully inside the slit,
        % then replace it by the maximum value possible pi/2
        i_complex = (imag(th_end_points) ~= 0);
        th_end_points(i_complex) = pi/2;
        %    th_end_points(i_complex) = [];
        % multiply q by the fraction theta/(pi/2) a quarter of a full rotation
        th_frac = th_end_points/(pi/2);
        q_temp = O.q_raw.*th_frac;
        
        eq_spaced_v = linspace(0,1,n_theta_divisions);
        th_mat = eq_spaced_v.*(th_end_points);
        new_x = new_r.*cos(th_mat);
        new_r = new_x(:);
        
        new_q = q_temp.*(ones(1,n_theta_divisions)./n_theta_divisions);
        new_q = new_q(:);
        
        new_z = z.*(ones(1,n_theta_divisions));
        new_z = new_z(:);
        
    else
        new_q = O.q_raw(:);
        new_z = z(:);
    end
    
    delay_time = (O.dtime + O.simulation_window - new_z)/O.c_cm*1e12; % ps
    t_simulation_window = (O.simulation_window)/O.c_cm*1e12; % ps
    
    trim_window_size = (max(O.z) - min(O.z))/O.c_cm*1e12;
    chargematrix = zeros(length(trans_lims) - 1, ceil(trim_window_size/binsize));
    
    for rr = 1:(length(trans_lims)-1)
        ir = (new_r >= trans_lims(rr)) & (new_r < trans_lims(rr+1));
        if sum(ir) == 0; continue; end
        
        % select only those particles inside the transverse limits
        q_r = new_q(ir);
        delay_time_r = delay_time(ir);
        [N,edges,ind_bin] = histcounts(delay_time_r,0:binsize:trim_window_size);
        A = accumarray(ind_bin+1,delay_time_r);
        % bin
        ind_bin = ceil(delay_time_r/binsize);
        
        % sort the indeces of the binning and then sort the
        % elements in q
        [ind_sort,ind_order] = sort(ind_bin);
        q_sort = q_r(ind_order);
        
        % say how many there is in each bin in a cumulative way,
        % to use that as indices
        % for q, to quickly build up the indexed density along z
        ind_sum = [0,0,histcounts(ind_sort,1:1:max(ind_bin),'Normalization','cumcount')];
        
        % go through each element in bz and arrange that charge
        % according to the bins, using the indeces from histconts
        for bz = 1:max(ind_bin)
            chargematrix(rr,bz) = sum(q_sort(ind_sum(bz)+1:ind_sum(bz+1)));
        end
        O.progress_dump(['building density',slit_text],rr,length(trans_lims)-1)
        
    end
    
    chargematrix = fliplr(chargematrix)*charge_norm_factor;
    r_in = trans_lims(1:end-1);
    r_ex = trans_lims(2:end);
    
    if slit_flag
        densitymatrix = chargematrix;
    else
        
        ringvolume = binsize*pi*(r_ex.^2-r_in.^2);
        densitymatrix = (chargematrix./ringvolume');
    end
    
    tbin = linspace(trim_window_size,0,ceil(trim_window_size/binsize));
    xlims = [tbin(1),tbin(end)];
    ylims = [-trans_lims(end),trans_lims(end)];
    P.plot_field_density('density_plot',densitymatrix,'r_plot',10*ylims,'z_plot',xlims);
    
    
    P.plot_name = ['n',num2str(O.dump),slit_text];
    P.fig_handle = gcf;
    P.save_plot();
    save_name = ['n',num2str(O.dump),'_m',num2str(measurement_point),slit_text,'.mat'];
    if ~isfolder(save_dir); mkdir(save_dir); end
    if save_flag
        save([save_dir,'/',save_name],'chargematrix','densitymatrix','xlims','ylims')
    end
    int_trans_lims = (r_in + r_ex)/2;
    charge = O.cylindrical_integration(int_trans_lims,fliplr(tbin),densitymatrix,'sum')/initial_protons
    
end % for dump