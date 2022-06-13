%________________________________________________________________________
% Charge or density transverse profile along z.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/10/2020
%________________________________________________________________________


% data directory
datadirs    = {'gm20'};
dataformat  = 'mat';
useAvg      = false;
dump_list   = 0:5:100;

% save directory
plots_dir       = ['gradsim_paper/charge_trans_evolution/',datadirs{1}];
plot_name_suffix= [''];
save_format     = {'png','eps','fig'};
save_plot       = 0;

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';


% analysis
plasma_radius   = 0.15; % cm
xi_range        = [21 0.74];
scan_type       = 'slice'; % slice, cumulative
on_axis         = 'int'; % int, sum, intw, lineout

% calculated parameters
trans_lims = 2*3.9598e-04:2*3.9598e-04:plasma_radius;

% initialize variables
charge_in_slice = ones(length(trans_lims),length(dump_list));
prop_distance_m = zeros(1,length(dump_list));

AC = AwakeFFT(...
    'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'xi_range',xi_range,'trans_lims',trans_lims,...
    'scan_type',scan_type,'on_axis',on_axis);

P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'save_format',save_format,'save_flag',save_plot);

for d = 1:length(datadirs)
    AC.datadir = datadirs{d};

    for n = 1:length(dump_list)
        AC.dump = dump_list(n);
        AC.fft_dataload();
        prop_distance_m(n) = AC.propagation_distance/100; % propagation distance in m

        switch AC.property
            case 'density'
                charge_in_slice(:,n) = AC.dz*sum(AC.fft_densitymatrix,2);
        end % switch property 
        AC.progress_dump('dump',n,length(dump_list));
    end %  for dump list
    
    fct = figure(1);
    [X,Y] = meshgrid(prop_distance_m,trans_lims*10);
    waterfall(X',Y',charge_in_slice');
    
    ylabel('r (mm)');
    xlabel('z (m)');
    zlabel('charge density (protons/cm)');
    
    drawnow;
    P.plot_name = ['charge_trans',AC.datadir,'n',num2str(dump_list(end))];
    P.fig_handle = fct;
    P.save_plot();
    
end % for datadirs







