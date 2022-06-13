%________________________________________________________________________
% 
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 
%________________________________________________________________________

clear;
% close all;

% EXPERIMENTAL INFO:

% parameters
datadir = 'r2l_302_c_noe';
dump_list = 0:1:10;
dataformat = 'mat';

% save parameters
save_format = {'png'};
save_flag = 1;
save_plot_flag = 0;
plots_dir = ['r2plot/transp/',datadir];

% load classes
O = OsirisDenormalizer('plasmaden',2e14,'datadir',datadir,...
    'property','density',...
    'dataformat',dataformat,'species','proton_beam');
O.getdata(); 
O.property = 'raw';

P = Plotty('datadir',datadir,'plasmaden',2e14,'property_plot','wakefields',...
    'save_format',save_format,'save_flag',save_plot_flag,'plots_dir',plots_dir);

r_lims = linspace(O.nr(1),O.nr(end),length(O.nr)+1); % 1/kp
r_edges = [-inf,r_lims(2:end-1),inf];
z_lims = linspace(O.nz(1),O.nz(end),length(O.nz)+1);
z_edges = [-inf,z_lims(2:end-1),inf];

% normalize to charge
% analytical charge fraction

for n = 1:length(dump_list)
    close all
    O.dump = dump_list(n);

    % get data
    O.raw_dataset = 'x'; O.direction = 'z'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'p'; O.direction = 'z'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'p'; O.direction = 'r'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
% 
%     [~,~,~,zbins,rbins] = histcounts2(O.ntime + ...
%         O.n_simulation_window - O.nz_raw,O.nr_raw,z_edges,r_edges);
    [~,~,~,zbins,rbins] = histcounts2(abs(O.nz_raw),O.nr_raw,z_edges,r_edges);
    pr = [O.npr_raw./O.npz_raw;0;0];
    nan_ind = isnan(pr);
    pr(nan_ind) = 0;
    rbins0 = [rbins;1;length(r_edges)-1];
    zbins0 = [zbins;1;length(z_edges)-1]; 
    phasematrix = accumarray([rbins0,zbins0],pr); 
    
    phasespace = fliplr(phasematrix);
    
    xlims = [z_lims(1),z_lims(end)];
    ylims = [-r_lims(end),r_lims(end)];
   % P.plot_field_density('field_plot',phasematrix,'r_plot',10*ylims,'z_plot',xlims);
    
    dump_char = sprintf('%06.6d',O.dump);
    save_name = ['PHASESPACE-proton_beam-p2-',dump_char,'','.mat'];
    partialpath_handle = what(O.datadir);
   
    x1_axis = O.nz; x2_axis = O.nr; 
    axis1 = [O.nz(1),O.nz(end)]; axis2 = [O.nr(1),O.nr(end)]; 
    time = O.ntime;

    if save_flag
        save_dir = [partialpath_handle(1).path,'/MAT/PHASESPACE','/proton_beam','/p2'];
        if ~isfolder(save_dir); mkdir(save_dir); end
        save([save_dir,'/',save_name],'phasespace','time','axis1','x1_axis','axis2','x2_axis')
    end

    O.progress_dump('building momentum: ',n,length(dump_list));
 
end % for dump