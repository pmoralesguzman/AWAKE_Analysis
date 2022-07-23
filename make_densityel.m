function [] = make_densityel(datadirs,dumplist,speciesname,xistepsize)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% close all;

% EXPERIMENTAL INFO:

% parameters
%datadir = 'fdr_26';
datadir = datadirs;
dump_list = dumplist;
dataformat = 'mat';
% species_name = '';
species_name = speciesname;
xi_step_size = xistepsize;

% save parameters
save_format = {'png'};
save_flag = 1;
save_plot_flag = 0;
plots_dir = ['r2plot/make_density/',datadir];

% load classes
O = OsirisDenormalizer('plasmaden',2e14,'datadir',datadir,...
    'property','fields',...
    'dataformat',dataformat,'species',species_name);
O.getdata();

O.property = 'raw';

P = Plotty('datadir',datadir,'plasmaden',2e14,'property_plot','density',...
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

    O.property = 'fields';

    O.getdata();

    O.property = 'raw';

    % get data
    O.raw_dataset = 'x'; O.direction = 'z'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();
    % O.raw_dataset = 'p'; O.direction = 'z'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'w'; O.getdata(); O.assign_raw();
    %
    %     [~,~,~,zbins,rbins] = histcounts2(O.ntime + ...
    %         O.n_simulation_window - O.nz_raw,O.nr_raw,z_edges,r_edges);
    [~,~,~,zbins,rbins] = histcounts2(abs(O.nz_raw),O.nr_raw,z_edges,r_edges);


    i_qp = O.q_raw > 0;
    i_qn = ~i_qp;

    r_in = r_lims(1:end-1);
    r_ex = r_lims(2:end);

        
    q = O.denorm_lcode_charge([O.w_raw;0;0],xi_step_size);
    rbins0 = [rbins;1;length(r_edges)-1];
    zbins0 = [zbins;1;length(z_edges)-1];
    chargematrix = accumarray([rbins0,zbins0],q);
    
    O.units = 'cm';
    binsize = O.denorm_distance(0.02);
    ringvolume = binsize*pi*...
        (O.denorm_distance(r_ex).^2 - O.denorm_distance(r_in).^2);
    densitymatrix = 1*(chargematrix./ringvolume')./O.plasmaden; %1836.15267343

    density = fliplr(densitymatrix);

    xlims = [z_lims(1),z_lims(end)];
    ylims = [-r_lims(end),r_lims(end)];
%     P.plot_field_density('density_plot',abs(densitymatrix),'r_plot',10*ylims,'z_plot',xlims);

    dump_char = sprintf('%06.6d',O.dump);
    save_name = ['charge-',species_name,'-',dump_char,'','.mat'];
    partialpath_handle = what(O.datadir);

    x1_axis = O.nz; x2_axis = O.nr;
    axis1 = [O.nz(1),O.nz(end)]; axis2 = [O.nr(1),O.nr(end)];
    time = O.ntime;

    if save_flag
        save_dir = [partialpath_handle(1).path,'/MAT/DENSITY','/',species_name,'/charge'];
        if ~isfolder(save_dir); mkdir(save_dir); end
        save([save_dir,'/',save_name],'density','time','axis1','x1_axis','axis2','x2_axis')
    end


    O.progress_dump('building density: ',n,length(dump_list));

end % for dump
end