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
datadir = 'DW_lcode_x20_pi';
dump_list = 100:1:100;
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

    % get particle data
    O.raw_dataset = 'x'; O.direction = 'z'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();
    % O.raw_dataset = 'p'; O.direction = 'z'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'p'; O.direction = 'r'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'p'; O.direction = 'z'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'p'; O.direction = '\theta'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
%     O.raw_dataset = 'ene'; O.getdata(); O.assign_raw();
    
    % get field data
    O.property = 'fields'; O.wakefields_direction = 'trans';
    O.getdata(); O.assign_fields();
    transfield_now = O.ntransfield;

%     O.dump = O.dump + 1;
    O.getdata(); O.assign_fields();
    transfield = O.ntransfield;

    i_p = (O.npr_raw > 0); 
%     i_z = abs(O.nz_raw) > 200 & abs(O.nz_raw) < 215;
    i_r = O.nr_raw < 0.5322;
    i_p = i_p & i_r;

    pr = O.npr_raw(i_p);
    pz = O.npz_raw(i_p);
    pt = O.npth_raw(i_p);
    p = sqrt(pr.^2 + pz.^2 + pt.^2);
    vr = sign(pr).*sqrt((pr/1836.15267).^2./(1 + (pr/1836.15267).^2));
%     vr = sqrt(2)*pr./(E + 1);
    t = 0.02./vr; 
    
    % beta = sqrt(1 - 1/gamma^2)
    % pr = gamma*beta*c*mass

    r = O.nr_raw(i_p);
    z = O.nz_raw(i_p);
    q = O.q_raw(i_p);

%     E = O.nE_raw(i_p);

    es = zeros(1,length(r));


    [~,~,~,zbins,rbins] = histcounts2(abs(z),r,z_edges,r_edges);
    idx = sub2ind(size(transfield_now), rbins, zbins);
    trf_q = transfield_now(idx).*1; % q
    max_r = length(O.nr);


    
    for pp = 1:max_r
        
        pr_n = pr + t.*trf_q;
        
        pr = pr_n;
        
        vr = sign(pr).*sqrt(((pr/1836.15267).^2)./(1 + (pr/1836.15267).^2));
%         vr = sqrt(2)*pr./(E + 1);

        t = 0.02./vr; 
        r = r + sign(pr).*0.02;

        [~,~,~,zbins,rbins] = histcounts2(abs(z),r,z_edges,r_edges);
%         [~,~,~,zbins,rbins] = histcounts2(O.ntime + ...
%         O.n_simulation_window - z,r,z_edges,r_edges);

        idx = sub2ind(size(transfield_now), rbins, zbins);
        trf_q = transfield_now(idx).*1; % q

        disp(pp)


    end

    a = 1;


% 
%     [~,~,~,zbins,rbins] = histcounts2(O.ntime + ...
%         O.n_simulation_window - O.nz_raw,O.nr_raw,z_edges,r_edges);

%     [~,~,~,zbins,rbins] = histcounts2(abs(O.nz_raw),O.nr_raw,z_edges,r_edges);
    pr_q = [O.q_raw.*O.npr_raw;0;0];
    rbins0 = [rbins;1;length(r_edges)-1];
    zbins0 = [zbins;1;length(z_edges)-1]; 
    phasematrix = accumarray([rbins0,zbins0],pr_q); 
    
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