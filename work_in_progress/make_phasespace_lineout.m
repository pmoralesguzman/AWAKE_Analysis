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
dump_list = 0:1:0;
dataformat = 'mat';

% save parameters
save_format = {'png'};
save_flag = 1;
save_plot_flag = 0;
plots_dir = ['r2plot/transp/',datadir];

lineout_direction = 'long';

% load classes
O = OsirisDenormalizer('plasmaden',2e14,'datadir',datadir,...
    'property','density',...
    'dataformat',dataformat,'species','proton_beam');
O.getdata(); 
O.property = 'raw';

P = Plotty('datadir',datadir,'plasmaden',2e14,'property_plot','phasespace',...
    'lineout_direction',lineout_direction,...
    'save_format',save_format,'save_flag',save_plot_flag,'plots_dir',plots_dir);

z_lims = linspace(O.nz(1),O.nz(end),length(O.nz)+1);
z_edges = [-inf,z_lims(2:end-1),inf];

legs = {'$p+$ (away from axis)','$p-$ (towards axis)'};
load("colororder_defaultblack.mat"); %corder_default
linestyles = {'-',':','--','-.','-','-','-','-.','--',':'};


% normalize to charge
% analytical charge fraction

for n = 1:length(dump_list)
    
    O.dump = dump_list(n);

    % get data
    O.raw_dataset = 'x'; O.direction = 'z'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();
    % O.raw_dataset = 'p'; O.direction = 'z'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'p'; O.direction = 'r'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
% 
%     [~,~,~,zbins,rbins] = histcounts2(O.ntime + ...
%         O.n_simulation_window - O.nz_raw,O.nr_raw,z_edges,r_edges);

    i_p = O.npr_raw > 0;
    i_n = O.npr_raw < 0;
    i_r = O.nr_raw < 0.5323;
    i_p = i_p & i_r;
    i_n = i_n & i_r;

    [~,~,zbins_p] = histcounts(abs(O.nz_raw(i_p)),z_edges);
    pr_p = [O.npr_raw(i_p);0;0];
    
    zbins0_p = [zbins_p;1;length(z_edges)-1]; 
    p2lineout_p = accumarray(zbins0_p,pr_p); 

    [~,~,zbins_n] = histcounts(abs(O.nz_raw(i_n)),z_edges);
    pr_n = [O.npr_raw(i_n);0;0];

    zbins0_n = [zbins_n;1;length(z_edges)-1]; 
    p2lineout_n = accumarray(zbins0_n,pr_n);

    z_plot = (z_lims(1:end-1) + z_lims(2:end))/2;
%     close all
%     P.plot_lineout2('lineout_plot',p2lineout_p,'z_plot',z_plot);
%     P.plot_handle.LineStyle = linestyles{1};
%     P.plot_handle.Color = corder_defblack(1,:);
%     P.plot_handle.LineWidth = 1;
% 
%     hold on
% 
%     P.plot_lineout2('lineout_plot',abs(p2lineout_n),'z_plot',z_plot);
%     P.plot_handle.LineStyle = linestyles{3};
%     P.plot_handle.Color = corder_defblack(3,:);
%     P.plot_handle.LineWidth = 1;
% 
%     hold off
% 
%     legend(legs{:},'location','best','interpreter','latex','FontSize',P.plot_fontsize-6)
% 
% 
%     drawnow; pause(0.5);

    dump_char = sprintf('%06.6d',O.dump);
    save_name = ['PHASESPACE-proton_beam-p2l-',dump_char,'','.mat'];
    partialpath_handle = what(O.datadir);
   
    x1_axis = O.nz; x2_axis = O.nr; 
    axis1 = z_plot; axis2 = [O.nr(1),O.nr(end)]; 
    time = O.ntime;

    if save_flag
        save_dir = [partialpath_handle(1).path,'/MAT/PHASESPACE','/proton_beam','/p2'];
        if ~isfolder(save_dir); mkdir(save_dir); end
        save([save_dir,'/',save_name],'p2lineout_p','p2lineout_n','time',...
            'axis1','x1_axis','axis2','x2_axis')
    end

    O.progress_dump('building momentum: ',n,length(dump_list));
 
end % for dump