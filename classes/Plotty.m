%__________________________________________________________________________
% Class with the main plotting function for the AWAKE Osiris Analysis
% Matlab Package
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/07/2020
%__________________________________________________________________________

% Input:
% - plots directory, plot name, saving instructions
%
% Output:
% - saved plot (no variable)
%
% Methods
% - save_plot: created saving directory and saves plot according to
% instrucions
%

classdef Plotty < handle & OsirisDenormalizer

    properties(Constant,Hidden)


    end % constant properties


    properties(Hidden)

        % waterfall input
        dump_list;
        do_plot;

        % waterfall output
        waterfall_xi;
        waterfall_z;

        % field density plot input
        save_movie_struct;

        % flags
        title_flag;
        make_pause;

        % plot scale
        plot_scale;

        % plot options
        fig_number;
        plot_fontsize = 18;


        % field density plot
        ylinepos;
        xlinepos;
        mirror_flag;

        ax_field;
        ax_density;
        ax_phasespace;
        tile_handle;
        include_lineout;


        ylineflag;
        xlineflag;
        lines_2D_flag;
        include_phasespace;

        % movie
        frame_counter = 0;

        % temp
        z_out;
        a = 0;

    end % hidden properties

    properties

        % input, description in parser
        plots_dir;
        plot_name;
        save_flag;
        save_format;
        exlegends;
        extitles;

        extradatadir;

        

        % plot field density input
        property_plot;
        denormalize_flag;
        create_movie;
        include_field_lineout;
        include_density_profile;
        include_phasespace_profile;
        field_geometry;

        % output
        waterfall_mat;
        waterfall_handle;
        fig_handle;plotty
        plot_handle;
        

        struct_movie;

        plot_density_lims;
        plot_field_lims;

    end % properties

    methods

        % --CONSTRUCTOR--

        function obj = Plotty(varargin)


            % parse input to load files
            p = inputParser;

            p.addParameter('extradatadir', '', @(x) ischar(x));

            % directory where the plot should be saved
            p.addParameter('plots_dir', 'plots', @(x) ischar(x));
            p.addParameter('plot_name', 'plot', @(x) ischar(x));


            % for the case of plots that need a dump list (ex.: waterfall)
            p.addParameter('dump_list', 1:1:1, @(x) isfloat(x) && all(x >= 0));

            % Specify which file formats to save (setting to empty ('')
            % will not save anything)
            p.addParameter('save_format', {'png','fig','eps'}, @(x) any(ismember(x,{'png','fig','eps','pdf',''})));
            % Specify if save or not
            p.addParameter('save_flag', true, @(x) islogical(x) || x == 0 || x == 1);
            % Specify if actually plot or not (the plot in question varies from function to function)
            p.addParameter('do_plot', true, @(x) islogical(x) || x == 0 || x == 1);
            % Specify if the units of the plot should be denormalized
            p.addParameter('denormalize_flag', true, @(x) islogical(x) || x == 0 || x == 1);
            % (plot_field_density) Specifiy which property to plot, or if both
            p.addParameter('property_plot', 'both', @(x) any(ismember(x,{'wakefields','density','both','phasespace'})));
            % ...
            p.addParameter('include_phasespace', false, @(x) islogical(x) || x == 0 || x == 1);
            % (plot_field_density) Specify if movie should be created
            p.addParameter('create_movie', false, @(x) islogical(x) || x == 0 || x == 1);
            % (plot_field_density) Specify if struct with movie frames
            % should be saved
            p.addParameter('save_movie_struct', false, @(x) islogical(x) || x == 0 || x == 1);

            % handle of the figure to save
            p.addParameter('fig_handle',[], @(x) ishghandle(x,'figure'));

            % flag to put a title in the figures
            p.addParameter('title_flag', true, @(x) islogical(x) || x == 0 || x == 1);

            % flag to make a pause in a series of plots
            p.addParameter('make_pause', false, @(x) islogical(x) || x == 0 || x == 1);

            % plot scale
            p.addParameter('plot_scale', 'linear', @(x) any(ismember(x,{'linear','log'})));
            p.addParameter('plot_density_lims', [-inf inf], @(x) isfloat(x));
            p.addParameter('plot_field_lims', [-inf inf], @(x) isfloat(x));

            % figure number
            p.addParameter('fig_number', 1, @(x) isfloat(x) && x > 0);
            p.addParameter('exlegends', {'',''}, @(x) iscell(x));
            p.addParameter('extitles', {'',''}, @(x) iscell(x));

            % field density plot options
            p.addParameter('mirror_flag', true, @(x) islogical(x) || x == 0 || x == 1);
            p.addParameter('field_geometry', 'cartesian', @ (x) any(ismember(x,{'cartesian','cylindrical'})) || x == 0 || x == 1);
            

            p.KeepUnmatched = true;
            p.parse(varargin{:});
            if isempty(fieldnames(p.Unmatched))
                unmatched = {};
            else
                unmatched = reshape([fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj = obj@OsirisDenormalizer(unmatched{:}); %Parse to superclass OsirisDenormalizer.m

            obj.extradatadir          = p.Results.extradatadir;
            obj.plots_dir             = p.Results.plots_dir;
            obj.plot_name             = p.Results.plot_name;
            obj.dump_list             = p.Results.dump_list;
            obj.save_format           = p.Results.save_format;
            obj.save_flag             = p.Results.save_flag;
            obj.do_plot               = p.Results.do_plot;
            obj.fig_handle            = p.Results.fig_handle;
            obj.denormalize_flag      = p.Results.denormalize_flag;
            obj.property_plot         = p.Results.property_plot;
            obj.create_movie          = p.Results.create_movie;
            obj.save_movie_struct     = p.Results.save_movie_struct;
            obj.title_flag            = p.Results.title_flag;
            obj.make_pause            = p.Results.make_pause;
            obj.plot_scale            = p.Results.plot_scale;
            obj.plot_density_lims     = p.Results.plot_density_lims;
            obj.plot_field_lims       = p.Results.plot_field_lims;
            obj.fig_number            = p.Results.fig_number;
            obj.exlegends             = p.Results.exlegends;
            obj.extitles              = p.Results.extitles;
            obj.mirror_flag           = p.Results.mirror_flag;
            obj.field_geometry        = p.Results.field_geometry;
            obj.include_phasespace    = p.Results.include_phasespace;

            if obj.field_geometry == 1
                obj.field_geometry = 'cylindrical';
            elseif obj.field_geometry == 0
                obj.field_geometry = 'cartesian';
            end % field geometry 1 

        end % constructor

        function obj = save_plot(obj,varargin)

            if nargin > 1
                obj.fig_handle = varargin{1};
            end

            plotsdir = 'plots/';
            if obj.save_flag
                if ~isfolder([plotsdir,obj.plots_dir]); mkdir([plotsdir,obj.plots_dir]); end


                % save as png image
                if any(ismember(obj.save_format,{'png'}))
                    if ~isfolder([plotsdir,obj.plots_dir,'/png/']); mkdir([plotsdir,obj.plots_dir,'/png/']); end
                    exportgraphics(obj.fig_handle,[plotsdir,obj.plots_dir,'/png/',obj.plot_name,'.png'])
                end

                % save as matlab figure
                if any(ismember(obj.save_format,{'fig'}))
                    if ~isfolder([plotsdir,obj.plots_dir,'/fig/']); mkdir([plotsdir,obj.plots_dir,'/fig/']); end
                    savefig(obj.fig_handle,[plotsdir,obj.plots_dir,'/fig/',obj.plot_name,'.fig'])
                end

                % save as vector image
                if any(ismember(obj.save_format,{'eps'}))
                    if ~isfolder([plotsdir,obj.plots_dir,'/eps/']); mkdir([plotsdir,obj.plots_dir,'/eps/']); end
                    exportgraphics(obj.fig_handle,[plotsdir,obj.plots_dir,'/eps/',obj.plot_name,'.eps'],'ContentType','vector')
                end

                % save as vector image
                if any(ismember(obj.save_format,{'pdf'}))
                    if ~isfolder([plotsdir,obj.plots_dir,'/pdf/']); mkdir([plotsdir,obj.plots_dir,'/pdf/']); end
                    exportgraphics(obj.fig_handle,[plotsdir,obj.plots_dir,'/pdf/',obj.plot_name,'.pdf'],'ContentType','vector')
                end

            end % save flag
        end % save_plot

        function obj = waterfall_plot(obj)

            if obj.fig_number > 0
                fig_waterfall = figure(obj.fig_number);
                %                 fig_waterfall = figure('visible','off');
            else
                fig_waterfall = figure;
            end

            obj.fig_handle = fig_waterfall;
            obj.waterfall_handle = imagesc(obj.waterfall_xi,...
                obj.waterfall_z/545.39872233,rot90(obj.waterfall_mat,2));
            switch obj.property
                case 'fields'
                    colormap(gca,bluewhitered);
                case 'density'
                    c = gray;
                    c = flipud(c);
                    colormap(c);
            end
            drawnow;
            ax = gca;
            ax.XDir = 'reverse';
            ax.YDir = 'normal';
            cbar = colorbar;
            switch obj.property
                %                 case 'fields'
                %                     switch obj.wakefields_direction
                %                         case 'long'
                %                             cbar.Label.String = 'E_z (MV/m)';
                %                         case 'trans'
                %                             cbar.Label.String = 'W_r (MV/m)';
                %                         otherwise
                %                             cbar.Label.String = 'E_z (MV/m)';
                %                     end

                case 'density'
                    cbar.Label.String = 'charge (e)';
            end
            ylim([obj.waterfall_z(1) obj.waterfall_z(2)]/545.39872233)
            ylabel('$z$ (m)','Interpreter','latex');
            xlabel('$\xi$ (cm)','Interpreter','latex');


        end % waterfall_plot

        function obj = plot_field(obj,varargin)

            % parse input to load files
            p = inputParser;

            p.addParameter('field_plot', 0, @(x) isfloat(x));
            p.addParameter('r_plot', 0, @(x) isfloat(x));
            p.addParameter('z_plot', 0, @(x) isfloat(x));
            p.addParameter('tile_number', 0, @(x) isfloat(x));

            p.parse(varargin{:});

            field_plot         = p.Results.field_plot;
            r_plot             = p.Results.r_plot;
            z_plot             = p.Results.z_plot;
            tile_number        = p.Results.tile_number;

            if length(field_plot) <= 1
                [field_plot,~,r_plot,z_plot] = obj.load_data_plot();
            elseif all(r_plot == 0) || all(z_plot == 0)
                error('Please also give r_plot and z_plot.')
            end

            % mirroring for a better image, due to cylindrical symmetry

            if strcmp(obj.field_geometry,'cartesian')
                mirror_factor = -1;
            elseif strcmp(obj.field_geometry,'cylindrical')
                mirror_factor = 1;
            end

            field_plot_mirrored = [mirror_factor*flip(field_plot);field_plot];

            % Establish the maximum field value
            % the standard deviation is used as a measure the avoid
            % noisy peaks that sets a wrong scale for the opaqueness
            % get 3 times the std deviation with no weights
            meanstd_field = 5*std(abs(field_plot),[],'all')+eps;


            if isempty(obj.tile_handle)
                obj.tile_handle = tiledlayout(1,1);
            end

            % set the limits to the axes
            obj.ax_field = axes('parent',obj.tile_handle,'NextPlot','add');
            obj.ax_field.Layout.Tile = tile_number;

            obj.ax_field.XLim = [min(z_plot),max(z_plot)];
            obj.ax_field.YLim = r_plot;
            obj.ax_field.XDir = 'reverse';
            imagesc(obj.ax_field,'XData',z_plot,'YData',r_plot,'CData',field_plot_mirrored,[-meanstd_field meanstd_field]);
            c_field = colorbar('location','eastoutside');
            switch obj.wakefields_direction
                case 'trans'
                    colorbar_string = ['$W_r$ (MV/m)'];
                case 'long'
                    colorbar_string = ['$E_z$ (MV/m)'];
            end

            c_field.Label.String = colorbar_string;
            c_field.Label.Interpreter = 'latex';
            colormap(obj.ax_field,bluewhitered);
        end % plot field

        function obj = plot_density(obj,varargin)

            % parse input to load files
            p = inputParser;

            p.addParameter('density_plot', 0, @(x) isfloat(x));
            p.addParameter('r_plot', 0, @(x) isfloat(x));
            p.addParameter('z_plot', 0, @(x) isfloat(x));
            p.addParameter('tile_number', 0, @(x) isfloat(x));
            p.addParameter('plot_number', 1, @(x) isfloat(x));
            p.addParameter('grad_color_triplet', 0, @(x) isfloat(x));

            p.parse(varargin{:});

            density_plot       = p.Results.density_plot;
            r_plot             = p.Results.r_plot;
            z_plot             = p.Results.z_plot;
            tile_number        = p.Results.tile_number;
            plot_number        = p.Results.plot_number;
            grad_color_triplet = p.Results.grad_color_triplet;

%             switch obj.species
%                 case {'positrons','electrons'}
%                     density_plot = density_plot - median(density_plot,'all');
%             end

            switch obj.plot_scale
                case 'log'
                    density_plot_l = log(density_plot+1);
                case 'linear'
                    density_plot_l = abs(density_plot);
            end % switch plot scale
            %  clear density_plot

            if obj.mirror_flag
                density_plot_mirrored = [flip(density_plot_l);density_plot_l];
            else
                density_plot_mirrored = density_plot_l;
            end
            clear density_plot_l

            % stablish the opaque index
            % the standard deviation is used as a measure the avoid
            % noisy peaks that sets a wrong scale for the opaqueness
            % get 3 times the std deviation with no weights
            meanstd_density = 3*std(density_plot_mirrored,[],'all') + eps;
            max_opaqueness = 1;
            ind_opaque = max_opaqueness*density_plot_mirrored;
            ind_opaque(density_plot_mirrored > meanstd_density) = max_opaqueness*meanstd_density;
            ind_opaque = ind_opaque/max(ind_opaque,[],'all');

            if isempty(obj.tile_handle)
                obj.tile_handle = tiledlayout(1,1);
            end
            
            if plot_number > 1
                old_ax_density = obj.ax_density;
            end

            obj.ax_density = axes('parent',obj.tile_handle,'NextPlot','add','color','none');
            obj.ax_density.Layout.Tile = tile_number;
            % set limits to the axis
            obj.ax_density.XLim = [min(z_plot),max(z_plot)];
            obj.ax_density.YLim = [min(r_plot),max(r_plot)];
            obj.ax_density.XDir = 'reverse';


            switch obj.plot_scale
                case 'log'
                    imagesc(obj.ax_density,'XData',z_plot,'YData',r_plot,'CData',density_plot_mirrored);
                case 'linear'
                    imagesc(obj.ax_density,'XData',z_plot,'YData',r_plot,'CData',double(density_plot_mirrored>0),'alphadata',ind_opaque);
                    %                     grad = colorGradient([1 1 1],[0 0 0],2);
                    switch plot_number
                        case 1 
                            grad = [1 1 1; 0 0 0];
                        case 2 
                            grad = [221,110,15]/255;
                            linkaxes([old_ax_density,obj.ax_density],'xy')
                            old_ax_density.XTickLabel = [];
                            old_ax_density.YTickLabel = [];
                        case 3
                            grad = [8, 146, 208]/255;
                            linkaxes([old_ax_density,obj.ax_density],'xy')
                            old_ax_density.XTickLabel = [];
                            old_ax_density.YTickLabel = [];
                    end
                    
                    if any(grad_color_triplet ~= 0)
                        grad = grad_color_triplet;
                    end

                    colormap(obj.ax_density,grad);

            end % switch plot scale
            
        end % plot density

        function obj = plot_phasespace(obj,varargin)

            % parse input to load files
            p = inputParser;

            p.addParameter('phasespace_plot', 0, @(x) isfloat(x));
            p.addParameter('r_plot', 0, @(x) isfloat(x));
            p.addParameter('z_plot', 0, @(x) isfloat(x));
            p.addParameter('tile_number', 0, @(x) isfloat(x));

            p.parse(varargin{:});

            phasespace_plot    = p.Results.phasespace_plot;
            r_plot             = p.Results.r_plot;
            z_plot             = p.Results.z_plot;
            tile_number        = p.Results.tile_number;

            if length(phasespace_plot) <= 1
                [phasespace_plot,~,r_plot,z_plot] = obj.load_data_plot();
            elseif all(r_plot == 0) || all(z_plot == 0)
                error('Please also give r_plot and z_plot.')
            end

            % mirroring for a better image, due to cylindrical symmetry

            if strcmp(obj.field_geometry,'cartesian')
                mirror_factor = -1;
            elseif strcmp(obj.field_geometry,'cylindrical')
                mirror_factor = 1;
            end

            phasespace_plot_mirrored = [mirror_factor*flip(phasespace_plot);phasespace_plot];

            % Establish the maximum field value
            % the standard deviation is used as a measure the avoid
            % noisy peaks that sets a wrong scale for the opaqueness
            % get 3 times the std deviation with no weights
            meanstd_phasesapce = 5*std(abs(phasespace_plot_mirrored),[],'all')+eps;


            if isempty(obj.tile_handle)
                obj.tile_handle = tiledlayout(1,1);
            end

            % set the limits to the axes
            obj.ax_phasespace = axes('parent',obj.tile_handle,'NextPlot','add');
            obj.ax_phasespace.Layout.Tile = tile_number;

            obj.ax_phasespace.XLim = [min(z_plot),max(z_plot)];
            obj.ax_phasespace.YLim = r_plot;
            obj.ax_phasespace.XDir = 'reverse';
            imagesc(obj.ax_phasespace,'XData',z_plot,'YData',r_plot,'CData',...
                phasespace_plot_mirrored,[-meanstd_phasesapce meanstd_phasesapce]);
            c_phasespace = colorbar('location','eastoutside');
            %%%%%% HARDCODED
            obj.phasespace_direction = 'trans';
            %%%%%% HARDCODED
            colorbar_string = [obj.phasespace_direction,'. momentum ($\mathrm{m_b}$c)'];
            
            c_phasespace.Label.Interpreter = 'latex';
            c_phasespace.Label.String = colorbar_string;
            
            colormap(obj.ax_phasespace,bluewhitered);
        end % plot phasespace


        function obj = plot_field_density(obj,varargin)

            % parse input to load files
            p = inputParser;

            % directory where the plot should be saved
            p.addParameter('field_plot', 0, @(x) isfloat(x));
            p.addParameter('density_plot', 0, @(x) isfloat(x));
            p.addParameter('r_plot', 0, @(x) isfloat(x));
            p.addParameter('z_plot', 0, @(x) isfloat(x));
            p.addParameter('mirror_flag', true, @(x) islogical(x) || x == 0 || x == 1);
            p.addParameter('ylineflag', false, @(x) islogical(x) || x == 0 || x == 1);
            p.addParameter('xlineflag', false, @(x) islogical(x) || x == 0 || x == 1);
            p.addParameter('ylinepos', [0,0], @(x) isfloat(x));
            p.addParameter('xlinepos', [0,0], @(x) isfloat(x));
            p.addParameter('include_lineout', 'no', @(x) any(ismember(x,{'no',...
                'both','field_lineout','density_profile','phasespace_lineout'})));
            p.addParameter('include_field_lineout', false, @(x) islogical(x) || x == 0 || x == 1);
            p.addParameter('include_density_profile', false, @(x) islogical(x) || x == 0 || x == 1);
            p.addParameter('include_phasespace_profile', false, @(x) islogical(x) || x == 0 || x == 1);

            p.parse(varargin{:});

            field_plot              = p.Results.field_plot;
            density_plot            = p.Results.density_plot;
            r_plot                  = p.Results.r_plot;
            z_plot                  = p.Results.z_plot;

            obj.mirror_flag         = p.Results.mirror_flag;
            obj.ylineflag           = p.Results.ylineflag;
            obj.xlineflag           = p.Results.xlineflag;
            obj.ylinepos            = p.Results.ylinepos;
            obj.xlinepos            = p.Results.xlinepos;
            obj.include_lineout     = p.Results.include_lineout;
            obj.include_field_lineout = ...
                                      p.Results.include_field_lineout;
            obj.include_density_profile = ...
                                      p.Results.include_density_profile;
            obj.include_phasespace_profile = ...
                                      p.Results.include_phasespace_profile;

            fontsize_all = 20;


            if obj.create_movie
                [~,video] = obj.setup_movie();
            end % create movie

            for n = 1:length(obj.dump_list)

                obj.dump = obj.dump_list(n);

                i_secondimagesc = 1;

                if (~strcmp(obj.extradatadir,'')) || obj.include_phasespace
                    i_secondimagesc = 2;
                end

                if obj.fig_number > 0
                    fig_double = figure(obj.fig_number);
                else
                    fig_double = figure;
                end

                % if there is an extradatadir, it brings another fielden
                % plot, and another phasespace plot (if required, barely
                % used), maximum possible should be 6 tiles
                total_number_of_plots = ...
                    (~strcmp(obj.extradatadir,'') + 1)*(1 + obj.include_phasespace)...
                    + obj.include_field_lineout + obj.include_density_profile;

                obj.tile_handle = tiledlayout(fig_double,total_number_of_plots,1);

                switch total_number_of_plots
                    case 1
                        plot_position = [0.1 0.1 0.8 0.5];
                    case 2
                        plot_position = [0.0471 0.2560 0.9081 0.5111];
                    case 3 
                        plot_position = [0.0346 0.1231 0.9008 0.7356];
                    case 4 
                        plot_position = [0.0471 0.0560 0.9081 0.94];
                    case 5
                        plot_position = [0.0471 0.0560 0.9081 0.94]; % not tested yet
                    case 6
                        plot_position = [0.0471 0.0560 0.9081 0.94]; % not tested yet
                end

                obj.tile_handle.TileSpacing = 'compact';
                obj.tile_handle.Padding = 'compact';

                if n == 1
                    fig_double.Units = 'normalized';
                    fig_double.OuterPosition = plot_position;
                end

                species_to_plot = {'electron_seed','antiproton_beam','electrons'};

                % cell initialization 
                field_plot_save = cell(2,1);
                density_plot_save = cell(2,1);
                density_plot_e_save = cell(2,1);
                r_plot_save = cell(2,1);
                z_plot_save = cell(2,1);
                z_lineout_den_save = cell(2,1);
                z_lineout_fld_save = cell(2,1);
                %density_plot = cell(length(species_to_plot),1);

                

                for i_ex = 1:i_secondimagesc %-------- extra datadir loop

                    if i_ex == 2 && (~obj.include_phasespace)%%-----------------new ex
                        field_plot_save{1} = field_plot;
                        density_plot_save{1} = density_plot;
                        if isfile(e_fullpath)
                            density_plot_e_save{1} = density_plot_e;
                        end
                        
                        r_plot_save{1} = r_plot;
                        z_plot_save{1} = z_plot;
                        z_lineout_den_save{1} = z_lineout_den;
                        z_lineout_fld_save{1} = z_lineout_fld;

                        field_plot = 0;
                        density_plot = 0;
                        datadir_save = obj.datadir;
                        obj.datadir = obj.extradatadir;
                    end %%--------------------------new ex

                    if i_ex == 2 && obj.include_phasespace
                        save_property_plot = obj.property_plot;
                        obj.property_plot = 'phasespace';
                    end

                    load_data_switch = ((numel(field_plot) == 1) && ...
                        (numel(density_plot) == 1)) && ...
                        (field_plot == 0 && density_plot == 0);

                    if load_data_switch
                        [field_plot,~,phasespace_plot,r_plot,z_plot,z_lineout_den,z_lineout_fld] = obj.load_data_plot();
                        obj.mirror_flag = true;
                    end

                    if i_ex == 2 && (~obj.include_phasespace)%%-----------------new ex
                        obj.ax_field.XTickLabel = [];
                        obj.ax_density.XTickLabel = [];
                        field_plot_save{2} = field_plot;
                        density_plot_save{2} = density_plot;
                        r_plot_save{2} = r_plot;
                        z_plot_save{2} = z_plot;
                        z_lineout_den_save{2} = z_lineout_den;
                        z_lineout_fld_save{2} = z_lineout_fld;
                        
                    end %%--------------------------new ex

                    include_quivers = 0;
                    if include_quivers
                        r_plot = r_plot/10;
                    end

                    if ismember(obj.property_plot,{'wakefields','both'})
                        obj.plot_field('field_plot',field_plot,'r_plot',r_plot,'z_plot',z_plot,'tile_number',i_ex);
                        obj.ax_field.FontSize = fontsize_all;
                    end

                    if ismember(obj.property_plot,{'density','both'})
                        %%-------------------------------------------------------------------
                        %%-------------------------------------------------------------------
                        %%-------------------------------------------------------------------
                        %%-------------------------------------------------------------------
                        %%-------------------------------------------------------------------
                        %%-------------------------------------------------------------------
                        save_property_plot = obj.property_plot;
                        obj.property_plot = 'density'; % later change load data plot to only load field or density
                        for i_density = 1:length(species_to_plot)
                            obj.species = species_to_plot{i_density};
                            [~,density_plot_temp,~,r_plot,z_plot,~,~] = obj.load_data_plot();
                            density_plot_cell{i_density} = density_plot_temp;
                            clear density_plot_temp;
                            obj.plot_density('density_plot',density_plot_cell{i_density},...
                                'r_plot',r_plot,'z_plot',z_plot,'tile_number',i_ex,'plot_number',i_density);
                        end
                        obj.property_plot = save_property_plot;


%                         obj.plot_density('density_plot',density_plot,'r_plot',r_plot,'z_plot',z_plot,'tile_number',i_ex);
%                         obj.ax_density.FontSize = fontsize_all;
%                         % here1
%                         obj.property = 'density';
%                         obj.species = 'electron_seed';
%                         obj.justPath = 1;
%                         obj.getdata();
%                         obj.justPath = 0;
%                         e_fullpath = which(obj.fullpath);
%                         if isfile(e_fullpath)
%                             obj.a = 1;
%                             obj.getdata(); obj.assign_density(); 
%                             density_plot_e = obj.denorm_density(obj.nelectron_seed);
%                             obj.plot_density('density_plot',density_plot_e,'r_plot',r_plot,'z_plot',z_plot,'tile_number',i_ex);
%                             obj.a = 0;
%                         end
% %                         obj.species = 'proton_beam';
%                        obj.species = 'electrons';
                        obj.ax_density.FontSize = fontsize_all;
                    end

                    if i_ex == 2 && (~obj.include_phasespace) && isfile(e_fullpath)%%-----------------new ex
                        density_plot_e_save{2} = density_plot_e;
                    end %%--------------------------new ex

                    if strcmp(obj.property_plot,'both')
                        obj.ax_field.Position = obj.ax_density.Position;
                        linkaxes([obj.ax_field,obj.ax_density],'xy')
                    end

                    if ismember(obj.property_plot,{'phasespace'})
                        obj.plot_phasespace('phasespace_plot',phasespace_plot,'r_plot',r_plot,'z_plot',z_plot,'tile_number',i_ex);
                        obj.ax_phasespace.FontSize = fontsize_all;
                    end

                    if i_ex == 2 && (~obj.include_phasespace)
                       
                        obj.datadir = datadir_save;
                        field_plot = field_plot_save{1};
                        density_plot = density_plot_save{1};
                        if isfile(which(e_fullpath))
                            density_plot_e = density_plot_e_save{1};
                        end
                        r_plot = r_plot_save{1};
                        z_plot = z_plot_save{1};
                        z_lineout_den = z_lineout_den_save{1};
                        z_lineout_fld = z_lineout_fld_save{1};
                    end


                    if include_quivers
                        property_save = obj.property;
                        formatsave = obj.dataformat;
                        obj.dataformat = 'mat';
                        obj.property = 'raw';
                        obj.direction = 'r';

                        obj.raw_dataset = 'x'; obj.getdata(); obj.assign_raw();
                        obj.raw_dataset = 'p'; obj.getdata(); obj.assign_raw();

                        obj.direction = 'z';
                        obj.raw_dataset = 'x'; obj.getdata(); obj.assign_raw();

                        %                     ixi = ((obj.dtime + obj.simulation_window - obj.denorm_distance(obj.n_simulation_window+obj.nz_raw)) < obj.xi_range(1)) &...
                        %                         (((obj.dtime + obj.simulation_window - obj.denorm_distance(obj.n_simulation_window+obj.nz_raw)) > obj.xi_range(2)));
                        ixi = (abs(obj.denorm_distance(obj.nz_raw)) < obj.xi_range(1)) &...
                            (abs(obj.denorm_distance(obj.nz_raw)) > obj.xi_range(2));


                        %                     xi_temp = obj.dtime + obj.simulation_window - obj.denorm_distance(obj.nz_raw(ixi));
                        xi_temp = abs(obj.denorm_distance(obj.nz_raw(ixi)));
                        r_temp = obj.denorm_distance(obj.nr_raw(ixi));
                        pr_temp = obj.npr_raw(ixi);

                        k = 1000;
                        is = randperm(length(pr_temp),k);

                        xi = xi_temp(is);
                        r = r_temp(is);
                        pr = pr_temp(is);


                        ip = pr > 0;

                        xi_p = xi(ip);
                        r_p = r(ip);
                        pr_p = pr(ip);

                        xi_n = xi(~ip);
                        r_n = r(~ip);
                        pr_n = pr(~ip);



                        scale = 0.5;
                        hold on
                        qp1 = quiver(xi_p,r_p,zeros(length(pr_p),1),pr_p,scale,'color',[0.5391,0.1680,0.8828],...
                            'ShowArrowHead',1,'LineWidth',1);
                        qp2 = quiver(xi_p,-r_p,zeros(length(pr_p),1),-pr_p,scale,'color',[0.5391,0.1680,0.8828],...
                            'ShowArrowHead',1,'LineWidth',1);
                        qn1 = quiver(xi_n,r_n,zeros(length(pr_n),1),pr_n,scale,'color',[0.1328,0.5430,0.1328],...
                            'ShowArrowHead',1,'LineWidth',1);
                        qn2 = quiver(xi_n,-r_n,zeros(length(pr_n),1),-pr_n,scale,'color',[0.1328,0.5430,0.1328],...
                            'ShowArrowHead',1,'LineWidth',1);
                        hold off

                        obj.property = property_save;
                        obj.dataformat = formatsave;
                    end


                    if obj.ylineflag
                        hold on
                        yline(obj.ylinepos*10,'r','LineWidth',1)
                        hold off
                    end

                    if obj.xlineflag
                        hold on
                        xline(obj.xlinepos,'color',[0.99,0.27,0],'LineWidth',1)
                        hold off
                    end

                    ylabel('x (mm)','Interpreter','Latex');
                    %%%____________________________________________________
                    %                 ylim([-2.1,2.1]);

                    if i_secondimagesc == 2 && (~obj.include_phasespace)
                       % title(obj.extitles{i_ex},'Interpreter','Latex');
                    end

                   
                    if i_ex == 2 && obj.include_phasespace
                        obj.ax_phasespace.XTickLabel = [];
                        obj.property_plot = save_property_plot;
                    end

                end % --------------- for i_ex


                if obj.title_flag
                    title(obj.tile_handle,['z = ',num2str(obj.propagation_distance/100,3),' m',''],'Interpreter','Latex','FontSize',fontsize_all)
                end

                if include_quivers %return values of r_plot if quivers were used
                    r_plot = r_plot*10;
                end

                 any_lineout_flag = obj.include_field_lineout | ...
                        obj.include_density_profile | obj.include_phasespace_profile;

                    if any_lineout_flag
                        obj.ax_density.XTickLabel = [];
                        obj.ax_field.XTickLabel = [];
                        obj.ax_phasespace.XTickLabel = [];
                    end

                    if ismember(obj.include_lineout,{'density_profile','field_lineout','both'})
                       obj.ax_density.XTickLabel = [];
                        
                       obj.ax_field.XTickLabel = [];
                    end


                if obj.include_density_profile

                    r_lineplot = linspace(0,max(r_plot),size(density_plot,1));
                    ir = (r_lineplot < obj.ylinepos(2)*10) & (r_lineplot > obj.ylinepos(1)*10);

                    int_option = 'sum'; % just sum = density
                    % long_profile = obj.cylindrical_radial_integration(r_lineplot(ir)/10,density_plot(ir,:),int_option); %just sum
                    %                     long_profile = density_plot(obj.lineout_point,:);
                    % HARD CODED
                    long_profile = density_plot(10,:);
                    
                    %long_profile = smooth(long_profile,obj.plasma_wavelength/13.62/8,'loess');

                    ax_longprofile = nexttile;

                    plongprofile = plot(z_lineout_den,long_profile,'k');

                    if isfile(e_fullpath)
                        hold on
                        long_profile_e = obj.cylindrical_radial_integration(r_lineplot(ir)/10,abs(density_plot_e(ir,:)),int_option); %just sum
                        % HARD CODED
                        long_profile_e = density_plot_e(10,:);
                        %long_profile_e = smooth(long_profile_e,obj.plasma_wavelength/13.62/8,'loess');
                        if length(z_lineout_den) == length(long_profile_e)
                            % nothing
                        else
                            z_lineout_den = [0,z_lineout_den];
                        end
                        plongprofile_e = plot(z_lineout_den,long_profile_e,'color',[221,110,15]/255);
                        hold off
                        %legend({'proton bunch','electron bunch'},'interpreter','latex');
                        % here2
                    end


                    if i_secondimagesc == 2 && (~obj.include_phasespace) %----------------------------- new ex

                        density_plot2 = density_plot_save{2};
                        density_plot_e2 = density_plot_e_save{2};
                        z_lineout2 = z_lineout_den_save{2};
                        r_plot2 = r_plot_save{2};

                        r_lineplot2 = linspace(0,max(r_plot2),size(density_plot2,1));
                        ir2 = (r_lineplot2 < obj.ylinepos(2)*10) & (r_lineplot2 > obj.ylinepos(1)*10);

                        long_profile2 = obj.cylindrical_radial_integration(r_lineplot2(ir2)/10,density_plot2(ir2,:),'sum'); %just sum
                        long_profile2 = smooth(long_profile2,obj.plasma_wavelength/13.62/8,'loess');
                        hold on
                        plongprofile2 = plot(z_lineout2,long_profile2,...
                            'LineStyle','-','Color',[255 36 0]/256);
                        hold off

                        if isfile(which(obj.fullpath))
                        hold on
                        long_profile_e2 = obj.cylindrical_radial_integration(r_lineplot(ir)/10,abs(density_plot_e2(ir,:)),int_option); %just sum
                        long_profile_e2 = smooth(long_profile_e2,obj.plasma_wavelength/13.62/8,'loess');
                        plongprofile_e2 = plot(z_lineout_den,long_profile_e2,'color',[221,110,15]/255);
                        hold off
                        %legend({'proton bunch','electron bunch'},'interpreter','latex');
                        % here2
                        end
%                    obj.species = 'proton_beam';
                        obj.species = 'electrons';
                    obj.justPath = 0;

                        leg_handle = legend([plongprofile,plongprofile2],obj.extitles,'Interpreter','Latex');
                        leg_handle.AutoUpdate = 'off';
                    end %----------------------------- new ex

                    ax_longprofile.XDir = 'reverse';
                    ax_longprofile.YTick = [];

                    if strcmp(obj.include_lineout,'density_profile')
                        %  ax_longprofile.XTick = [z_lineplot];
                        %  ax_longprofile.XTickLabel = [z_lineplot];
                        xlabel('$\xi$ (cm)','Interpreter','Latex');
                    end

                    if obj.xlineflag
                        hold on
                        xline(obj.xlinepos,'color', [0.99,0.27,0], 'LineWidth', 1)
                        hold off
                    end

                    xlim(obj.ax_density.XLim)
                    ylim(obj.plot_density_lims);
                    if ~strcmp(int_option,'just_sum')
                        ylabel({'charge','(arb. units)'},'Interpreter','Latex');
                    else
                        ylabel({'density','(arb. units)'},'Interpreter','Latex');
                    end
                    ax_longprofile.FontSize = fontsize_all;

                end % if include density profile

                if obj.include_phasespace_profile
                    ax_longprofile.XTickLabel = [];
                    obj.load_lineout = 1;
                    obj.property = 'phasespace';
                    obj.getdata();
                    ax_phasespace_lineout = nexttile;
                    x = load("colororder_defaultblack.mat"); %corder_default
                    %corder_defblack = x.corder_defblack;
                    corder_defblack = {'r','b'};
                    linestyles = {'-',':','--','-.','-','-','-','-.','--',':'};

                    obj.title_flag = 0;

                    obj.plot_lineout2('lineout_plot',obj.p2lineout_p,'z_plot',obj.denorm_distance(obj.nz));
                    obj.plot_handle.LineStyle = linestyles{1};
                    obj.plot_handle.Color = corder_defblack{1};
                    obj.plot_handle.LineWidth = 1;
                    
                    hold on
                    obj.plot_lineout2('lineout_plot',abs(obj.p2lineout_n),'z_plot',obj.denorm_distance(obj.nz));
                    obj.plot_handle.LineStyle = linestyles{3};
                    obj.plot_handle.Color = corder_defblack{2};
                    obj.plot_handle.LineWidth = 1;

                    hold off
                    legs = {'$p+$ (away from axis)','$p-$ (towards axis)'};
                    legend(legs{:},'location','best','interpreter','latex',...
                        'FontSize',obj.plot_fontsize-6);
                    ylabel({'trans.','momen. ($\mathrm{m_b c}$)'},'fontsize',fontsize_all,'interpreter','latex')

                    ax_phasespace_lineout.XDir = 'reverse';
                    ax_phasespace_lineout.FontSize = fontsize_all;
                    ax_phasespace_lineout.YTickLabel = [];
                    xlim(obj.ax_field.XLim);

                    obj.title_flag = 1;
                    obj.load_lineout = 0;
                end
                
                include_alpha = 0;
                if include_alpha

                end % incldue alpha

                %if ismember(obj.include_lineout,{'field_lineout','both'})
                if obj.include_field_lineout

                    ax_longprofile.XTickLabel = [];
                    ax_phasespace_lineout.XTickLabel = [];
                    ax_phasespace_lineout.XLabel = [];

                    ax_lineout = nexttile;
                    if ischar(obj.lineout_point)
                        lineout_point_char = obj.lineout_point;
                        r_lineout1 = linspace(0,max(r_plot),size(field_plot,1))/10;
                        [~,minloc] = min(abs(r_lineout1 - str2double(obj.lineout_point)));
                        obj.lineout_point = minloc;
                    end

                    obj.lineout = field_plot(obj.lineout_point,:);

                    plineout = plot(z_lineout_fld,obj.lineout,'color',[0 0.4470 0.7410]);
                    yline(0);
                    % 0.929 0.694 0.125

                    if i_secondimagesc == 2 && (~obj.include_phasespace) %----------------------------- new ex
                        obj.lineout_point = lineout_point_char;

                        field_plot2 = field_plot_save{2};
                        r_plot2 = r_plot_save{2};

                        if ischar(obj.lineout_point)
                            r_lineout2 = linspace(0,max(r_plot2),size(field_plot2,1))/10;
                            [~,minloc] = min(abs(r_lineout2 - str2double(obj.lineout_point)));
                            obj.lineout_point = minloc;
                        end


                        lineout2 = field_plot2(obj.lineout_point,:);

                        hold on
                        plineout2 = plot(z_lineout_fld,lineout2,...
                            'LineStyle','-','Color',[0.929 0.694 0.125]);

                        hold off

                        leg_handle = legend([plineout,plineout2],obj.exlegends,'Interpreter','Latex');
                        leg_handle.AutoUpdate = 'off';

                        obj.lineout_point = lineout_point_char;

                    end %----------------------------- new ex

                    ax_lineout.XDir = 'reverse';
                    ax_lineout.FontSize = fontsize_all;

                    xlim(obj.ax_field.XLim);

                    if strcmp(obj.include_lineout,'both')

                        %                         ax_longprofile.XTickLabel = [];
                        %                         ax_longprofile.XTick = [z_lineplot];

                    end
                    ylim(obj.plot_field_lims);

                    switch obj.wakefields_direction
                        case 'long'
                            ylabel({'$E_z$ (MV/m)'},'Interpreter','Latex','fontsize',fontsize_all);
                        case 'trans'
                            ylabel({'$W_r$ (MV/m)'},'Interpreter','Latex','fontsize',fontsize_all);

                    end

                    xlabel('$\xi$ (cm)','Interpreter','Latex');

                end % if include long profile

                drawnow;

                if obj.make_pause && (n < length(obj.dump_list))
                    pause(1);
                end

                obj.fig_handle = fig_double;

                %                 if strcmp(obj.plot_name,'plot')
                obj.plot_name = [obj.datadir,obj.property_plot,'n',num2str(obj.dump),...
                    'xi',num2str(round(obj.xi_range(1))),'xi',...
                    num2str(round(obj.xi_range(2))),'t',num2str(round(obj.trans_range(1))),...
                    't',num2str(round(obj.trans_range(2)))];
                %                 else
                %                     obj.plot_name = [obj.datadir,'n',num2str(obj.dump)];
                %                 end % if plot name

                obj.save_plot();


                if obj.create_movie
                    obj.struct_movie(obj.frame_counter+1) = getframe(gcf);
                    obj.frame_counter = obj.frame_counter + 1;
                end % create movie

                if n < length(obj.dump_list)
                    clf
                    density_plot = 0;
                    field_plot = 0;
                    r_plot = 0;
                    z_plot = 0;
                end % clear figure
                obj.progress_dump('Plotting 2D',n,length(obj.dump_list))
            end % length dump list

            if obj.create_movie
                writeVideo(video,obj.struct_movie);
                % uncomment when the structure needs saving to optimize
                % videos (15/03/2021)
                %                 struct_movie_save = obj.struct_movie;
                %                 save(['save_files/field_density/struct',obj.datadir,...
                %                     '_',obj.wakefields_direction,'.mat'],'struct_movie_save')
                close(video);
            end % if create movie

        end % plot_field_density

        function [obj,video] = setup_movie(obj)

            %             movie_dir = ['movies/field_density',obj.include_lineout,'/'];
            movie_dir = ['movies/r2a_collab','/'];

            struct_dir = ['save_files/field_density',obj.include_lineout,'/'];
            if ~isfolder(movie_dir)
                mkdir(movie_dir);
            end
            if ~isfolder(struct_dir)
                mkdir(struct_dir);
            end
            video = VideoWriter([movie_dir,...
                obj.wakefields_direction,obj.datadir,obj.property_plot,...
                'xi',num2str(round(obj.xi_range(1))),'xi',...
                num2str(round(obj.xi_range(2))),'t',num2str(round(obj.trans_range(1))),...
                't',num2str(round(obj.trans_range(2))),'x.avi']);
            video.FrameRate = 4;
            open(video);
            struct_movie_temp(length(obj.dump_list)) = struct('cdata',[],'colormap',[]);
            obj.struct_movie = struct_movie_temp;
        end % setup movie

        function [field_plot,density_plot,phasespace_plot,r_plot,z_plot,z_lineout_den,z_lineout_fld] = load_data_plot(obj)
            %initialize variables to make life simpler afterwards (less
            %switches and ifs)
            field_plot = 0;
            density_plot = 0;
            phasespace_plot = 0;
            z_lineout_fld = 0;
            z_lineout_den = 0;

            if obj.include_phasespace
                obj.property = 'phasespace';
                obj.getdata(); obj.trim_data();
                phasespace_plot = obj.ndataOut;
            end

            switch obj.property_plot
                case 'wakefields'
                    obj.property = 'fields';
                    obj.getdata();
                    if obj.denormalize_flag
                        obj.trim_data();
                        obj.denorm_Efield();
                        field_plot = obj.([obj.wakefields_direction,'field']);
                    else
                        field_plot = obj.ndataOut;
                    end
                    z_lineout_fld = obj.dtime + obj.simulation_window - obj.z;
                    z_lineout = z_lineout_fld;

                case 'density'
                    obj.property = 'density';
                    obj.justPath = 1;
                    obj.getdata();
                    obj.justPath = 0;
                    species_fullpath = which(obj.fullpath);
                    if isfile(species_fullpath)
                        obj.justPath = 0;
                        obj.getdata();
                        if obj.denormalize_flag
                            obj.trim_data();
                            obj.denorm_density();
                            density_plot = obj.(obj.species);
                        else
                            density_plot = obj.ndataOut;
                        end % denorm flag
                        z_lineout_den = obj.dtime + obj.simulation_window - obj.z;
                        z_lineout = z_lineout_den;
                    else
                        warning(['No ',obj.species,' file found, skipping for plot.'])
                    end 

                case 'both'
                    obj.property = 'density';
                    obj.getdata();

                    if obj.denormalize_flag
                        obj.trim_data();
                    else
                        density_plot = obj.ndataOut;
                    end
                    z_lineout_den = obj.dtime + obj.simulation_window - obj.z;
                    

                    obj.property = 'fields';
                    obj.getdata();
                    if obj.denormalize_flag
                        obj.trim_data();
                    else
                        field_plot = obj.ndataOut;
                    end
                    z_lineout_fld = obj.dtime + obj.simulation_window - obj.z;
                    z_lineout = z_lineout_fld;

                    if obj.denormalize_flag
                        obj.denorm_Efield();
                        obj.denorm_density(); % look here
                        field_plot = obj.([obj.wakefields_direction,'field']);
                        % field_plot = obj.(['n',obj.wakefields_direction,'field']);

                        density_plot = obj.(obj.species);
                    end % denorm flag
            end % switch property_plot

            r_plot = [-max(obj.r),max(obj.r)]*10; % in mm
            % z_plot = ([min(obj.z),max(obj.z)]-obj.simulation_window)/100; % in m
            z_lineout_fld = obj.dtime + obj.simulation_window - obj.z;

            z_plot = [z_lineout(1),z_lineout(end)];

        end % load data field density plot

        % ------------------------------------------------------------------------------
        function obj = plot_lineout(obj)

            fontsize_all = 20;

            if obj.create_movie
                obj.frame_counter = 0;
                [~,video] = obj.setup_movie();
            end % create movie

            for n = 1:length(obj.dump_list)
                %initialize variables to make life simpler afterwards (less
                %switches and ifs)
                lineout_plot = 0;
                obj.dump = obj.dump_list(n);

                switch obj.property_plot
                    case 'wakefields'
                        obj.property = 'fields';
                        obj.getlineout();
                        if obj.denormalize_flag
                            obj.denorm_Efield();
                            lineout_plot = obj.lineout;
                        else
                            lineout_plot = obj.nlineout;
                        end
                    case 'density'
                        obj.property = 'density';
                        obj.getdata();
                        if obj.denormalize_flag
                            obj.trim_data();
                            obj.denorm_density();
                            density_plot = obj.(obj.species);
                        else
                            density_plot = obj.ndataOut;
                        end % denorm flag
                        lineout_plot = obj.cylindrical_radial_integration(obj.r,density_plot,'sum');
                        lineout_plot = smooth(lineout_plot,obj.plasma_wavelength/13.62,'loess');

                end % switch property_plot

                % select the axis

                if obj.denormalize_flag
                    z_plot = (obj.dtime + obj.simulation_window - obj.z); % cm
                    r_plot = obj.r;
                else
                    z_plot = (obj.ntime + obj.n_simulation_window - obj.nz)/(2*pi);
                    r_plot = obj.nr;
                end

                switch obj.lineout_direction
                    case 'trans'
                        zr_plot = r_plot;
                    case 'long'
                        zr_plot = z_plot;
                end

                % begin the plot

                fig_lineout = figure(obj.fig_number);
                if n == 1
                    fig_lineout.Units = 'normalized';
                    fig_lineout.OuterPosition = [0.1 0.3 0.8 0.5]; %[0.1 0.3 0.8 0.5]
                end

                obj.plot_handle = plot(zr_plot,lineout_plot,'k',...
                    'LineWidth',1);
                obj.lineout = lineout_plot;
                obj.z_out = z_plot;

                xlim(([min(zr_plot),max(zr_plot)]));
                ax_lo = gca;
                ax_lo.FontSize = obj.plot_fontsize;
                if strcmp(obj.lineout_direction,'long')
                    ax_lo.XDir = 'reverse';
                end

                obj.species = 'electron_seed';
                obj.justPath = 1;
                obj.getdata();
                obj.justPath = 0;
                e_fullpath = which(obj.fullpath);
                if isfile(e_fullpath) && strcmp(obj.property,'density')
                    obj.a = 1;
                    obj.getdata(); %obj.assign_density(); 
                    obj.trim_data();
                    density_plot_e = obj.denorm_density(obj.nelectron_seed);
                    lineout_plot = obj.cylindrical_radial_integration(obj.r,density_plot_e,'sum');
                    lineout_plot = smooth(lineout_plot,obj.plasma_wavelength/13.62,'loess');
                    hold on
                    plot(zr_plot,lineout_plot,'k','LineWidth',1);
                    hold off

                    obj.a = 0;
                end
               % obj.species = 'proton_beam'; 
                obj.species = 'electrons';
 
                if obj.denormalize_flag
                    if obj.title_flag
                        title([obj.datadir,'z = ',num2str(obj.propagation_distance/100,2),' m',''],'Interpreter','Latex','FontSize',fontsize_all)
                    end
                    ylabel([obj.wakefields_direction,'. fields (MV/m)'], 'FontSize', obj.plot_fontsize,'Interpreter','Latex')
                    xlabel('\xi (cm)', 'FontSize', obj.plot_fontsize,'Interpreter','Latex');
                else
                    if obj.title_flag
                        title(['z = ',num2str(obj.n_propagation_distance,2),''], 'FontSize', obj.plot_fontsize,'Interpreter','Latex')
                    end
                    ylabel([obj.wakefields_direction,'. fields'],'Interpreter','Latex')
                    xlabel('$\xi (\lambda_p)$','Interpreter','Latex');

                end

                switch obj.property_plot
                    case 'wakefields'
                        ylim(obj.plot_field_lims);
                    case 'density'
                        ylim(obj.plot_density_lims);
                end

                drawnow;

                if obj.make_pause && (n < length(obj.dump_list))
                    pause;
                end

                obj.fig_handle = fig_lineout;

                %                 if strcmp(obj.plot_name,'plot')
                obj.plot_name = [obj.datadir,obj.property_plot,'n',num2str(obj.dump),...
                    'xi',num2str(round(obj.xi_range(1))),'xi',...
                    num2str(round(obj.xi_range(2))),'t',num2str(round(obj.trans_range(1))),...
                    't',num2str(round(obj.trans_range(2)))];
                %                 else
                %                     obj.plot_name = [obj.plot_name,'n',num2str(obj.dump)];
                %                 end % if plot name

                obj.save_plot();

                if obj.create_movie
                    obj.struct_movie(obj.frame_counter+1) = getframe(gcf);
                    obj.frame_counter = obj.frame_counter + 1;
                end % create movie

                if n < length(obj.dump_list)
                    clf
                end % clear figure
                obj.progress_dump(['Plotting lineout ',obj.datadir],n,length(obj.dump_list))
            end % length dump list

            if obj.create_movie
                writeVideo(video,obj.struct_movie);
                struct_movie_save = obj.struct_movie;
                save(['save_files/lineout/struct',obj.datadir,...
                    '_',obj.wakefields_direction,'.mat'],'struct_movie_save')
                close(video);
            end % if create movie

        end % lineout plot

        function obj = plot_lineout2(obj,varargin)

            % parse input to load files
            p = inputParser;

            % directory where the plot should be saved
            p.addParameter('lineout_plot', 0, @(x) isfloat(x));
            p.addParameter('r_plot', 0, @(x) isfloat(x));
            p.addParameter('z_plot', 0, @(x) isfloat(x));

            p.parse(varargin{:});

            lineout_plot       = p.Results.lineout_plot;
            r_plot             = p.Results.r_plot;
            z_plot             = p.Results.z_plot;


            switch obj.lineout_direction
                case 'trans'
                    zr_plot = r_plot;
                case 'long'
                    zr_plot = z_plot;
            end

            % begin the plot

            fig_lineout = figure(obj.fig_number);
            fig_lineout.Units = 'normalized';
            %fig_lineout.OuterPosition = [0.1 0.3 0.8 0.5]; %[0.1 0.3 0.8 0.5]

            obj.plot_handle = plot(zr_plot,lineout_plot,'k',...
                'LineWidth',1);

            obj.lineout = lineout_plot;
            obj.z_out = z_plot;

            xlim(([min(zr_plot),max(zr_plot)]));
            ax_lo = gca;
            ax_lo.FontSize = obj.plot_fontsize;
            if strcmp(obj.lineout_direction,'long')
                ax_lo.XDir = 'reverse';
            end

            if obj.denormalize_flag
                if obj.title_flag
                    title(['z = ',num2str(obj.propagation_distance/100,2),' m',''],'Interpreter','Latex')
                end
                switch obj.property_plot
                    case 'fields'
                        ylabel([obj.wakefields_direction,'. fields (MV/m)'], 'FontSize', obj.plot_fontsize,'Interpreter','Latex')
                    case 'density'
                        ylabel(['density (arb. units)',''], 'FontSize', obj.plot_fontsize,'Interpreter','Latex')
                    case 'phasespace'
                        ylabel(['momentum (mc)',''], 'FontSize', obj.plot_fontsize,'Interpreter','Latex')
                end % switch property plot
                xlabel('$\xi$ (cm)', 'FontSize', obj.plot_fontsize,'Interpreter','Latex');
            else
                if obj.title_flag
                    title(['z = ',num2str(obj.n_propagation_distance,2),''], 'FontSize', obj.plot_fontsize,'Interpreter','Latex')
                end
                ylabel([obj.wakefields_direction,'. fields'],'Interpreter','Latex')
                xlabel('$\xi (\lambda_p)$','Interpreter','Latex');
            end

            drawnow;

            obj.fig_handle = fig_lineout;

            %                 if strcmp(obj.plot_name,'plot')
            obj.plot_name = [obj.datadir,obj.property_plot,'n',num2str(obj.dump),...
                'xi',num2str(round(obj.xi_range(1))),'xi',...
                num2str(round(obj.xi_range(2))),'t',num2str(round(obj.trans_range(1))),...
                't',num2str(round(obj.trans_range(2)))];
            %                 else
            %                     obj.plot_name = [obj.plot_name,'n',num2str(obj.dump)];
            %                 end % if plot name

            obj.save_plot();

%             if obj.create_movie
%                 obj.struct_movie(obj.frame_counter+1) = getframe(gcf);
%                 obj.frame_counter = obj.frame_counter + 1;
%             end % create movie
       

%         if obj.create_movie
%             writeVideo(video,obj.struct_movie);
%             struct_movie_save = obj.struct_movie;
%             save(['save_files/lineout/struct',obj.datadir,...
%                 '_',obj.wakefields_direction,'.mat'],'struct_movie_save')
%             close(video);
%         end % if create movie



    end % plot lineout2

end % ordinary methods


methods(Static)


end % static methods

end % classdef