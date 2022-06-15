%__________________________________________________________________________
% Superclass which loads the data from the Osiris HDF5 or MAT files and
% creates the correspoding variables for the other classes.
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/07/2020
%
% To-do: decide which properties should be hidden
%__________________________________________________________________________

% Input: options to get the specific data (specified below)
%
% Output: data matrix, ntime, sim. window length, and axes
%
% Methods
% - getdata: gets the data
%
%

% pointer to add more species to lcode2mat: 2x3g56jh

classdef OsirisLoader < handle

    properties(Constant, Hidden)

        % Physical constants
        c_m                 = physconst('LightSpeed'); % lightspeed, m/s
        c_cm                = physconst('LightSpeed')*100; % lightspeed in cm , cm/s
        e_charge_C          = 1.60217662e-19; % electron charge, C
        e_mass_kg           = 9.1093837015e-31; % electron mass, kg
        p_mass_kg           = 1.67262192369e-27; % proton mass, kg
        e_mass_eV           = 0.51099895; % electron mass in MeV/c^2
        p_mass_eV           = 938.2720813; % proton mass in MeV/c^2

        % note: dimensional analysis says F = s^2C^2 / m^2 kg
        permittivity        = 8.8541878128e-12; % permittivity of space, F/m
        permittivity_cm     = 8.8541878128e-14; % permittivity of space, F/cm
        permittivity_cmx    = 8.8541878128e-18; % permittivity of space, s^2C^2 / cm^2 kg / cm

        % values from: https://physics.nist.gov/cuu/Constants/index.html

    end % constant properties

    properties(Hidden)
        % hidden input, description in parser
        openrange; justPath; partialpath; fullpath;

        % h52mat
        downsample_z;

        % step size (normalized)
        ndr; ndz; % 2*beampop*e*c/IA / (dz/k) / N *np.ones(N)

        % track data
        tracks_time;
        tracks_q;
        tracks_ene;
        tracks_z;
        tracks_r;
        tracks_pz;
        tracks_pr;
        tracks_pth;
        trackfile_suffix;

        % lcode2mat
        lcode_timestep;
        lcode_dump_factor;
        lcode_output_time_period;
        lcode_time;
        ans_path = 1;
        load_lineout = 0;

        % momentum special temporal fix
        p2lineout_p;
        p2lineout_n;



    end

    properties

        % input, description in parser
        datadir; dataformat; useAvg;
        property; species; field; raw_dataset; track_dataset; direction;
        dump; wakefields_direction; lineout_point; lineout_direction;
        plasmaden; phasespace_direction;

        % output (all in norm. units (n prefix))
        ndataOut; % matrix with the desired data
        ntime; % ntime of the dump in norm. units
        nz; nr; % long. and trans. axes
        n_simulation_window; % simulation window
        n_propagation_distance; % propagation distance
        n_first_microbunch_position; % first microbunch position in normalized units
        first_microbunch_position; % first microbunch position in cm
        first_microbunch_position_px; % first microbunch position in pixels (for experimental data)

        % fields (normalized)
        nlongfield; % longitudinal wakefields (Ez)
        ntransfield; % transverse wakefields (Er-Bth)

        % species density (normalized)
        nproton_beam; % proton bunch density
        nplasma_electrons; % plasma electrons density
        nelectron_bunch; % accelerated electron bunch density
        nelectron_seed; % accelerated electron bunch density
        nplasma_positrons; % plasma positron density
        nelectrons;
        ndensity_feature;
        nantiproton_beam;

        %phasespace
        npr;

        % more species to come ...

        % raw_data
        q_raw; % weights of each particle in the raw data (osiris), normalized charge of particles (lcode)
        tag; % tag of particle
        w_raw; % weight of each particle (lcode)

        % raw_data (normalized)
        nz_raw; % long. pos. of particle
        nr_raw; % trans. pos. of particle
        npz_raw; % momentum in z of particle
        npr_raw; % momentum in r of particle
        npth_raw; % momentum in theta of particle
        nE_raw; % energy of particle
        nlineout; % lineout with desired specifications

    end % properties

    methods

        % --CONSTRUCTOR--

        function obj = OsirisLoader(varargin)

            % obj.CreatePath()

            % parse input to load files
            p = inputParser;

            % Directory where the MAT or MS folder is located
            p.addParameter('datadir', 'baseline', @(x) ischar(x));
            % Data format: h5 or mat
            p.addParameter('dataformat', 'mat', @(x) ismember(x,{'h5','mat'}));
            % Downsamples the data when transforming to  mat
            p.addParameter('downsample_z', 1, @(x) isfloat(x))
            % Flag to return only the full path where the file is located and not open
            % any of its contents
            p.addParameter('justPath', false, @(x) islogical(x) || x == 0 || x == 1);
            % Flag to use the avg of the files given by Osiris (if it is available)
            p.addParameter('useAvg', false, @(x) islogical(x) || x == 0 || x == 1);
            % Work in progress, ignore this one
            p.addParameter('openrange', [1 1 inf inf], @(x) isfloat(x));

            % What folder from MS or MAT to open: fields, density or raw
            p.addParameter('property', 'fields', @(x) ismember(x,{'fields',...
                'density','raw','phasespace','tracks'}));
            % For density or raw, specify the species name
            p.addParameter('species', 'proton_beam', @(x) ismember(x,{'proton_beam',...
                'electrons','electron_bunch','electron_seed','electron_beam',...
                'antiproton_seed','proton_beamfront','plasma_electrons','plasma_positrons',...
                'density_feature','antiproton_beam'}));

            % For the fields, specify magnetic (b) or electrin (e)
            p.addParameter('field', 'e', @(x) ismember(x,{'e','b'}));
            % For the raw dataset, specify which data: energy - ene, momentum - p,
            % position - x, charge - , tag - tag
            p.addParameter('raw_dataset', '', @(x) ismember(x,{'ene','p','x','q','tag'}));
            % For the track data, specify which data: enersgy - ene, momentum - p,
            % position - x, charge - , tag - tag
            p.addParameter('track_dataset', '', @(x) ismember(x,{'ene','p','x','q','tag'}));
            % For the fields or the positions and momenta in the raw data, specify
            % direction
            p.addParameter('direction', 'z', @(x) ismember(x,{'r','z','\theta','f'}));
            % Specify dump number
            p.addParameter('dump', 0, @(x) isfloat(x) && x >=0);
            % load wakefields going in the long. or trans. direction
            % (mainly do automatically the ndataOut = Er - Bth for trans fields)
            % empty means don't do that
            p.addParameter('wakefields_direction', '', @(x) ismember(x,{'long','trans',''}));

            % for lineouts. direction
            p.addParameter('lineout_direction', 'long', @(x) ismember(x,{'long','trans'}));
            % depending on direction, point away from axis (long), or away from front of simulation window (trans)
            p.addParameter('lineout_point', 2, @(x) (isfloat(x) && x > 0) || str2double(x) > 0);

            % Plasma density
            p.addParameter('plasmaden', 1, @(x) isfloat(x) && x >=0);

            % Tracks file
            p.addParameter('trackfile_suffix', '', @(x) ischar(x));

            % Plasma density
            p.addParameter('lcode_dump_factor', 1, @(x) isfloat(x) && x >=0);


            p.parse(varargin{:});

            obj.datadir        = p.Results.datadir;
            obj.dataformat     = p.Results.dataformat;
            obj.downsample_z   = p.Results.downsample_z;
            obj.justPath       = p.Results.justPath;
            obj.useAvg         = p.Results.useAvg;
            obj.openrange      = p.Results.openrange;
            obj.property       = p.Results.property;
            obj.species        = p.Results.species;
            obj.field          = p.Results.field;
            obj.raw_dataset    = p.Results.raw_dataset;
            obj.track_dataset  = p.Results.track_dataset;
            obj.direction      = p.Results.direction;
            obj.dump           = p.Results.dump;
            obj.plasmaden      = p.Results.plasmaden;
            obj.wakefields_direction...
                = p.Results.wakefields_direction;
            obj.lineout_direction...
                = p.Results.lineout_direction;
            obj.lineout_point  = p.Results.lineout_point;
            obj.trackfile_suffix ...
                = p.Results.trackfile_suffix;
            obj.lcode_dump_factor ...
                = p.Results.lcode_dump_factor;

        end % constructor

        function obj = getdata(obj)

            % Names of files (with correct numbering)
            dump_char = sprintf('%06.6d',obj.dump);

            % override parameters if property is fields and wakefields
            % direction is set
            if strcmp(obj.property,'fields')

                switch obj.wakefields_direction
                    case 'trans'
                        if ~strcmp(obj.field,'e') || ~strcmp(obj.direction,'r')
                            warning('wakefields_direction is set. Overriding field and direction.')
                            obj.field = 'e';
                            obj.direction = 'r';
                        end
                    case 'long'
                        if ~strcmp(obj.field,'e') || ~strcmp(obj.direction,'z')
                            warning('wakefields_direction is set. Overriding field and direction.')
                            obj.field = 'e';
                            obj.direction = 'z';
                        end
                end % switch wakefields dir
            end % property == fields

            switch obj.direction
                case 'z'
                    temp_direction = '1';
                case 'r'
                    temp_direction = '2';
                case '\theta'
                    temp_direction = '3';
            end

            avg = ''; if obj.useAvg; avg = '-savg'; end

            switch obj.dataformat
                case 'h5'
                    format_dir = 'MS';
                case 'mat'
                    format_dir = 'MAT';
            end

            switch obj.species
                case 'plasma_electrons'
                    species_name = 'electrons';
                case 'plasma_positrons'
                    species_name = 'positrons';
                otherwise
                    species_name = obj.species;

            end

            switch obj.property
                case 'fields'
                    property_name_file = [obj.field,temp_direction];
                    obj.partialpath = [obj.datadir,'/',format_dir,'/FLD/',property_name_file,avg,'/'];
                    obj.fullpath = [obj.partialpath,property_name_file,avg,'-',dump_char,'.',obj.dataformat];

                case {'density','raw','phasespace'}

                    switch obj.property
                        case 'density'
                            property_name_file = 'charge';
                            filename = ['charge',avg,'-',species_name];
                            obj.partialpath = [obj.datadir,'/',format_dir,'/DENSITY/',species_name,'/charge',avg,'/'];
                        case 'raw'
                            if ismember(obj.raw_dataset,{'ene','q','tag','w'})
                                temp_direction = '';
                            end
                            property_name_file = [obj.raw_dataset,temp_direction];
                            filename = ['RAW-',obj.species];
                            obj.partialpath = [obj.datadir,'/',format_dir,'/RAW/',obj.species,'/'];
                            if strcmp(obj.dataformat,'mat')
                                obj.partialpath = [obj.partialpath,property_name_file,'/'];
                                filename = [filename,'-',property_name_file];
                            end
                        case 'phasespace'
                            property_name_file = 'phasespace';
                            l = '';
                            if obj.load_lineout; l = 'l'; end
                            filename = ['PHASESPACE',avg,'-',species_name,'-','p2',l];

                            obj.partialpath = [obj.datadir,'/',format_dir,'/PHASESPACE/',species_name,'/p2',avg,'/'];

                    end % switch property

                    obj.fullpath = [obj.partialpath,filename,'-',dump_char,'.',obj.dataformat];
                case 'tracks'
                    obj.partialpath = [obj.datadir,'obj.include_phasespace_p/',format_dir,'/TRACKS/'];
                    obj.fullpath = [obj.partialpath,obj.species,'-tracks-repacked',obj.trackfile_suffix,'.',obj.dataformat];
            end % switch property


            % Return just the path without reading or loading the mat or hdf5 files
            if obj.justPath

                return;
            end


            switch format_dir
                case 'MS'
                    MS_parameters = h5info(obj.fullpath);

                    % Apply 90 rotation counterclockwise and flip
                    switch obj.property
                        case {'fields','density','phasespace'}
                            obj.ntime = h5readatt(obj.fullpath,'/','TIME');
                            MS_size(1) = MS_parameters.Datasets.Dataspace.Size(1);
                            MS_size(2) = MS_parameters.Datasets.Dataspace.Size(2);
                            obj.ndataOut = h5read(obj.fullpath,['/',property_name_file],obj.openrange(1:2),obj.openrange(3:4));
                            if strcmp(obj.wakefields_direction,'trans') && strcmp(obj.property,'fields')
                                obj.fullpath = [obj.datadir,'/',format_dir,'/FLD/','b3',avg,'/','b3',avg,'-',dump_char,'.',obj.dataformat];
                                Bth = h5read(obj.fullpath,['/','b3'],obj.openrange(1:2),obj.openrange(3:4));
                                obj.ndataOut = obj.ndataOut - Bth;
                            end
                            z_startend = double(h5read(obj.fullpath,'/AXIS/AXIS1'))';
                            r_startend = double(h5read(obj.fullpath,'/AXIS/AXIS2'))';
                            r_startend = r_startend + abs(r_startend(1));
                            % for some reason, the first element in the r axis is not 0, but a negative value of 1/2 cell or so, so this sets it to 0.
                            % https://osirisdoc.wimpzilla.ist.utl.pt/dev/index.php/Reference_Guide:_Space
                            z_temp = linspace(z_startend(1),z_startend(2),MS_size(1));
                            obj.ndz = z_temp(2) - z_temp(1);
                            obj.nz = z_temp(obj.openrange(1):min(size(obj.ndataOut,1),obj.openrange(3)));
                            r_temp = linspace(r_startend(1),r_startend(2),MS_size(2));
                            obj.ndr = r_temp(2) - r_temp(1); % Actually r_temp(2) should be correct, but ok, just to be sure.
                            obj.nr = r_temp(obj.openrange(2):min(size(obj.ndataOut,2),obj.openrange(4)));
                            obj.ndataOut = obj.ndataOut.';
                            obj.n_simulation_window = diff(z_startend);
                            obj.n_propagation_distance = mean(z_startend) - obj.n_simulation_window;
                        case 'raw'
                            obj.ntime = h5readatt(obj.fullpath,'/','TIME');
                            obj.ndataOut = h5read(obj.fullpath,['/',property_name_file]);
                        case 'tracks'
                            tracks_data_temp = h5read(obj.fullpath,'/data');
                            tracks_iter_temp = h5read(obj.fullpath,'/itermap');

                            nsave = max(tracks_iter_temp(2,:)); % data was dumped each niter*nsave iterations (50)
                            itersave = unique(tracks_iter_temp(3,:)); % first iteration at which the dump was saved (1x20)
                            ntracks_per_particle = length(itersave)*nsave; % track points per particle (1000)
                            npar = max(tracks_iter_temp(1,:));
                            iter_length = length(tracks_iter_temp);

                            % initialize tables to memory
                            tracks_in_order = zeros(npar,8,ntracks_per_particle);
                            ind = ones(npar,2); ind(:,2) = zeros(npar,1);
                            begin = 1; final = 0;

                            %                             indtable = ones(npar,2*ntracks_per_particle/nsave);
                            %                             indtable(:,2) = zeros(npar,1);
                            %                             par_occurance_table = zeros(npar,1);

                            for ii = 1:iter_length

                                par = tracks_iter_temp(1,ii);

                                final = final + tracks_iter_temp(2,ii);
                                ind(par,2) = ind(par,2) + tracks_iter_temp(2,ii);

                                %                                 indtable(par,par_occurance_table(par)*2+2) = ...
                                %                                     indtable(par,(par_occurance_table(par)-1)*2+2) + tracks_iter_temp(2,ii);


                                tracks_in_order(par,:,ind(par,1):ind(par,2))...
                                    = tracks_data_temp(:,begin:final);

                                begin = final + 1;
                                ind(par,1) = ind(par,2) + 1;

                                %                                 par_occurance_table(par) = par_occurance_table(par) + 1;
                                %                                 indtable(par,par_occurance_table(par)*2+1) = ...
                                %                                     indtable(par,(par_occurance_table(par)-1)*2+2) + 1;

                                if mod(ii,round(0.01*iter_length)) == 0
                                    fprintf('PROGRESS REARRANGING TRACKS: %.d %% \n',round(ii/iter_length*100))
                                end % if fprintf
                            end % for iter length

                            if strcmp(obj.fullpath,'gm20/MS/TRACKS/proton_beam-tracks-repacked_n130repeated.h5')
                                tracks_in_order(:,:,251:300) = [];
                            end

                            clear tracks_data_temp tracks_iter_temp

                            tracks_obj.partialpath_handle = what(obj.partialpath);
                            tracks_partial_save_filename = [tracks_obj.partialpath_handle.path(1:end-9),'MAT\TRACKS'];
                            tracks_full_save_filename = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked',obj.trackfile_suffix,'.mat'];
                            tracks_full_save_filename_time = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked-time',obj.trackfile_suffix,'.mat'];
                            tracks_full_save_filename_q = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked-q',obj.trackfile_suffix,'.mat'];
                            tracks_full_save_filename_ene = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked-ene',obj.trackfile_suffix,'.mat'];
                            tracks_full_save_filename_z = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked-z',obj.trackfile_suffix,'.mat'];
                            tracks_full_save_filename_r = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked-r',obj.trackfile_suffix,'.mat'];
                            tracks_full_save_filename_pz = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked-pz',obj.trackfile_suffix,'.mat'];
                            tracks_full_save_filename_pr = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked-pr',obj.trackfile_suffix,'.mat'];
                            tracks_full_save_filename_pth = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked-pth',obj.trackfile_suffix,'.mat'];

                            if isfile(tracks_full_save_filename)
                                obj.tracks_time  = squeeze(tracks_in_order(:,1,:));
                                obj.tracks_q = squeeze(tracks_in_order(:,2,:));
                                obj.tracks_ene = squeeze(tracks_in_order(:,3,:));
                                obj.tracks_z = squeeze(tracks_in_order(:,4,:));
                                obj.tracks_r = squeeze(tracks_in_order(:,5,:));
                                obj.tracks_pz = squeeze(tracks_in_order(:,6,:));
                                obj.tracks_pr = squeeze(tracks_in_order(:,7,:));
                                obj.tracks_pth = squeeze(tracks_in_order(:,8,:));
                            else

                                % distribute data intro matrices
                                tracks_time  = squeeze(tracks_in_order(:,1,:)); %#ok<PROP>
                                tracks_q = squeeze(tracks_in_order(:,2,:)); %#ok<PROP>
                                tracks_ene = squeeze(tracks_in_order(:,3,:)); %#ok<PROP>
                                tracks_z = squeeze(tracks_in_order(:,4,:)); %#ok<PROP>
                                tracks_r = squeeze(tracks_in_order(:,5,:)); %#ok<PROP>
                                tracks_pz = squeeze(tracks_in_order(:,6,:)); %#ok<PROP>
                                tracks_pr = squeeze(tracks_in_order(:,7,:)); %#ok<PROP>
                                tracks_pth = squeeze(tracks_in_order(:,8,:)); %#ok<PROP>

                                if ~isfolder(tracks_partial_save_filename); mkdir(tracks_partial_save_filename); end
                                save(tracks_full_save_filename,'tracks_time','tracks_q',...
                                    'tracks_ene','tracks_z','tracks_r',...
                                    'tracks_pz','tracks_pr','tracks_pth');
                                %                                 save(tracks_full_save_filename_time,'tracks_time');
                                %                                 save(tracks_full_save_filename_q,'tracks_q');
                                %                                 save(tracks_full_save_filename_ene,'tracks_ene');
                                %                                 save(tracks_full_save_filename_z,'tracks_z');
                                %                                 save(tracks_full_save_filename_r,'tracks_r');
                                %                                 save(tracks_full_save_filename_pz,'tracks_pz');
                                %                                 save(tracks_full_save_filename_pr,'tracks_pr');
                                %                                 save(tracks_full_save_filename_pth,'tracks_pth');

                                obj.tracks_time  = squeeze(tracks_in_order(:,1,:));
                                obj.tracks_q = squeeze(tracks_in_order(:,2,:));
                                obj.tracks_ene = squeeze(tracks_in_order(:,3,:));
                                obj.tracks_z = squeeze(tracks_in_order(:,4,:));
                                obj.tracks_r = squeeze(tracks_in_order(:,5,:));
                                obj.tracks_pz = squeeze(tracks_in_order(:,6,:));
                                obj.tracks_pr = squeeze(tracks_in_order(:,7,:));
                                obj.tracks_pth = squeeze(tracks_in_order(:,8,:));

                            end


                    end % if switch property

                case 'MAT'
                    data = matfile(obj.fullpath);

                    if obj.load_lineout
                        data.Properties.Writable = true;
                        data.phasespace = [];
                        obj.p2lineout_p = data.p2lineout_p;
                        obj.p2lineout_n = data.p2lineout_n;
                    end % load lineout

                    switch obj.property
                        case {'fields','density','phasespace'}
                            obj.ndataOut = data.(obj.property);
                            if strcmp(obj.wakefields_direction,'trans') && strcmp(obj.property,'fields')
                                obj.fullpath = [obj.datadir,'/',format_dir,'/FLD/','b3',avg,'/','b3',avg,'-',dump_char,'.',obj.dataformat];
                                Bth = matfile(obj.fullpath);
                                obj.ndataOut = obj.ndataOut - Bth.fields;
                            end
                            z_startend = data.axis1;
                            obj.nz = data.x1_axis;
                            obj.ndz = obj.nz(2) - obj.nz(1);
                            nr_temp = data.x2_axis;
                            obj.nr = nr_temp + abs(nr_temp(1));
                            obj.ndr = obj.nr(2) - obj.nr(1); % Actually obj.nr(2) should be correct, but ok, just to be sure.
                            % for some reason, the first element in the r axis is not 0, but a negative value of 1/2 cell, so this sets it to 0.
                            % https://osirisdoc.wimpzilla.ist.utl.pt/dev/index.php/Reference_Guide:_Space
                            obj.n_simulation_window = diff(z_startend);
                            obj.n_propagation_distance = mean(z_startend) - obj.n_simulation_window;
                            obj.ntime = data.time;
                        case 'raw'
                            obj.ndataOut = data.rawdata;
                            obj.ntime = data.time;
                        case 'tracks'
                            tracks_partialpath_handle = what(obj.partialpath);
                            tracks_partial_save_filename = tracks_partialpath_handle.path;
                            tracks_full_save_filename = [tracks_partial_save_filename,'\',obj.species,'-tracks-repacked.mat'];
                            if ~isfile(tracks_full_save_filename); error('Please try running the code once with hdf5'); end

                            switch obj.track_dataset
                                case 'q'
                                    obj.tracks_q = data.tracks_q;
                                case 'ene'
                                    obj.ene = data.ene;
                                case 'tag'
                                    obj.tag = data.tag;
                                case 'x'
                                    switch obj.direction
                                        case 'z'
                                            obj.tracks_z = data.tracks_z;
                                        case 'r'
                                            obj.tracks_r= data.tracks_r;
                                    end % switch direction
                                case 'p'
                                    switch obj.direction
                                        case 'z'
                                            obj.tracks_pz = data.tracks_pz;
                                        case 'r'
                                            obj.tracks_pr = data.tracks_pr;
                                        case '\theta'
                                            obj.tracks_pth = data.tracks_pth;
                                    end % switch direction
                            end % switch raw dataset

                            %                             obj.tracks_time  = data.tracks_time;
                            %                             obj.tracks_q = data.tracks_q;
                            %                             obj.tracks_ene = data.tracks_ene;
                            %                             obj.tracks_z = data.tracks_z;
                            %                             obj.tracks_r = data.tracks_r;
                            %                             obj.tracks_pz = data.tracks_pz;
                            %                             obj.tracks_pr = data.tracks_pr;
                            %                             obj.tracks_pth = data.tracks_pth;

                    end % switch property


            end % switch formatdir
        end % getdata

        function obj = getnlineogetdatout(obj)
            % get lineout of the data trimming is not needed
            obj.getdata();

            % AXIS
            switch obj.lineout_direction
                case 'long'
                    obj.nlineout = obj.ndataOut(obj.lineout_point,:);
                case 'trans'
                    obj.nlineout = obj.ndataOut(:,obj.lineout_point);
            end

        end % getnlineout

        function obj = assign_longfields(obj)
            obj.nlongfield = obj.ndataOut;
            obj.ndataOut = [];
        end % assign_longfields

        function obj = assign_transfields(obj)
            obj.ntransfield = obj.ndataOut;
            obj.ndataOut = [];
        end % assign_transfields

        function obj = assign_fields(obj)
            obj.(['n',obj.wakefields_direction,'field']) = obj.ndataOut;
            obj.ndataOut = [];
        end

        function obj = assign_density(obj)
            switch obj.species
                case {'proton_beam','proton_bunch','density_feature'}
                    obj.nproton_beam = abs(obj.ndataOut);
                    %                     obj.(['n',obj.species]) = obj.ndataOut;
                case 'plasma_electrons'
                    obj.nplasma_electrons = abs(obj.ndataOut);
                case 'plasma_positrons'
                    obj.nplasma_positrons = obj.ndataOut;
                otherwise
                    obj.(['n',obj.species]) = abs(obj.ndataOut);
            end %switch
            obj.ndataOut = [];
        end % assign_density

        function obj = assign_raw(obj)
            switch obj.raw_dataset
                case 'q'
                    obj.q_raw = obj.ndataOut;
                case 'ene'
                    obj.nE_raw = obj.ndataOut;
                case 'tag'
                    obj.tag = obj.ndataOut;
                case 'w'
                    obj.w_raw = obj.ndataOut;
                case 'x'
                    switch obj.direction
                        case 'z'
                            obj.nz_raw = obj.ndataOut;
                        case 'r'
                            obj.nr_raw = obj.ndataOut;
                    end % switch direction
                case 'p'
                    switch obj.direction
                        case 'z'
                            obj.npz_raw = obj.ndataOut;
                        case 'r'
                            obj.npr_raw = obj.ndataOut;
                        case '\theta'
                            obj.npth_raw = obj.ndataOut;
                    end % switch direction
            end % switch raw dataset
            obj.ndataOut = [];
        end %assign_raw

        function h52mat(obj)
            obj.getdata();

            avg = '';
            if obj.useAvg; avg = '-savg'; end

            data_save = obj.ndataOut(:,1:obj.downsample_z:end);
            nx1 = size(data_save,2);
            nx2 = size(data_save,1);

            partialpath_handle = what(obj.partialpath);
            save_partialpath = strrep(partialpath_handle.path,'MS','MAT');

            switch obj.direction
                case 'z'; temp_direction = '1';
                case 'r'; temp_direction = '2';
                case '\theta'; temp_direction = '3';
            end



            switch obj.property
                case 'raw'
                    if strcmp(obj.raw_dataset,'ene') || strcmp(obj.raw_dataset,'q') || strcmp(obj.raw_dataset,'tag')
                        temp_direction = '';
                    end
                    attfilename = ['RAW','-',obj.species,'-','attribute','.mat'];
                    if ~isfile([save_partialpath,attfilename])
                        attval = h5readatt(obj.fullpath,'/','NAME');
                        NAME = attval{1};
                        attval = h5readatt(obj.fullpath,'/','TYPE');
                        TYPE = attval{1};
                        SELECT_GAMMA_LIMIT = h5readatt(obj.fullpath,'/','SELECT GAMMA LIMIT');
                        SELECT_FRACTION = h5readatt(obj.fullpath,'/','SELECT FRACTION');
                        attval = h5readatt(obj.fullpath,'/','SELECT MATH EXPR');
                        SELECT_MATH_EXPR = attval{1};
                        DT = h5readatt(obj.fullpath,'/','DT');
                        attval = h5readatt(obj.fullpath,'/','TIME UNITS');
                        TIME_UNITS = attval{1};
                        PERIODIC = h5readatt(obj.fullpath,'/','PERIODIC');
                        MOVE_C = h5readatt(obj.fullpath,'/','MOVE C');
                        PAR_NODE_CONF = h5readatt(obj.fullpath,'/','PAR_NODE_CONF');
                        PAR_NX_X1 = h5readatt(obj.fullpath,'/','PAR_NX_X1');
                        PAR_NX_X2 = h5readatt(obj.fullpath,'/','PAR_NX_X2');
                    end
                    attpropertyfilename = ['RAW','-',obj.species,'-',obj.raw_dataset,'-attribute','.mat'];
                    if (~strcmp(obj.raw_dataset,'tag')) && (~isfile([save_partialpath,obj.raw_dataset,'/',attpropertyfilename]))
                        attval = h5readatt(obj.fullpath,['/',obj.raw_dataset,temp_direction],'UNITS');
                        UNITS = attval{1};
                    end

                otherwise

                    ITER = h5readatt(obj.fullpath,'/','ITER');
                    axis1 = h5read(obj.fullpath,'/AXIS/AXIS1')';
                    axis2 = h5read(obj.fullpath,'/AXIS/AXIS2')';
                    x1_axis = linspace(axis1(1),axis1(2),nx1);
                    x2_axis = linspace(axis2(1),axis2(2),nx2);
                    switch obj.property
                        case 'fields'
                            attfilename = [obj.property,avg,'-','attribute','.mat'];
                        case 'density'
                            attfilename = [obj.property,avg,'-',obj.species,'-','attribute','.mat'];
                            temp_direction = ''; obj.field = 'charge'; % a bit confusing, but less code, it is just to find the file.
                    end

                    if ~isfile([save_partialpath,attfilename])
                        attval = h5readatt(obj.fullpath,'/','NAME');
                        NAME = attval{1};
                        attval = h5readatt(obj.fullpath,'/','TYPE');
                        TYPE = attval{1};
                        attval = h5readatt(obj.fullpath,'/','TIME UNITS');
                        TIME_UNITS = attval{1};
                        PERIODIC = h5readatt(obj.fullpath,'/','PERIODIC');
                        MOVE_C = h5readatt(obj.fullpath,'/','MOVE C');
                        PAR_NODE_CONF = h5readatt(obj.fullpath,'/','PAR_NODE_CONF');
                        PAR_NX_X1 = h5readatt(obj.fullpath,'/','PAR_NX_X1');
                        PAR_NX_X2 = h5readatt(obj.fullpath,'/','PAR_NX_X2');
                        attval = h5readatt(obj.fullpath,['/',obj.field,temp_direction],'UNITS');
                        CHARGE_UNITS = attval{1};
                        attval = h5readatt(obj.fullpath,'/AXIS/AXIS1','UNITS');
                        AXIS1_UNITS = attval{1};
                        attval = h5readatt(obj.fullpath,'/AXIS/AXIS1','TYPE');
                        AXIS1_TYPE = attval{1};
                        attval = h5readatt(obj.fullpath,'/AXIS/AXIS2','UNITS');
                        AXIS2_UNITS = attval{1};
                        attval = h5readatt(obj.fullpath,'/AXIS/AXIS2','TYPE');
                        AXIS2_TYPE = attval{1};
                    end
            end
            time = obj.ntime;
            save_fullpath = strrep(which(obj.fullpath),'MS','MAT');
            save_fullpath = strrep(save_fullpath,'h5','mat');

            switch obj.property
                case {'fields', 'density', 'phasespace'}
                    if ~isfolder(save_partialpath)
                        mkdir(save_partialpath)
                    end
                    switch obj.property % don't erase this one, should be like this
                        case 'fields'; fields = data_save;
                        case 'density'; density = data_save;
                        case 'phasespace'; phasespace = data_save;
                    end
                    save(save_fullpath,...
                        obj.property,'ITER','time','axis1','axis2','x1_axis','x2_axis','-v6');
                    if ~isfile([save_partialpath,attfilename])
                        save([save_partialpath,attfilename],'NAME','TYPE',...
                            'TIME_UNITS','PERIODIC','MOVE_C','PAR_NODE_CONF','PAR_NX_X1','PAR_NX_X2',...
                            'CHARGE_UNITS','AXIS1_UNITS','AXIS1_TYPE','AXIS2_UNITS','AXIS2_TYPE');
                    end
                case 'raw'
                    savepath = [save_partialpath,'/',obj.raw_dataset,temp_direction,'/'];
                    if ~isfolder(savepath)
                        mkdir(savepath)
                    end
                    dump_char = sprintf('%06.6d',obj.dump);
                    filename = ['RAW','-',obj.species,'-',obj.raw_dataset,temp_direction,'-',dump_char,'.mat'];
                    rawdata = data_save;
                    save([savepath,filename],'rawdata','time','-v6');
                    if ~isfile([save_partialpath,attfilename])
                        save([save_partialpath,attfilename],'NAME','TYPE','SELECT_GAMMA_LIMIT',...
                            'DT','TIME_UNITS','PERIODIC','MOVE_C','PAR_NODE_CONF','PAR_NX_X1','PAR_NX_X2');
                    end
                    if (~strcmp(obj.raw_dataset,'tag')) && (~isfile([save_partialpath,obj.raw_dataset,'/',attpropertyfilename]))
                        save([save_partialpath,obj.raw_dataset,'/',attpropertyfilename],'UNITS');
                    end
            end


        end % h52mat

        function lcode2mat(obj)

            % input
            n_ave = 1;

            % xi_lcode = 600;
            % r_lcode = 3;

            % read file
            lcode_input = fileread([obj.datadir,'.cfg']);
            lcode_timestep_temp = regexp(lcode_input,'time-step = (\d+\.?\d*)','tokens');
            obj.lcode_timestep = str2double(lcode_timestep_temp{1});
            lcode_output_time_period_temp = regexp(lcode_input,'output-time-period = (\d+\.?\d*)','tokens');
            if ~isempty(lcode_output_time_period_temp) && (str2double(lcode_output_time_period_temp{1}) > obj.lcode_timestep)
                obj.lcode_output_time_period = str2double(lcode_output_time_period_temp{1});
            else
                obj.lcode_output_time_period = obj.lcode_timestep;
            end
            lcode_xi_temp = regexp(lcode_input,'window-length = (\d+\.?\d*)','tokens');
            lcode_xi = str2double(lcode_xi_temp{1});
            lcode_dxi_temp = regexp(lcode_input,'xi-step = (\d+\.?\d*)','tokens');
            lcode_dxi = str2double(lcode_dxi_temp{1});
            nx1 = lcode_xi/lcode_dxi;
            lcode_r_temp = regexp(lcode_input,'window-width = (\d+\.?\d*)','tokens');
            lcode_r = str2double(lcode_r_temp{1});
            lcode_dr_temp = regexp(lcode_input,'r-step = (\d+\.?\d*)','tokens');
            lcode_dr = str2double(lcode_dr_temp{1});
            nx2 = lcode_r/lcode_dr+1;

            obj.lcode_time = round(round(obj.lcode_output_time_period*obj.dump/obj.lcode_dump_factor/obj.lcode_timestep)*obj.lcode_timestep - eps*round(obj.lcode_output_time_period*obj.dump/obj.lcode_dump_factor/obj.lcode_timestep)*obj.lcode_timestep);

            openID = obj.get_openID();

            if openID == -1
                obj.lcode_time = obj.lcode_time + 1;
                openID = obj.get_openID();
            end

            if openID == -1
                obj.lcode_time = obj.lcode_time - 2;
                openID = obj.get_openID();
            end

            if strcmp(obj.property,'raw')
                obj.ndataOut = fread(openID,[8,inf],'double')';
                % z,r,pz,pr,pa,q,w,id
            else
                data_temp = sscanf(fread(openID, [1, inf], '*char'), '%f');
                obj.ndataOut  = single(fliplr(reshape(data_temp,[nx2,nx1])));
                data_save = obj.ndataOut(1:end-1,1:obj.downsample_z:end);

                if n_ave > 1
                    data_save4D = reshape(data_save,2,size(data_save,1)/2,...
                        2,size(data_save,2)/2);
                    data_save_temp = sum(sum(data_save4D, 1), 3)*0.25;
                    data_save = reshape(data_save_temp, size(data_save,1)/2, size(data_save,2)/2);
                end
            end % si raw

            fclose(openID);

            partialpath_handle = what(obj.datadir);

            %             if (length(partialpath_handle) > 1) && (obj.ans_path == 1)
            %                 obj.ans_path = inputdlg('Two locations of the directory, which one is correct? (input number)');
            %             else
            %                 obj.ans_path = 1;
            %             end
            %             save_partialpath = [partialpath_handle(str2double(obj.ans_path{1})).path,'/MAT/'];

            save_partialpath = [partialpath_handle.path,'/MAT/'];

            switch obj.direction
                case 'z'; temp_direction = '1';
                case 'r'; temp_direction = '2';
                case 'f'; temp_direction = '3';
            end

            if strcmp(obj.raw_dataset,'ene') || strcmp(obj.raw_dataset,'q') || strcmp(obj.raw_dataset,'tag')
                temp_direction = '';
            end

            % Names of files (with correct numbering)
            dump_char = sprintf('%06.6d',obj.dump);
            avg = ''; if obj.useAvg; avg = '-savg'; end

            axis1 = obj.lcode_time + [0,lcode_xi];
            axis2 = [0,lcode_r];
            x1_axis = linspace(axis1(1),axis1(2),nx1/n_ave);
            x2_axis = linspace(axis2(1),axis2(2),(nx2 - 1)/n_ave);
            time = obj.lcode_time;

            switch obj.property
                case 'fields'
                    property_name_file = [obj.field,temp_direction];
                    obj.partialpath = [save_partialpath,'FLD/',property_name_file,avg,'/'];
                    obj.fullpath = [obj.partialpath,property_name_file,avg,'-',dump_char,'.mat'];
                    fields = data_save;
                    if ~isfolder(obj.partialpath); mkdir(obj.partialpath); end
                    save(obj.fullpath,...
                        'fields','time','axis1','axis2','x1_axis','x2_axis'); % '-v6'

                case 'density'
                    filename = ['charge',avg,'-',obj.species];
                    obj.partialpath = [save_partialpath,'DENSITY/',obj.species,'/charge',avg,'/'];
                    obj.fullpath = [obj.partialpath,filename,'-',dump_char,'.mat'];
                    density = data_save;
                    if ~isfolder(obj.partialpath); mkdir(obj.partialpath); end
                    save(obj.fullpath,...
                        'density','time','axis1','axis2','x1_axis','x2_axis'); %'-v6'

                case 'raw' %

                    %{'ene','p','x','q','w','tag'}
                    % z,r,pz,pr,pa,q,w,id
                    % 1,2, 3, 4, 5,6,7, 8

                    species_s = {'proton_beam','electron_seed','antiproton_beam',''};

                    if ~ismember(obj.species,species_s)
                        species_s = obj.species;
                    end % if selecting another species loke electron_bunch, to have a corret name 
                    

                    raw_dataset_s = {'x','p','q','w','tag'};
                    % proton electron antiproton indices 
                    proton_index = (obj.ndataOut(:,7) > 0) & (obj.ndataOut(:,6) < 1);
                    electron_index = (obj.ndataOut(:,7) < 0) & (obj.ndataOut(:,6) > 0.99); % avoid rounding errors
                    antiproton_index = (obj.ndataOut(:,7) < 0) & (obj.ndataOut(:,6) < 1);

                    

                    for ss = 1:(length(species_s)-1)
                        obj.species = species_s{ss};
                        
                        % 2x3g56jh
                        switch obj.species
                            case 'proton_beam'
                                if ~any(proton_index); continue; end 
                                raw_data_temp = obj.ndataOut(proton_index,:);
                            case {'electron_seed','electron_bunch'}
                                if ~any(electron_index); continue; end 
                                raw_data_temp = obj.ndataOut(electron_index,:);
                            case 'antiproton_beam'
                                if ~any(antiproton_index); continue; end 
                                raw_data_temp = obj.ndataOut(antiproton_index,:);
                        end

                        ds = 1;

                        for raw = 1:4
                            obj.raw_dataset = raw_dataset_s{raw};


                            switch obj.raw_dataset
                                case {'ene','q'}
                                    temp_direction_s = {''};
                                case 'x'
                                    temp_direction_s = {'1','2'};
                                case 'p'
                                    temp_direction_s = {'1','2','3'};
                                case 'w'
                                    temp_direction_s = {''};
                                case 'tag'
                                    temp_direction_s = {''};
                            end % switch raw data set

                            for dd = 1:length(temp_direction_s)
                                temp_direction = temp_direction_s{dd};
                                data_save = raw_data_temp(:,ds);
                                ds = ds + 1;

                                obj.partialpath = [save_partialpath,'RAW/',obj.species,...
                                    '/',obj.raw_dataset,temp_direction,'/'];

                                if ~isfolder(obj.partialpath); mkdir(obj.partialpath); end

                                filename = ['RAW','-',obj.species,'-',obj.raw_dataset,temp_direction];
                                obj.fullpath = [obj.partialpath,filename,'-',dump_char,'.mat'];
                                rawdata = data_save;
                                save(obj.fullpath,'rawdata','time'); %'-v6'

                            end % for temp direction
                        end % for raW
                    end % for species

                case 'tracks' % solve later
                    obj.partialpath = [obj.datadir,'/',format_dir,'/TRACKS/'];
                    obj.fullpath = [obj.partialpath,obj.species,'-tracks-repacked',obj.trackfile_suffix,'.',obj.dataformat];
            end % switch property

        end %lcode2mat

        function read_lcode(obj)

            % input
            n_ave = 1;

            % xi_lcode = 600;
            % r_lcode = 3;

            % read file
            lcode_input = fileread([obj.datadir,'.cfg']);
            lcode_timestep_temp = regexp(lcode_input,'time-step = (\d+\.?\d*)','tokens');
            obj.lcode_timestep = str2double(lcode_timestep_temp{1});
            lcode_output_time_period_temp = regexp(lcode_input,'output-time-period = (\d+\.?\d*)','tokens');
            if ~isempty(lcode_output_time_period_temp) && (str2double(lcode_output_time_period_temp{1}) > obj.lcode_timestep)
                obj.lcode_output_time_period = str2double(lcode_output_time_period_temp{1});
            else
                obj.lcode_output_time_period = obj.lcode_timestep;
            end
            lcode_xi_temp = regexp(lcode_input,'window-length = (\d+\.?\d*)','tokens');
            lcode_xi = str2double(lcode_xi_temp{1});
            lcode_dxi_temp = regexp(lcode_input,'xi-step = (\d+\.?\d*)','tokens');
            lcode_dxi = str2double(lcode_dxi_temp{1});
            nx1 = lcode_xi/lcode_dxi;
            lcode_r_temp = regexp(lcode_input,'window-width = (\d+\.?\d*)','tokens');
            lcode_r = str2double(lcode_r_temp{1});
            lcode_dr_temp = regexp(lcode_input,'r-step = (\d+\.?\d*)','tokens');
            lcode_dr = str2double(lcode_dr_temp{1});
            nx2 = lcode_r/lcode_dr+1;

            obj.lcode_time = round(round(obj.lcode_output_time_period*obj.dump/obj.lcode_dump_factor/obj.lcode_timestep)*obj.lcode_timestep - eps*round(obj.lcode_output_time_period*obj.dump/obj.lcode_dump_factor/obj.lcode_timestep)*obj.lcode_timestep);

            openID = obj.get_openID();

            if openID == -1
                obj.lcode_time = obj.lcode_time + 1;
                openID = obj.get_openID();
            end

            if openID == -1
                obj.lcode_time = obj.lcode_time - 2;
                openID = obj.get_openID();
            end

            if strcmp(obj.property,'raw')
                obj.ndataOut = fread(openID,[8,inf],'double')';
                % z,r,pz,pr,pa,q,w,id
            else
                data_temp = sscanf(fread(openID, [1, inf], '*char'), '%f');
                obj.ndataOut  = single(fliplr(reshape(data_temp,[nx2,nx1])));
                data_save = obj.ndataOut(1:end-1,1:obj.downsample_z:end);

                if n_ave > 1
                    data_save4D = reshape(data_save,2,size(data_save,1)/2,...
                        2,size(data_save,2)/2);
                    data_save_temp = sum(sum(data_save4D, 1), 3)*0.25;
                    data_save = reshape(data_save_temp, size(data_save,1)/2, size(data_save,2)/2);
                end
            end

            fclose(openID);

             temp_path = what(obj.datadir);
             obj.partialpath = temp_path.path;


            % Names of files (with correct numbering)
            lcode_dump_char = sprintf('%05.5d',obj.lcode_time);


            switch obj.property
                case 'fields'
                    property_name_file = [obj.field,temp_direction];
                    obj.partialpath = [save_partialpath,'FLD/',property_name_file,avg,'/'];
                    obj.fullpath = [obj.partialpath,property_name_file,'-',dump_char,'.mat'];
                    fields = data_save;
                    if ~isfolder(obj.partialpath); mkdir(obj.partialpath); end
                    save(obj.fullpath,...
                        'fields','time','axis1','axis2','x1_axis','x2_axis'); % '-v6'

                case 'density'
                    filename = ['charge',avg,'-',obj.species];
                    obj.partialpath = [save_partialpath,'DENSITY/',obj.species,'/charge',avg,'/'];
                    obj.fullpath = [obj.partialpath,filename,'-',dump_char,'.mat'];
                    density = data_save;
                    if ~isfolder(obj.partialpath); mkdir(obj.partialpath); end
                    save(obj.fullpath,...
                        'density','time','axis1','axis2','x1_axis','x2_axis'); %'-v6'

                case 'raw' %

                    %{'ene','p','x','q','tag'}
                    % z,r,pz,pr,pa,q,w,id
                    % 1,2, 3, 4, 5,6,7, 8


                    if ~isfolder(obj.partialpath); mkdir(obj.partialpath); end

                    txt = 'tb';
                    obj.fullpath = [obj.partialpath,'/',txt,lcode_dump_char,'.swp'];
                    ndataOut_check = obj.ndataOut;
                    obj.ndataOut(:,4) = -1*obj.ndataOut(:,4);

                    ebeam = obj.ndataOut;
                    openID2 = fopen(['r2l_302_c_e550_dr_hd2','/','beamfile.bin'],'r');
                    pbeam = fread(openID2,[8,inf],'double')';
                    fclose(openID2);
                    ebeam(end,:) = [];
                    ind_ebeam = (length(pbeam)+1):1:(length(pbeam)+length(ebeam));
                    ebeam(:,8) = ind_ebeam';
                    two_beams = [ebeam;pbeam];
                    fileID = fopen(['beamfile.bin'],'w');
                    fwrite(fileID,two_beams,'single');
                    fclose(fileID);

            end % switch property

        end %lcode2mat

        function openID = get_openID(obj)

            lcode_dump_char = sprintf('%05.5d',obj.lcode_time);

            switch obj.property
                case 'fields'
                    openID = fopen([obj.datadir,'/',obj.field,obj.direction,lcode_dump_char,'.swp'],'r');
                case {'density'}
                    switch obj.species
                        case 'electrons'
                            species_letter = 'e';
                        case 'ions'
                            species_letter = 'i';
                        otherwise % proton_beam, electron_seed, electron_bunch, ...
                            species_letter = 'b';
                    end
                    openID = fopen([obj.datadir,'/n',species_letter,lcode_dump_char,'.swp'],'r');
                case 'raw'
                    switch obj.species
                        case 'electrons'
                            txt = 'pl';
                        case 'ions'
                            % nothing
                        otherwise % proton_beam, electron_seed, electron_bunch, ...
                            txt = 'tb';
                    end
                    openID = fopen([obj.datadir,'/',txt,lcode_dump_char,'.swp'],'r');

            end % switch property


        end % get_openID


    end % ordinary methods


    methods(Static)

        function progress_dump(ptext,varargin)

            if length(varargin) == 2
                n = varargin{1};
                total = varargin{2};
                fprintf([ptext,': %.d / %.d \n'],n,total);
            elseif length(varargin) == 4
                n1 = varargin{1};
                total1 = varargin{2};
                n2 = varargin{3};
                total2 = varargin{4};
                fprintf([ptext,': %.d / %.d ; %.d / %.d \n'],n1,total1,n2,total2);
            else
                error('Incorrect number of arguments (2 or 4 after text)')

            end % if length varargin


        end % progress_dump

        function intout = cylindrical_integration(r,z,data,varargin)
            if nargin == 3
                integral_type = 'sum';
            else
                integral_type = varargin{1};
            end % if nargin

            switch integral_type
                case 'sum'
                    dr = r(2) - r(1); dz = z(2) - z(1);
                    intout = dz*sum(2*pi*dr*sum((r').*data));
                case 'trapz'
                    intout = trapz(z,2*pi*trapz(r',(r').*data));
                case 'simpsons'
                    intout = simpsons(z,2*pi*simpsons(r',(r').*data));
            end % switch integral type
        end % cylindrical_integration

        function intout = radial_integration(r,z,data,varargin)
            if nargin == 3
                integral_type = 'sum';
            else
                integral_type = varargin{1};
            end % if nargin

            switch integral_type
                case 'sum'
                    dr = r(2) - r(1); dz = z(2) - z(1);
                    intout = dz*sum(2*pi*dr*((r').*data),2);
                case 'trapz'
                    intout = trapz(z,2*pi*trapz(r',(r').*data));
                case 'simpsons'
                    intout = simpsons(z,2*pi*simpsons(r',(r').*data));
            end % switch integral type
        end % cylindrical_integration

        function intout = cylindrical_radial_integration(r,data,varargin)
            if nargin == 2
                integral_type = 'sum';
            else
                integral_type = varargin{1};
            end % if nargin

            switch integral_type
                case 'sum'
                    dr = r(2) - r(1);
                    intout = 2*pi*dr*sum((r').*data);
                case 'trapz'
                    intout = 2*pi*trapz(r',(r').*data);
                case 'simpsons'
                    intout = 2*pi*simpsons(r',(r').*data);
                case 'just_sum'
                    dr = r(2) - r(1);
                    intout = dr*sum(data);
            end % switch integral type
        end % cylindrical_integration

        function new_r = charge_pusher(OD,prop_distance)

            norm_time = OD.norm_distance(prop_distance);

            vel = sqrt(2)*OD.npr_raw./(OD.nE_raw + 1);
            %             vel = sqrt(OD.npr_raw.^2 + OD.npth_raw.^2)./(OD.nE_raw + 1);
            new_r = abs(OD.nr_raw + vel.*norm_time);

        end

    end % static methods

end % classdef