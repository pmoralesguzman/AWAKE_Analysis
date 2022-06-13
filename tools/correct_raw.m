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
datadir = 'r2l_202_c_noe';
dump_list = 0:1:0;
dataformat = 'mat';
species_name = 'proton_beam';

% save parameters
save_format = {'png'};
save_flag = 1;
save_plot_flag = 0;

% load classes
O = OsirisDenormalizer('plasmaden',2e14,'datadir',datadir,...
    'property','fields',...
    'dataformat',dataformat,'species',species_name);
O.getdata();

O.property = 'raw';

% normalize to charge
% analytical charge fraction

for n = 1:length(dump_list)
    close all
    O.dump = dump_list(n);

    % get data
    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();

    partialpath_handle = what(O.datadir);
    save_partialpath = [partialpath_handle.path,'/MAT/'];
    % Names of files (with correct numbering)
    dump_char = sprintf('%06.6d',O.dump);
    partialpath = [save_partialpath,'RAW/','proton_beam','/','q','/'];

    filename = ['RAW','-','proton_beam','-','q'];
    fullpath = [partialpath,filename,'-',dump_char,'.mat'];
    rawdata = O.q_raw*1836.15267343;
    time = O.time;
    save(fullpath,'rawdata','time'); %'-v6'


    O.progress_dump('correcting q: ',n,length(dump_list));

end % for dump