function [] = lcode2matel(datadir,propertiess,species,dumplist,lcodedumpfactor,plasmaden)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
datadirs=datadir;
properties=propertiess;
speciess=species;
dump_list=dumplist;
lcode_dump_factor=lcodedumpfactor;


%datadirs = {'fdr_26e'};

%properties = {'raw'}; % 'raw','density','fields'
%speciess = {'proton_beam'};
%speciess = {'electron_seed'};
species_list = {''}; % species list to get the raw data, empty means use only the one indicated in speciess

%dump_list = [0:1:200];

%lcode_dump_factor = 1;


O = OsirisDenormalizer('plasmaden',plasmaden,'lcode_dump_factor',lcode_dump_factor,...
    'species_list',species_list);


for ii = 1:length(datadirs)
    % specify the name of the directories where the data is in cell form
    % ex.: {'AW_ionizationfront18HD','AW_onlyproton','AW_apseed_HD256','AW_apseedx5_HD256','AW_apseed20_HD256'}
    
    O.datadir = datadirs{ii};
    
    for jj = 1:length(properties)
        
        O.property = properties{jj};
        switch O.property
            case 'density'
                for kk = 1:length(speciess)
                    O.species = speciess{kk};
                    
                    for n = 1:length(dump_list)
                        O.dump = dump_list(n);
                        O.lcode2mat();
                        O.progress_dump('lcode2mat density',n,length(dump_list))
                    end
                end
            case 'fields'
                for kk = 1:3
                    switch kk
                        case 1; O.field = 'e'; O.direction = 'z';
                        case 2; O.field = 'e'; O.direction = 'r';
                        case 3; O.field = 'b'; O.direction = 'f';
                    end
                    
                    for n = 1:length(dump_list)
                        O.dump = dump_list(n);
                        O.lcode2mat();
                        O.progress_dump(['lcode2mat fields ',O.field,O.direction],n,length(dump_list))
                    end % n dumps
                end % kk fields
            case 'raw'

                for n = 1:length(dump_list)
                    O.species = speciess{1};
                    O.dump= dump_list(n);
                    O.lcode2mat();
                    O.progress_dump('lcode2mat raw',n,length(dump_list))
                end % n dumps

        end
    end % length properties

end % length datadirs


end

