%________________________________________________________________________
% Script to convert lcode output files into AWAKE_Class readable files.
%
% LCODE
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 22/11/2021
%________________________________________________________________________

% input
%location = 'C:/LCODE_MPP'; % C:/LCODE_HP, E:/LCODE_MPP, /home/iwsatlas1/pmorales/LCODE_MPP/
% datadirs = {'DW_lcode_x5_t4','DW_lcode_x5_r4','DW_lcode_x5_z4','DW_lcode_x5_m4'};
datadirs = {'testcu'};

properties = {'raw','fields','density'}; % raw 'fields'
% properties = {'fields','raw','density'}; % raw 'fields'
speciess = {'electron_seed', 'proton_beam'};  %proton_beam
% speciess = {'electron_seed'};

dump_list = [0:200];

lcode_dump_factor = 1;


O = OsirisDenormalizer('plasmaden',2e14,'lcode_dump_factor',lcode_dump_factor);


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

