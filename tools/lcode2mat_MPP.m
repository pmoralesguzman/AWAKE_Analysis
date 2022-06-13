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
datadirs = {'etest2'};
% properties = {'fields','raw','density'}; % raw
properties = {'fields','raw'}; % raw
speciess = {'proton_beam'};
field_k = 3; % 1 up to ez, 2 up to er, 3 up to bth (all)

dump_list = (0:1:100); % 4:2:100
lcode_dump_factor = 1;


O = OsirisLoader('plasmaden',2e14,'lcode_dump_factor',lcode_dump_factor);


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
                for kk = 1:field_k
                    switch kk
                        case 1; O.field = 'e'; O.direction = 'z';
                        case 2; O.field = 'e'; O.direction = 'r';
                        case 3; O.field = 'b'; O.direction = 'f';
                    end
                    
                    for n = 1:length(dump_list)
                        O.dump = dump_list(n);
                        O.lcode2mat();
                        O.progress_dump('lcode2mat fields',n,length(dump_list))
                    end % n dumps
                end % kk fields
            case 'raw'
                for ll = 1:length(speciess)
                    O.species = speciess{ll};
                    for kk = 1:1
                        switch kk
                            case 1; O.raw_dataset = 'q';
                            case 2; O.raw_dataset = 'x'; O.direction = 'z';
                            case 3; O.raw_dataset = 'p'; O.direction = 'z';
                            case 8; O.raw_dataset = 'tag';
                            case 5; O.raw_dataset = 'ene';
                            case 6; O.raw_dataset = 'x'; O.direction = 'r';
                            case 7; O.raw_dataset = 'p'; O.direction = 'r';
                            case 4; O.raw_dataset = 'p'; O.direction = '\theta';
                        end
                        
                        for n = 1:length(dump_list)
                            O.dump= dump_list(n);
                            O.lcode2mat();
                            O.progress_dump('lcode2mat raw',n,length(dump_list))
                        end % n dumps
                    end %% ll speciess
                end % kk raw
        end
    end % length properties
    
end % length datadirs

