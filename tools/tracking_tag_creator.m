%__________________________________________________________________________
% Example for the classes in the AWAKE Osiris Analysis Matlab Package
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 04/06/2020
%__________________________________________________________________________


datadir = 'gp10';
plasmaden = 1.81e14;
dump_list = 133:1:133;

dataformat = 'h5';
useAvg = false;


trans_range = [0 0.16]; 
xi_range    = [21 0];
ntag = 50000; % number of tags (macroparticles) to track

sz = 4; % scatter size

save_filename = ['protons',datadir,'test'];% save filename

OPT = OsirisParticleTracking('datadir',datadir,'plasmaden',plasmaden,...
    'dataformat',dataformat,...
    'useAvg',useAvg,...
    'trans_range',trans_range,'xi_range',xi_range,...
    'property','raw','ntag',ntag);

% initialize tags variables
tags = zeros(ntag*(length(dump_list)),2);
temp_tags = cell(length(dump_list),1);
for n = 1:length(dump_list)
    
    OPT.dump = dump_list(n);
    OPT.select_tags();
    temp_tags{n} = OPT.selected_tags';
    tags((n-1)*ntag + (1:ntag),:) = temp_tags{n};
    OPT.progress_dump('getting tags',n,length(dump_list));
    
end

tags(tags(:,1) == 0,:) = [];

tags = unique(tags,'rows');

number_of_tags = length(tags);

csvwrite([save_filename,'.tags'],[number_of_tags,0;tags])


%% check by plotting

for n = 1:length(dump_list)
    OPT.in_taglist = temp_tags{n};
    OPT.dump = dump_list(n);
    OPT.raw_dataset = 'x'; OPT.direction = 'r';
    OPT.getdata(); OPT.assign_raw();
    OPT.direction = 'z'; OPT.getdata(); OPT.assign_raw();
    
    OPT.find_tags();
    
    
    r = OPT.denorm_distance(OPT.nr_raw(OPT.ind_tag));
    z = OPT.denorm_distance(OPT.nz_raw(OPT.ind_tag));
    
    hold on 
    scatter(z,r,sz,'k','filled')
    hold off  
    xlabel('z (cm)');
    ylabel('r (cm)');
    
    OPT.progress_dump('plotting tags',n,length(dump_list));
end

