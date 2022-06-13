%________________________________________________________________________
% Script to calculate, from the waterfall data, the dephasing of the fields
% close to some selected xi.
% Especially developed for the APS plots.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/07/2020
%________________________________________________________________________

clear;
close all;
datadirs = {'DWdc3_lcode_x2','DWdc3_lcode_x4','DWdc3_lcode_x5','DWdc3_lcode_x10'};

plasmaden = 2e14;
property = 'fields';
dump_list = 0:2:98;
useAvg = 0;
dataformat = 'mat';
save_format = {'png'};

trans_range = [0 0.02]; % if density
% dephasing_xis = [1,1.5,2,3,4,5,6,7]; % cm
dephasing_xis = [7.5];
dephasing_search = '0x'; % 0x, max
force_waterfall = true;


for ph = 1:length(dephasing_xis)
    dephasing_xi = dephasing_xis(ph);
    close all;
    for d = 1:length(datadirs)
        
        datadir = datadirs{d};
        
        
        OPA = OsirisPhaseAnalysis('datadir',datadir,...
            'property',property,'species','proton_beam',...
            'direction','z',...
            'wakefields_direction','long',...
            'trans_range',trans_range,...
            'plasmaden',plasmaden,...
            'dump_list',dump_list,'useAvg',useAvg,...
            'dataformat',dataformat,...
            'dephasing_xi',dephasing_xi,...
            'dephasing_search',dephasing_search,...
            'force_waterfall',force_waterfall,...
            'lineout_point',3);
        OPA.dephasing();
        
    end
    
end



