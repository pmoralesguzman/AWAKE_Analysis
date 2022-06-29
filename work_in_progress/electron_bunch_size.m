
%________________________________________________________________________
% Plot the 2D wakefields together with the proton bunch, or each
% individually
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress 
%
% P. I. Morales Guzman
% Last update: 16/03/2021
%________________________________________________________________________

% close all;
clear;

% file location variables
% datadir = 'DW_lcode_x5_pi_rz2_th_m2_w3'; %''; DW_lcode_justelectron DW_justelectron
% extradatadir = 'DW_lcode_x5_pi_rz2_th_m2'; %'DWdc3_lcode_x1_pi'; DW_lcode_nofront

datadir = 'dr_26'; %''; DW_lcode_justelectron DW_justelectron
extradatadir = 'dr_26_o'; %'DWdc3_lcode_x1_pi'; DW_lcode_nofront

dataformat = 'mat';
useAvg = 0;
dump_list = 0:5:200;
% dump_list = 15:1:15;

% saving plot
save_flag = 1;
save_format = {'png'}; % eps, fig

% plasma properties
plasmaden = 2e14; % !!!!!! 2e14, 7e14 

% choose species density to plot
species = 'electron_seed';    %proton_beam, electrons,  electron_seed

% choose limits (in cm, must denormalize)
trans_range = [0 0.11];
xi_range = [30 0];

% choose between normalized and denormalized units
denormalize_flag = 1; % true, false

% directory to save the plots
plots_dir = ['test/fielden/',datadir,'x'];


O1 = OsirisDenormalizer('datadir',datadir,'dataformat',dataformat,...
    'useAvg',useAvg,...
    'property','raw','raw_dataset','x','direction','r',...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'species',species);

O2 = OsirisDenormalizer('datadir',extradatadir,'dataformat',dataformat,...
    'useAvg',useAvg,...
    'property','raw','raw_dataset','x','direction','r',...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'species',species);

xx = 1;
for n = 0:5:90
    O1.dump = n;
    O2.dump = n;

    O1.getdata(); O1.assign_raw();
    O2.getdata(); O2.assign_raw();
    O2.raw_dataset = 'q'; 
    O2.getdata(); O2.assign_raw();
    O2.raw_dataset = 'x';

    r_raw1 = O1.nr_raw;
    r_raw2 = O2.nr_raw;

    r1(xx) = mean(O1.nr_raw);
    r2(xx) = sum(O2.nr_raw.*O2.q_raw)./sum(O2.q_raw);
    xx = xx + 1;
    O1.progress_dump('radius',n,100)


end

plot(r1)
hold on
plot(r2)
hold off

legend('lcode','osiris')

ylabel('electron bunch size')
xlabel('z')




