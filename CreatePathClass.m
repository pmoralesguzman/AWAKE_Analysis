%---------------------------------------------------------------------
% Tells MATLAB where the directories are.
% Run before any other script.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
%
% P. I. Morales Guzman
% Last update: 07/01/2021
%---------------------------------------------------------------------

function CreatePathClass(varargin)
% Create path where all studies and files are.

if nargin == 0
    force = false;
elseif nargin == 1
    force = true;
end

pathCell = regexp(path, pathsep, 'split');

paths = {'../AWAKE_Analysis',... all directories with analysis code
    'E:/AWAKE_Data',... for the external hard drive (must be set to E:/)
    'D:/AWAKE_Data_HD',... for the SSD in the work laptop
    'D:/LCODE_CERN',... for the hard drive desktop PC at CERN
    'E:/AWAKE_Data_SSD',... for the portable SSD
    'C:/AWAKE_Data_laptop',... for the SSD in the personal laptop
    '../simulations',... for Cobra
    '../lcode_simulations',... for Cobra
    '../../../../../local/scratch/pmorales/real_studies/',...  for the hard drive desktop PC at MPP
    '../../../../../../../../Volumes/SanDisk/HIWI/test_runs',... for SSD in E's personal laptop
    '../../../../../../../../Volumes/SanDisk/HIWI/AWAKE_Data',... for SSD in E's personal laptop
    '../../../../../local/scratch/elainedv'}; % for the hard drive desktop PC at MPP


for pp = 1:length(paths)
    s = what(paths{pp});
    if ispc && ~isempty(s) && ~any(strcmpi(s.path, pathCell)) || force  % Windows is not case-sensitive
        addpath(genpath(paths{pp}))
    elseif (~ispc) && ~isempty(s) && ~any(strcmp(s.path, pathCell)) || force 
        addpath(genpath(paths{pp}))
    end
end % paths

end

