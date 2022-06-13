fileID = fopen('../LCODE_Sims/awake-baseline/beamfile.bin');
fileID = fopen('../LCODE_Sims/equidistant-bunch-train/beamfile.bin');
A = fread(fileID,[8 2000],'double');