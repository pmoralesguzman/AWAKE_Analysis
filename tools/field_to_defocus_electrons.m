
O = OsirisLoader();

c = physconst('lightspeed'); % m/s
m = O.e_mass_eV; % MeV/c^2
me = O.e_mass_kg; % kg
K = 18.89; % MeV
l = K/m;
q = O.e_charge_C; % C

Dy = 2e-4; % m
Dz = 0.3; % m


E = (c^2*me*l*Dy)/(q*Dz^2)/1e6

