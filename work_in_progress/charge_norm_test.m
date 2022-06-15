a = sum(q)*(4*pi*P.permittivity*P.denorm_distance(0.02)/100*P.e_mass_kg*P.c_m^2)...
    ./(2*P.e_charge_C^1)

% w = 2*beampop*e*c / (IA*(dz/k)*N)
% N = (IA*(dz/k)*w) / (2*beampop*e*c)

% IA = 4*pi*e0*m*c**3/e

% N = (4*pi*e0*m*c^2*(dz/k)*w) / (2*e^2)

