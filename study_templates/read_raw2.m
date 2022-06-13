%%%%%%


O = OsirisDenormalizer('plasmaden',1.81e14,'datadir','gm20',...
    'property','raw','raw_dataset','q',...
    'dump',52,'dataformat','h5','species','proton_beam');

O.raw_dataset = 'p'; O.direction = 'r';

O.getdata(); O.assign_raw();

O.raw_dataset = 'x'; O.direction = 'r';

O.getdata(); O.assign_raw();

O.raw_dataset = 'ene';

O.getdata(); O.assign_raw();

O.raw_dataset = 'q';

O.getdata(); O.assign_raw();

ir = O.denorm_distance(O.nr_raw) > 0.15;


r_new = O.denorm_distance(O.charge_pusher(O,1350-520));

r_roi = r_new(ir);

histogram(r_roi)



% 
% x = 1;
% for n = 110:200
%     O.dump = n;
%     O.raw_dataset = 'q';
%     O.getdata(); O.assign_raw();
%     
%     O.raw_dataset = 'ene';
%     O.getdata(); O.assign_raw();
%     
% %     energy(x) = sum(O.nE_raw.*O.q_raw)./sum(O.q_raw)/initial_energy;
%     energy(x) = max(O.nE_raw)/initial_energy;
%     x = x + 1;
%     O.progress_dump('hola',n,200);
% end
% 
%     plot(linspace(11,20,length(energy)),energy,'LineWidth',2);
% 
%     xlabel('z');
%     ylabel('max energy / average initial energy (a.u.)');
%     title('gap = 0 m')