% Calculate seed wakefields from a bunch distribution
% Compare initial seed wakefields for a cosine profile and for a Gaussian
% cut at different positions

clear;
close all;

plasmaden = 2e14;

savedir = 'save_files/test/MAT/FLD/';
savedir_e1 = [savedir,'e1-savg/'];
savedir_e2 = [savedir,'e2-savg/'];
savedir_b3 = [savedir,'b3-savg/'];

% Simulation grid parameters
Pt = Plotty('plasmaden',plasmaden,'save_flag',0);

nx1 = 17300/4;
nx2 = 600/4;

L1 = 173;
L2 = 6;

% full size
xmax = 173;


xi_range = Pt.denorm_distance([L1,0]);
trans_range = Pt.denorm_distance([0,L2]);

Ps = Plotty('plasmaden',plasmaden,'save_flag',0,...
    'dataformat','mat',...
    'xi_range',xi_range,'trans_range',trans_range);


dx1 = L1/nx1;
dx2 = L2/nx2;

xi = xmax - L1 + (0:dx1:(L1-dx1));
r = dx2:dx2:L2;

[XI,R] = meshgrid(xi,r);

% Physical parameters

sigzs = 203982; % sigma z as in input.f
sigr = 0.566616; % sigma r as in input.f
xis = xmax - 2.0;
xicg1 = 184.215; % center


nbn0 = 0.0001; %0.0198418;


ind1 = find(xi>xis,1);
% LC(:,ind1:nx1) = 0.0;

% Long Gaussian bunch
LG1 = exp(-(XI-xicg1).^2./sigzs).*exp(-R.^2./sigr);
LG1(:,ind1:nx1) = 0.0;

besselk0_access = besselk(0,r);
besselk1_access = besselk(1,r);

besseli0_access = besseli(0,r);
besseli1_access = besseli(1,r);

r0Lb_a = (r.').*LG1.*besseli0_access.';
rinfLb_a = (r.').*LG1.*besselk0_access.';

% Calculate fields

wr_th_LG1 = zeros(nx2,nx1);
wz_th_LG1 = zeros(nx2,nx1);

xi_parcos = xi;
xi_parsin = xi;

mkdir counter


% for ii = 1:nx1
%
%     cos_a(ii,:) = cos(xi(ii) - [1:ii,xi_parcos(ii:nx1)]);
%     sin_a(ii,:) = sin(xi(ii) - [1:ii,xi_parsin(ii:nx1)]);
%
% end

for ii = 1:nx1

    cos_a = cos(xi(ii) - xi_parcos(ii:nx1));
    sin_a = sin(xi(ii) - xi_parsin(ii:nx1));

    for j = 1:nx2
        ri = j;
        xii = ii;

        wz_th_LG1(j,ii) = sum(cos_a.*(besselk0_access(ri)* ...
            sum(r0Lb_a(1:ri,xii:nx1),1) + ...
            besseli0_access(ri)*sum(rinfLb_a(ri:nx2,ii:nx1),1)),2);

        wr_th_LG1(j,ii) = sum(sin_a.*(besselk1_access(ri)*sum(...
            r0Lb_a(1:ri,xii:nx1),1) - besseli1_access(ri)*sum(...
            rinfLb_a(ri:nx2,ii:nx1),1)),2);

    end
    fclose(fopen(['counter/',num2str(ii)], 'w'));
    Pt.progress_dump('progress ',length(dir('counter'))-2,nx1);
end

pause(2)
fclose('all');
pause(1)
rmdir('counter','s')

wr_th_LG1 = nbn0*dx1*dx2*wr_th_LG1;
wz_th_LG1 = -nbn0*dx1*dx2*wz_th_LG1;


% figure(3)
wr = Pt.denorm_Efield(wr_th_LG1);
rr = Pt.denorm_distance(r);
xxi = Pt.denorm_distance(xi);
Pt.plot_field_density('field_plot',fliplr(wr),'r_plot',[min(rr),max(rr)],'z_plot',[min(xxi),max(xxi)]);
title('Wr')

% figure(4)
wz = Pt.denorm_Efield(wz_th_LG1);
Pt.plot_field_density('field_plot',fliplr(wz),'r_plot',[min(rr),max(rr)],'z_plot',[min(xxi),max(xxi)]);
title('Wz')

if false

    % figure(5)
    Ps.datadir = 'TFHb'; Ps.dump_list = 1; Ps.useAvg = 1;
    Ps.property = 'fields'; Ps.dataformat = 'h5';
    Ps.wakefields_direction = 'long';
    Ps.property_plot = 'wakefields';
    Ps.plot_field_density();

    % figure(6)
    Ps.wakefields_direction = 'trans';
    Ps.plot_field_density();

    x1_axis = xi; x2_axis = r;
    fields = wz_th_LG1;
    time = 0;
    axis1 = [min(xi),max(xi)];
    axis2 = [min(r),max(r)];

    if ~isfolder(savedir_e1)
        mkdir(savedir_e1)
    end
    save([savedir_e1,'e1-savg-000000.mat'],'axis1','axis2','fields','time','x1_axis','x2_axis')

    fields = wr_th_LG1;
    if ~isfolder(savedir_e2)
        mkdir(savedir_e2)
    end
    save([savedir_e2,'e2-savg-000000.mat'],'axis1','axis2','fields','time','x1_axis','x2_axis')

    fields = zeros(size(wr_th_LG1));
    if ~isfolder(savedir_b3)
        mkdir(savedir_b3)
    end
    save([savedir_b3,'b3-savg-000000.mat'],'axis1','axis2','fields','time','x1_axis','x2_axis')


end
