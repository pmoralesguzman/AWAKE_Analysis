% Calculate seed wakefields from a bunch distribution
% Compare initial seed wakefields for a cosine profile and for a Gaussian
% cut at different positions

clear;
close all;

plasmaden = 2e14;

savefiles_flag = 1;
savedir = 'save_files/electron_bunch/MAT/FLD/';
savedir_e1 = [savedir,'e1-savg/'];
savedir_e2 = [savedir,'e2-savg/'];
savedir_b3 = [savedir,'b3-savg/'];

savedir_den = 'save_files/electron_bunch/MAT/DENSITY/';
savedir_p = [savedir_den,'proton_beam/charge-savg/'];

saveplot_flag = 1;
saveplot_dir = 'DW/theory';
saveplot_name = 'electronbunch';

plotsim_flag = 0;

% Simulation grid parameters
Pt = Plotty('plasmaden',plasmaden,'save_flag',saveplot_flag,...
    'plots_dir',saveplot_dir,'save_format',{'png'});

nx1 = 1500/4;
nx2 = 300/4;

L1 = 15;
L2 = 3;

% full size
xmax = 15;


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

sigzs = 0.798401; % sigma z as in input.f
sigr = 0.71856; % sigma r as in input.f
xis = xmax - 2.0;
% xiback = xmax - ;
xicg1 = 9; % center


nbn0 = -0.05; 


ind1 = find(xi>xis,1);
ind2 = 1;
% LC(:,ind1:nx1) = 0.0;

% Long Gaussian bunch
LG1 = exp(-(XI-xicg1).^2./sigzs).*exp(-R.^2./sigr);
LG1(:,ind1:nx1) = 0.0;

besselk0_access = besselk(0,r);
besselk1_access = besselk(1,r);

besseli0_access = besseli(0,r);
besseli1_access = besseli(1,r);

r0Lb_a = (r.').*LG1.*besseli0_access.';
r0Lb_a_wz = r0Lb_a;
rinfLb_a = (r.').*LG1.*besselk0_access.';

% Calculate fields

wr_th_LG1 = zeros(nx2/2,nx1);
wz_th_LG1 = zeros(nx2/2,nx1);
% wr_th_LG2 = zeros(nx2/2,nx1);
% wz_th_LG2 = zeros(nx2/2,nx1);

mkdir counter

cos_a = zeros(nx1,nx1);
sin_a = zeros(nx1,nx1);

for ii = 1:nx1
    
    cos_a(ii,:) = cos(xi(ii) - [(xi(ii)-pi/2)*ones(1,ii-1),xi(ii:nx1)]); %.*([zeros(1,ii-1),ones(1,nx1-ii+1)]);
    sin_a(ii,:) = sin(xi(ii) - [xi(ii)*ones(1,ii-1),xi(ii:nx1)]); %.*([zeros(1,ii-1),ones(1,nx1-ii+1)]);

end

% r0Lb_sum = cumsum(r0Lb_a);
r0Lb_sum_wz = (besselk0_access').*cumsum(r0Lb_a);
r0Lb_sum_wr = (besselk1_access').*cumsum(r0Lb_a);

% r0Lb_sum_wz2 = r0Lb_sum_wz(nx2/2:end,:);
% r0Lb_sum_wr2 = r0Lb_sum_wr(nx2/2:end,:);

% rinfLb_sum = cumsum(rinfLb_a,'reverse');
rinfLb_sum_wz = (besseli0_access').*cumsum(rinfLb_a,'reverse');
rinfLb_sum_wr = (besseli1_access').*cumsum(rinfLb_a,'reverse');

% rinfLb_sum_wz2 = rinfLb_sum_wz(nx2/2:end,:);
% rinfLb_sum_wr2 = rinfLb_sum_wr(nx2/2:end,:);

counter_fun = @Pt.progress_dump;

for j = 1:nx2
    wz_th_LG1(j,:) = sum(cos_a.* ...
        (r0Lb_sum_wz(j,:) + ...
        rinfLb_sum_wz(j,:)),2)';

    wr_th_LG1(j,:) = sum(sin_a.* ...
        (r0Lb_sum_wr(j,:) - ...
        rinfLb_sum_wr(j,:)),2)';

%     wz_th_LG2(j,:) = sum(cos_a.* ...
%         (r0Lb_sum_wz2(j,:) + ...
%         rinfLb_sum_wz2(j,:)),2)';
% 
%     wr_th_LG2(j,:) = sum(sin_a.* ...
%         (r0Lb_sum_wr2(j,:) - ...
%         rinfLb_sum_wr2(j,:)),2)';

    fclose(fopen(['counter/',num2str(j)], 'w'));
    counter_fun('progress ',length(dir('counter'))-2,nx2);

end

pause(2)
fclose('all');
pause(1)
rmdir('counter','s')

wr_th_LG = nbn0*dx1*dx2*wr_th_LG1;
wz_th_LG = -nbn0*dx1*dx2*wz_th_LG1; 

% wr_th_LG = nbn0*dx1*dx2*[wr_th_LG1;wr_th_LG2];
% wz_th_LG = -nbn0*dx1*dx2*[wz_th_LG1;wz_th_LG2]; 


% figure(3)
wr = Pt.denorm_Efield(wr_th_LG);

rr = Pt.denorm_distance(r);
xxi = Pt.denorm_distance(xi);

Pt.plot_field_density('density_plot',fliplr(LG1),'field_plot',fliplr(wr),'r_plot',[min(rr),max(rr)],'z_plot',[min(xxi),max(xxi)]);
title('Wr');

Pt.plot_name = [saveplot_name,'long'];
Pt.save_plot();

% figure(4)
wz = Pt.denorm_Efield(wz_th_LG);
Pt.plot_field_density('density_plot',fliplr(LG1),'field_plot',fliplr(wz),'r_plot',[min(rr),max(rr)],'z_plot',[min(xxi),max(xxi)]);
title('Wz');

Pt.plot_name = [saveplot_name,'trans'];
Pt.save_plot();

if plotsim_flag

figure(5)
Ps.datadir = 'DWe'; Ps.dump_list = 1; Ps.useAvg = 1;
Ps.property = 'fields'; Ps.dataformat = 'h5';
Ps.wakefields_direction = 'long';
Ps.property_plot = 'wakefields';
Ps.plot_field_density();

figure(6)
Ps.wakefields_direction = 'trans';
Ps.plot_field_density();

end

if savefiles_flag

x1_axis = xi; x2_axis = r;
fields = wz_th_LG;
time = 0;
axis1 = [min(xi),max(xi)];
axis2 = [min(r),max(r)];

if ~isfolder(savedir_e1)
    mkdir(savedir_e1)
end
save([savedir_e1,'e1-savg-000000.mat'],'axis1','axis2','fields','time','x1_axis','x2_axis')

fields = wr_th_LG;
if ~isfolder(savedir_e2)
    mkdir(savedir_e2)
end
save([savedir_e2,'e2-savg-000000.mat'],'axis1','axis2','fields','time','x1_axis','x2_axis')

fields = zeros(size(wr_th_LG));
if ~isfolder(savedir_b3)
    mkdir(savedir_b3)
end
save([savedir_b3,'b3-savg-000000.mat'],'axis1','axis2','fields','time','x1_axis','x2_axis')

density = LG1;
if ~isfolder(savedir_p)
    mkdir(savedir_p)
end
save([savedir_p,'charge-savg-proton_beam-000000.mat'],'axis1','axis2','density','time','x1_axis','x2_axis')


end
