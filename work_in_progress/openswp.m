a = fopen('C:\LCODE_HP\DWdc3_lcode\ez26124.swp','r');
% a = fopen('/home/iwsatlas1/pmorales/LCODE_MPP/DWb_lcode/ez01000.swp','r');
%c = fscanf(a,'%f');
b = fscanf(a,'%f',[301,30000]);

imagesc(fliplr(b))

colormap(bluewhitered)

% b = fscanf(a,'%f');

% P = Plotty();

nx1 = 30000;
nx2 = 301;

L1 = 300;
L2 = 3;

% full size
xmax = 300;

plasmaden = 2e14;

Ps = Plotty('plasmaden',plasmaden,'save_flag',0,...
    'dataformat','mat');


dx1 = L1/nx1;
dx2 = L2/nx2;

xi = xmax - L1 + (0:dx1:(L1-dx1));
r = dx2:dx2:L2;



savedir = 'save_files/DWb_lcode/MAT/FLD/';
savedir_e1 = [savedir,'e1-savg/'];
savedir_e2 = [savedir,'e2-savg/'];
savedir_b3 = [savedir,'b3-savg/'];

%% save data

x1_axis = xi; x2_axis = r;
fields = fliplr(b);
time = 0;
axis1 = [min(xi),max(xi)];
axis2 = [min(r),max(r)];

if ~isfolder(savedir_e1)
    mkdir(savedir_e1)
end
save([savedir_e1,'e1-savg-000000.mat'],'axis1','axis2','fields','time','x1_axis','x2_axis')
