x = 0:0.5:2661.5*2;
r = 26.6125; % 1 cm
ramp = 0.5 - (-1*(x-2661.5)/(2*r))./sqrt((-1*(x-2661.5)/(2*r)).^2 + 0.25)/2;

P =  Plotty('plasmaden',2e14);
figure(1)
plot(x,ramp)
figure(2)
plot(P.denorm_distance(x),ramp)

point_five = num2str(0.5*ones(length(ramp)-1,1));
str_L = ' L '; 
ramp1 = num2str(ramp(1:end-1)');
ramp2 = num2str(ramp(2:end)');

str_L_mat = repmat(str_L,length(ramp)-1,1);
str_space_mat = repmat(' ',length(ramp)-1,1);

text_cfg = [point_five,str_space_mat,ramp1,str_L_mat,ramp2];

writematrix(text_cfg,'test_cfg.txt');