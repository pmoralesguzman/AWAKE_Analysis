
close all
N = 10000;
x = linspace(0,pi*8,N);
Fs = 1/(x(2)-x(1));
y1 = 0*sin(1.05*x).*sin(0.05*x+pi/2);
y2 = 2*(sin(x)+sin(1.015*x));
xx = linspace(0,pi*32,16/2+1);
[xb,yb] = stairs(xx,sin(0.05*xx + pi/2));
xb(2:2:end) = xb(2:2:end) - 0.000001;
yy = interp1(xb,yb,x,'linear',0);
y3 = 0*sin(1*x).*yy;

y = y1 + y2 + y3;

figure
plot(x,y)


[p,l] = findpeaks(y,x);

hold on
scatter(l,p,'filled')
hold off

figure
scatter(1:length(l),l,'filled')

figure

plot(1:(length(l)-1),diff(l))

NN = N*1;
ff = fft(y,NN);

P2 = abs(ff/NN);
P1 = P2(1:NN/2+1);
P1(2:end-1) = 2*P1(2:end-1);
ff_th = ff(1:NN/2+1);
ff_th(ff_th<max(ff_th)/100) = 0;
theta = angle(ff_th);
f = Fs*(0:(NN/2))/NN;
figure(6)
plot(f,P1)
title('FFT of toy model')




fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0.0,0.5,0],...
               'Upper',[5,1.5,2*pi],...
               'StartPoint',[2 0.9 1],...
               'MaxFunEvals',5000,'Display','iter',...
               'DiffMinChange',1e-4,...
               'DiffMaxChange',0.25);
ft = fittype('a*sin(b*x + c)','options',fo);

[curve2,gof2] = fit(x(:),y(:),ft)

figure
hold on
plot(x,y)
plot(x,curve2.a*sin(curve2.b*x + curve2.c))
hold off









