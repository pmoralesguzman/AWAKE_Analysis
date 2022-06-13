

z = 0:10;

g = 2;

n = 1.81e14*(1 + g/100*z);


n2 = 1.81e14*(1 + 0/100*z);

hold on
plot(z,n,'k')
plot(z,n2,'--r')
hold off


xlabel('z (m)')
ylabel('plasma density (cm^{-3})')

ylim([0 2.5e14])
xlim([0 10])