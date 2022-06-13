close all

load('density_vs_phasediff.mat')

densities_plot = densities(1:2:end);
densities_plot = [0,1,1.5,2]';


% phasediff_plot = (abs(diff_list(1:2:end)-180)/180 + abs(diff_list(2:2:end)-360*0.75)/(360*0.75))/2; 

phasediff_plot_1([1,3,5,6]) = abs(diff_list([1,3,5,6])-180)/180;
phasediff_plot_1([4,7]) = abs(diff_list([4,7])-270)/270;
phasediff_plot_1([2]) = abs(diff_list(2)-90)/90;

phasediff_plot(4) = mean(phasediff_plot_1([6,7]));
phasediff_plot(1) = mean(phasediff_plot_1([1,2])); 
phasediff_plot(2) = mean(phasediff_plot_1([3,4])); 
phasediff_plot(3) = mean(phasediff_plot_1(5)); 
 

fig1 = figure(1);

scatter(densities_plot(1:4),phasediff_plot(1:4),'filled')

xlabel('distance from density cut ($\sigma_z$)','Interpreter','latex');
ylabel('average phase difference','Interpreter','latex');

P = Plotty('save_format',{'png'});
P.fig_handle = fig1;
P.plot_name = 'density_vs_phasediff';
P.plots_dir = 'DW_PEB';

P.save_plot();

