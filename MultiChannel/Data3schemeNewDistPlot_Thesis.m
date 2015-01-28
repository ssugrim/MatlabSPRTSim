load('I:\My Documents\Academics\research\ChannelSurfing\Simulation\SampleGen\MultiChannel\Data3SchemeNewDist-Nov04.2013-04.02CURRENT.mat');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 1600 1200]);
set(0,'defaultlinelinewidth',2)

f = figure(1);

plot(Channels,Data_SPRT_channels_found_H0,'-^g',Channels,Data_Tree_channels_found_H0,'-*r',Channels,Data_Simple_channels_found_H0,'-ob');
legend('greedy-SPRT, \lambda = 0.3, \delta = 0.15', 'Tree, \psi = 4, \rho = 0.5', 'Simple, L/N = 0.97');
xlim([20,1050]);
xlabel('N, number of Channels');
ylabel('Number of channels classifed as H_0');
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',2);


f = figure(2);
plot(Channels,Data_SPRT_error_2,'-^g',Channels,Data_Tree_error_2,'-*r',Channels,Data_Simple_error_2,'-ob');
legend('greedy-SPRT, \lambda = 0.3, \delta = 0.15', 'Tree, \psi = 4, \rho = 0.5', 'Simple, L/N = 0.97');
xlim([20,1050]);
xlabel('N, number of Channels');
ylabel('Probablity of missclasification as H_0');
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);

f = figure(3);
plot(Channels,Data_SPRT_CPU_time,'-^g',Channels,Data_Tree_CPU_time,'-*r',Channels,Data_Simple_CPU_time,'-ob');
xlim([20,1050])
xlabel('N, number of Channels');
ylabel('CPU time in seconds');
legend('greedy-SPRT, \lambda = 0.3, \delta = 0.15', 'Tree, \psi = 4, \rho = 0.5', 'Simple, L/N = 0.97');
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);

f = figure(4);
plot(Channels,Data_SPRT_channels_found_H0,'-^g',Channels,Data_Tree_channels_found_H0,'-*r');
legend('greedy-SPRT, \lambda = 0.3, \delta = 0.15', 'Tree, \psi = 4, \rho = 0.5');
xlim([20,1050]);
xlabel('N, number of Channels');
ylabel('Number of channels classifed as H_0');
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);

f = figure(5);
plot(Channels,Data_SPRT_error_2,'-^g',Channels,Data_Tree_error_2,'-*r');
legend('greedy-SPRT, \lambda = 0.3, \delta = 0.15', 'Tree, \psi = 4, \rho = 0.5');
xlim([20,1050]);
xlabel('N, number of Channels');
ylabel('Probablity of missclasification as H_0');
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);

f = figure(6);
plot(Channels,Data_SPRT_mean_channels_passed_on,'-^g')
xlim([20,1050])
xlabel('N, number of Channels')
ylabel('Number of Channels passed on');
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);

f = figure(7);
foo = Data_SPRT_mean_samples_to_decide - Data_SPRT_mean_samples_make_decision;
plot(Channels,Data_SPRT_mean_samples_make_decision,'--xm',Channels,foo,'-.sk')
legend('Mean number of samples to complete test','Mean number of samples to find channel');
xlabel('N, number of Channels')
ylabel('Samples');
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);



