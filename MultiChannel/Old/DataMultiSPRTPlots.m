Ratios = Channels ;

figure(1); hold;


plot(Ratios,mean(squeeze(Data_SPRT_error_2(5,:,:))),'--b');
plot(Ratios,mean(squeeze(Data_SPRT_error_2(4,:,:))),'.-r');
plot(Ratios,mean(squeeze(Data_SPRT_error_2(3,:,:))),'x-g');
plot(Ratios,mean(squeeze(Data_SPRT_error_2(2,:,:))),'o-m');
plot(Ratios,mean(squeeze(Data_SPRT_error_2(1,:,:))),'*-k');

xlim([min(Ratios),max(Ratios)]);
legend('\lambda = 0.15','\lambda = 0.1250','\lambda = 0.1','\lambda = 0.0750','\lambda = 0.005')
xlabel('Channels Searched');
ylabel('Type 2 error');
title('Type 2 error as a function of Channels Searched');

set(0,'DefaultAxesFontSize',20);
set(gca, 'FontSize', 20);
set(gcf, 'Position', [0 0 1200 1200]);
saveas(gcf, './MultiLambdaType2.png');


figure(2); hold;

plot(Ratios,mean(squeeze(Data_SPRT_error_1(5,:,:))),'--b');
plot(Ratios,mean(squeeze(Data_SPRT_error_1(4,:,:))),'.-r');
plot(Ratios,mean(squeeze(Data_SPRT_error_1(3,:,:))),'x-g');
plot(Ratios,mean(squeeze(Data_SPRT_error_1(2,:,:))),'o-m');
plot(Ratios,mean(squeeze(Data_SPRT_error_1(1,:,:))),'*-k');

xlim([min(Ratios),max(Ratios)]);
legend('\lambda = 0.15','\lambda = 0.1250','\lambda = 0.1','\lambda = 0.0750','\lambda = 0.005')
xlabel('Channels Searched');
ylabel('Type 1 error');
title('Type 1 error as a function of Channels Searched');

set(0,'DefaultAxesFontSize',20);
set(gca, 'FontSize', 20);
set(gcf, 'Position', [0 0 1200 1200]);
saveas(gcf, './MultiLambdaType1.png');

figure(3); hold;

plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_1(5,:,:))),'--b');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_1(4,:,:))),'.-r');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_1(3,:,:))),'x-g');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_1(2,:,:))),'o-m');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_1(1,:,:))),'*-k');

xlim([min(Ratios),max(Ratios)]);
legend('\lambda = 0.15','\lambda = 0.1250','\lambda = 0.1','\lambda = 0.0750','\lambda = 0.005')
xlabel('Channels Searched');
ylabel('Mean of the estimates placed in bin 1');
title('Mean of parameters place in bin 1 as a function of Channels Searched');

set(0,'DefaultAxesFontSize',20);
set(gca, 'FontSize', 20);
set(gcf, 'Position', [0 0 1200 1200]);
saveas(gcf, './MultiLambdaMeanBin1.png');

figure(4); hold;

plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_2(5,:,:))),'--b');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_2(4,:,:))),'.-r');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_2(3,:,:))),'x-g');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_2(2,:,:))),'o-m');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_2(1,:,:))),'*-k');

xlim([min(Ratios),max(Ratios)]);
legend('\lambda = 0.15','\lambda = 0.1250','\lambda = 0.1','\lambda = 0.0750','\lambda = 0.005')
xlabel('Channels Searched');
ylabel('Mean of the estimates placed in bin 2');
title('Mean of parameters place in bin 2 as a function of Channels Searched');

set(0,'DefaultAxesFontSize',20);
set(gca, 'FontSize', 20);
set(gcf, 'Position', [0 0 1200 1200]);
saveas(gcf, './MultiLambdaMeanBin2.png');

figure(5); hold;

plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_3(5,:,:))),'--b');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_3(4,:,:))),'.-r');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_3(3,:,:))),'x-g');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_3(2,:,:))),'o-m');
plot(Ratios,mean(squeeze(Data_SPRT_mean_bin_3(1,:,:))),'*-k');

xlim([min(Ratios),max(Ratios)]);
legend('\lambda = 0.15','\lambda = 0.1250','\lambda = 0.1','\lambda = 0.0750','\lambda = 0.005');
xlabel('Channels Searched');
ylabel('Mean of the estimates placed in bin 3');
title('Mean of parameters place in bin 3 as a function of Channels Searched');

set(0,'DefaultAxesFontSize',20);
set(gca, 'FontSize', 20);
set(gcf, 'Position', [0 0 1200 1200]);
saveas(gcf, './MultiLambdaMeanBin3.png');

figure(6); hold;

plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_between_decsion(5,:,:))),'--b');
plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_between_decsion(4,:,:))),'.-r');
plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_between_decsion(3,:,:))),'x-g');
plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_between_decsion(2,:,:))),'o-m');
plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_between_decsion(1,:,:))),'*-k');

xlim([min(Ratios),max(Ratios)]);
legend('\lambda = 0.15','\lambda = 0.1250','\lambda = 0.1','\lambda = 0.0750','\lambda = 0.005')
xlabel('Channels Searched');
ylabel('Mean Samples takes between decisions');
title('How many samples on average are spent finding one channel.');

set(0,'DefaultAxesFontSize',20);
set(gca, 'FontSize', 20);
set(gcf, 'Position', [0 0 1200 1200]);
saveas(gcf, './MultiLambdaMeanSamplesBetween.png');

figure(7); hold;

plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_make_decision(5,:,:))),'--b');
plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_make_decision(4,:,:))),'.-r');
plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_make_decision(3,:,:))),'x-g');
plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_make_decision(2,:,:))),'o-m');
plot(Ratios,mean(squeeze(Data_SPRT_mean_samples_make_decision(1,:,:))),'*-k');

xlim([min(Ratios),max(Ratios)]);
legend('\lambda = 0.15','\lambda = 0.1250','\lambda = 0.1','\lambda = 0.0750','\lambda = 0.005');
xlabel('Channels Searched');
ylabel('Mean Samples to complete a Channel');
title('Of the samples required to find a channel, how many are actually spent classifying it.');

set(0,'DefaultAxesFontSize',20);
set(gca, 'FontSize', 20);
set(gcf, 'Position', [0 0 1200 1200]);
saveas(gcf, './MultiLambdaMeanSamplesNeededToClassify.png');

figure(8); hold;

plot(Ratios,mean(squeeze(Data_SPRT_percent_classified(5,:,:)) .* 100),'--b');
plot(Ratios,mean(squeeze(Data_SPRT_percent_classified(4,:,:)) .* 100),'.-r');
plot(Ratios,mean(squeeze(Data_SPRT_percent_classified(3,:,:)) .* 100),'x-g');
plot(Ratios,mean(squeeze(Data_SPRT_percent_classified(2,:,:)) .* 100),'o-m');
plot(Ratios,mean(squeeze(Data_SPRT_percent_classified(1,:,:)) .* 100),'*-k');

xlim([min(Ratios),max(Ratios)]);
legend('\lambda = 0.15','\lambda = 0.1250','\lambda = 0.1','\lambda = 0.0750','\lambda = 0.005');
xlabel('Channels Searched');
ylabel('Percentage of channels Calssified');
title('Percentage of channels classifed vs Channels searched');

set(0,'DefaultAxesFontSize',20);
set(gca, 'FontSize', 20);
set(gcf, 'Position', [0 0 1200 1200]);
saveas(gcf, './MultiLambdaMeanSamplesPercentClassifed.png');



