% Generates a Graph of K value vs ARL for fixed hypothesis and Actual

load('./ARK_Kvary_data_mu_lt_mu0.mat');

%debug output
ARL

% for saving jpg
figure, close;

%plot and hold for other data
plot(K,ARL);

%plot formatting and labeling
axis([0 0.5 0 10000])
leglbl = legend('ARL');
xlbl = xlabel('K value');
ylbl = ylabel('ARL');
titlelbl = title(sprintf('ARL with variable K, mu = %5.3f, mu_0 = %5.3f',mu,mu_0));

%set page resolution, assume square
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 1600 1200]);

%set font size
set(leglbl, 'fontsize', font);
set(xlbl, 'fontsize', font);
set(ylbl, 'fontsize', font);
set(titlelbl, 'fontsize', font);

%save to file (jpg)
print(gcf, '-djpeg', '../../SimulationFigures/ARL_k_vary_mu_lt_mu0_10k.jpg')
