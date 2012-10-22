% Generates a Graph of K value vs ARL for fixed hypothesis and Actual

%Jpeg params
res = 1600;
font = 25;

%actual
mu = 0.50;

%hypothesis
mu_0 = 0.50;

%sample length
N = 10000;

%slack value
K= 0:.01:.5;
%head start
ustart=0;
lstart=0;

%threshold to graph
H = 3;

%number of trials to run
trials = 51;

%plot params
ymax= 1000;

%generate the samples compute the ARL, store in ARL variable
ARL = zeros(size(K));
for i = 1:length(K)
	data = zeros(1,trials);
	for t = 1:trials
		samples = Sample1dim(mu,N);
		cusum = Gen_tab_x(samples,mu_0,K(i),ustart,lstart);
		data(t) = Threshold(cusum, H, 0);
	end
	ARL(i) = mean(data);
end

%debug output
ARL

save('./ARK_Kvary_data.mat');

% for saving jpg
figure, close;

%plot and hold for other data
plot(K,ARL);

%plot formatting and labeling
axis([0 0.5 -ymax ymax])
leglbl = legend('ARL');
xlbl = xlabel('K value');
ylbl = ylabel('ARL');
titlelbl = title(sprintf('ARL with variable K, mu = %5.3f, mu_0 = %5.3f',mu,mu_0));

%set page resolution, assume square
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 res res]);

%set font size
set(leglbl, 'fontsize', font);
set(xlbl, 'fontsize', font);
set(ylbl, 'fontsize', font);
set(titlelbl, 'fontsize', font);

%save to file (jpg)
print(gcf, '-djpeg', '../SimulationFigures/ARL_k_vary.jpg')
