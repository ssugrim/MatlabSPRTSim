% Generates data and plots the CUSUM under given parameters

%plotting params
res = 1600;
font = 25;


mu = 0.89;
mu_0 = 0.95;
%sample length
N = 10000;
%slack value
K=.2;
%head start
ustart=0;
lstart=0;

%graph max
Y_max = 5;
X_max = 100;

%threshold to graph
H = 3;

%generate the samples and CUSUM via my functions
samples = Sample1dim(mu,N);
cusum = Gen_tab_x(samples,mu_0,K,ustart,lstart);

% for saving jpg
figure, close;

%plot and hold for other data
plot(cusum');
hold on;

%threhsold lines
topline = H * ones(1,N);
botline = -H * ones(1,N);
plot(0:N-1,topline,'red',0:N-1,botline,'black');

%plot formatting and labeling
leglbl = legend('C^{+}_{i}','C^{+}_{i}','y=3','y=-3');
axis([0 X_max -Y_max Y_max])
xlbl = xlabel('step index, i');
ylbl = ylabel('C_i');
mu_0_str = sprintf('%5.2f',mu_0);
mu_str = sprintf('%5.2f',mu);
n_str = sprintf(' N = %6.0f,',N);
k_str = sprintf(' K = %5.2f',K);
titlelbl = title(strcat('Tabular CUSUM for \mu_0 =',mu_0_str,' and \mu =',mu_str,'; with parameters:',n_str,k_str));
hold off;


%set page resolution, assume square
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 res res]);

%set font size

%old refrennce code font size code
%textobj = findobj('type', 'text');
%set(textobj, 'fontunits', 'points');
%set(textobj, 'fontsize', 20);
%set(textobj, 'fontweight', 'bold');

set(leglbl, 'fontsize', font);
set(xlbl, 'fontsize', font);
set(ylbl, 'fontsize', font);
set(titlelbl, 'fontsize', font);

%save to file (jpg)
print(gcf, '-djpeg', '../SimulationFigures/TabCusum.jpg')
