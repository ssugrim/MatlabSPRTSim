% Generates a plot of Max E(n) v Delta

%needed when running on linux.
%addpath '../lib'

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%tuneing parameters - Center and width
p_prime = 0.3;
gran= 0.01;

delta = 0.01:0.01:0.25;

%Max curve for 50% split
p_0_split = 0.5;
E_n_05_max = zeros(size(delta));

for j = 1:1:length(delta)
    p_zero = p_prime - delta(j)*p_0_split;
    p_one = p_prime + delta(j)*(1-p_0_split);
    [p,L,E_n] = GenOCTheoCurve(alpha, beta,p_one,p_zero,gran);
    E_n_05_max(j)=max(E_n);
end

%Max curve for 40% split
p_0_split = 0.4;
E_n_04_max = zeros(size(delta));

for j = 1:1:length(delta)
    p_zero = p_prime - delta(j)*p_0_split;
    p_one = p_prime + delta(j)*(1-p_0_split);
    [p,L,E_n] = GenOCTheoCurve(alpha, beta,p_one,p_zero,gran);
    E_n_04_max(j)=max(E_n);
end

%Max curve for 30% split
p_0_split = 0.3;
E_n_03_max = zeros(size(delta));

for j = 1:1:length(delta)
    p_zero = p_prime - delta(j)*p_0_split;
    p_one = p_prime + delta(j)*(1-p_0_split);
    [p,L,E_n] = GenOCTheoCurve(alpha, beta,p_one,p_zero,gran);
    E_n_03_max(j)=max(E_n);
end

%Max curve for 20% split
p_0_split = 0.2;
E_n_02_max = zeros(size(delta));

for j = 1:1:length(delta)
    p_zero = p_prime - delta(j)*p_0_split;
    p_one = p_prime + delta(j)*(1-p_0_split);
    [p,L,E_n] = GenOCTheoCurve(alpha, beta,p_one,p_zero,gran);
    E_n_02_max(j)=max(E_n);
end


save('./MaxEnVdelta.mat');

