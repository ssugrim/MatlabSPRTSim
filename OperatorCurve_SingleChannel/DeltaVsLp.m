%Generates the data for a plot of \delta vs \beta

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%tuneing parameters
p_prime = 0.9;
delta = 0.001:0.001:0.09;
tail_thres = 0.0001;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%dummy parameter
h = -500:.0001:500;

%L has no p_i dependance, only p
L = ((A.^h) - ones(size(h))) ./ (A.^h - B.^h);

%data holder
beta_data = zeros(size(delta));

for i = 1:1:length(delta)
  %compute p
  p_one = p_prime + delta(i);
  p_zero = p_prime - delta(i);
  p_upper = (1 - p_one) / (1 - p_zero);
  p_lower = (p_one) / (p_zero);
  p = (ones(size(h)) - p_upper.^h) ./ (p_lower.^h - p_upper.^h);
  
  %Isolate the delcare H_1 set
  A = p( p > p_prime);
  B = L( p > p_prime);
  
  %Get rid of trailing zeros
  A_san = A(B > tail_thres);
  B_san = B(B > tail_thres);
  
  %interpolate for regular spacing
  p_interp = p_prime:0.001:1;
  L_interp = interp1(A_san,B_san,p_interp);
  
  %dump any discontinities (resultsing from h values 
  %that are close to p_prime)
  L_interp(isnan(L_interp)) = 0 ;
  
  %store
  beta_data(i)=trapz(p_interp,L_interp);
end


