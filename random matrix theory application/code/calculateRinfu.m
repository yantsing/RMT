% This function is to calculate the Stiljes trasform of -\rho.


function R_inf = calculateRinfu(l, rho, K, c,sigma2)
m = calculateg(l, rho, K, c);
dm = calculatedg(l, rho, K, c);
A = m;
B = m + rho*dm;
temp1 = (sigma2 * (1+A)^2)/(1 + sigma2*(1+A)^2);
temp2 = (B*sigma2*(1+A)^2*A^2)/(B + B*sigma2*(1+A)^2)^2;
R_inf = K * max(log2(temp1 + temp2),0);

end