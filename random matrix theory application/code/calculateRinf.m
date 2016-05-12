% This function is to calculate the Stiljes trasform of -\rho.


function R_inf = calculateRinf(l, rho, K, c,sigma2)
m = calculateM(l, rho, K, c);
dm = calculatedMdrho(l, m, K, c);
A = m;
B = m + rho*dm;
% temp1 = (sigma2 * (1+A)^2)/(1 + sigma2*(1+A)^2);
% temp2 = (B*sigma2*(1+A)^2*A^2)/(B + B*sigma2*(1+A)^2)^2;
% R_inf = K * max(log2(temp1 + temp2),0);


numerator = (B * sigma2 * (1 + A)^2) * (B + B * sigma2 * (1 + A)^2 + A^2);
denomenator = (B + B * sigma2 * (1 + A)^2)^2;
R_inf = K * max(log2(numerator/denomenator),0);

end