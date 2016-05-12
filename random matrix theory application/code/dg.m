% This function is to calculate the Stiljes trasform of -\rho.


function result = dg(l, rho, K, c)


temp1 = (1-c)^2/rho^2+2*(1+c)/rho +1;
temp2 = 1/4*temp1^(-1/2)*(-2*(1-c)^2/rho^3-2*(1+c)/rho^2);
result = temp2 + 1/2*c/rho^2

end