% This function is to calculate the Stiljes trasform of -\rho.


function result = g(l, rho, K, c)

result = 1/2*(sqrt((1-c)^2/rho^2+2*(1+c)/rho +1) +(1-c)/rho -1);

end