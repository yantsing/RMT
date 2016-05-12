% This function is to calculate the Stiljes trasform of -\rho.


function m = calculateM(l, rho, K, c)
m = 1;

termiation_condition = 0.0001;

previousM = m + 6*termiation_condition;
count = 0;
while abs(previousM - m) > termiation_condition
    temp_sum = 0;
    for k = 1 : K
        temp_sum = temp_sum + l(k)/(1 + l(k) * m);
    end
    previousM = m;
    m = (rho + (c/K)*temp_sum)^(-1);
    count = count + 1;
end

end