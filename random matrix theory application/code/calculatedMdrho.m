% This function is to calculate the Stiljes trasform of -\rho.


function dm = calculatedMdrho(l,m, K, c)

% The following should be wrong!
% temp_sum = 0;
% for k = 1 : K
%     temp_sum = temp_sum + l(k)^2/(1 + l(k) * m)^2;
% end
% dm = 1/(-m^2 + (c/K)*temp_sum);


%%%%%%%%%%%%%%%%I am not sure whether it is correct %%%%%%%%%%%%%%%%%%%%%%
dm = 0;

termiation_condition = 0.0001;

previousdM = dm + 6*termiation_condition;
count = 0;
while abs(previousdM - dm) > termiation_condition
    temp_sum = 0;
    for k = 1 : K
        temp_sum = temp_sum + l(k)^2*dm/((1 + l(k) * m)^2);
    end
    previousdM = dm;
    dm = -(1 - (c/K)*temp_sum)*m^2;
    count = count + 1;
end

% delta = 0.00001;
% m = calculateM(l, rho, K, c);
% mPlus = calculateM(l, rho+delta/2, K, c);
% mMinus = calculateM(l, rho-delta/2, K, c);
% dm = (mPlus - mMinus)/delta;

end