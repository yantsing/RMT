% Give the channels H, noise variance \sigma^2, and regularization
% parameter \alpha, calculate the true secrecy sum-rate using RCI precoding. 


function R = calculateR(H, sigma2, alpha)
DEBUG_A_B = 0
K = size(H,1);
N = size(H,2);
I_K = eye(K);
I_N = eye(N);
W = H' * inv(H*H' + alpha * I_K);
zeta = real(trace(W' * W))
gamm = real(trace(H'*H*(H'*H + alpha * I_N)^(-2)))

X = (H'*H + alpha * I_N)^(-1);



for k = 1 : K
    h_k = H(k,:);
%     temp = 0
%     for j = 1 : K
%         if j ~= k
%             h_j =  H(j,:);
%             temp = temp + (abs(h_k * X * h_j'))^2;
%         end
%     end
    
    
    
    
    if k == 1
        H_k = H(2:K, :);
    elseif k == K
        H_k = H(1:K-1, :);
    else
        index = [1:k-1,k+1:K]
        H_k = H(index,:);
    end
    
    temp = real(h_k* X * H_k' * H_k * X* h_k');
    
    SINR(k) = (abs(h_k * X * h_k'))^2/(zeta * sigma2 + temp);

    SINR_tilde(k) = (norm(H_k*X*h_k'))^2/(zeta * sigma2);
end

R = 0;
for i = 1 : k
    R = R + max(log2(SINR(k))-log2(SINR_tilde(k)),0);
end


for k = 1 : K
   if k == 1
        H_k(:,:,k) = H(2:K, :);
    elseif k == K
        H_k (:,:,k)= H(1:K-1, :);
    else
        index = [1:k-1,k+1:K];
        H_k(:,:,k) = H(index,:);
    end  
end

for k = 1 : K
    Z = H_k(:,:,k)' * H_k(:,:,k);
    I = size(Z, 1);
    
    
    
    A(k)= H(k,:) * inv(Z + alpha  * I) * H(k,:)';
    if DEBUG_A_B == 1
    test1 = H(k,:) * inv(H' * H + alpha * I) * H(k,:)'
    test2 = A(k)/ (1 + A(k))
    end
    
    test3 = real(H(k,:) * inv(H' * H + alpha * I) * H_k(:,:,k)' * H_k(:,:,k) * inv(H' * H + alpha * I)* H(k,:)')
    
    B(k) = H(k,:) * inv(Z + alpha  * I) * Z* inv(Z + alpha  * I) * H(k,:)';
    if DEBUG_A_B == 1
        test3 = real(H(k,:) * inv(H' * H + alpha * I) * H_k(:,:,k)' * H_k(:,:,k) * inv(H' * H + alpha * I)* H(k,:)')
        test4 =  B(k)/(1+ A(k))^2
    end
    SINR1(k) = A(k)^2/(B(k) + zeta * sigma2 * (1 + A(k))^2);
    SINR1_tilde(k) = B(k)/(zeta * sigma2 * (1 + A(k))^2);
    
    test5 = (1 + SINR1(k))/(1+ SINR1_tilde(k))
    numerator = (zeta * sigma2 * (1 + A(k))^2) * (B(k) + zeta * sigma2 * (1 + A(k))^2 + A(k)^2);
    denomenator = (B(k) + zeta * sigma2 * (1 + A(k))^2)^2;
    test6  = real(numerator/denomenator);
    Rtemp(k) = max(log2(real(numerator/denomenator)),0);
    Rtemp1(k) = max(log2(real(test5)),0)
end

R1 = 0;
for i = 1 : k
    R1 = R1 + Rtemp(k) ;
end

R = R1

% R = 0;
% for i = 1 : k
%     R = R + max(log2(1 + SINR(k))-log2(1 + SINR_tilde(k)),0);
% end



end