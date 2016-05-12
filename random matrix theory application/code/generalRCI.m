% This program simulate the CoMP conference paper.

clear all
DEBUG = 0
UTEST = 1

t0 = clock;
%Parameter: # of MS, # of BS, # of receiving antenna, # of transmitting
%antenna.
N = 4; 
K = 4;

sigma2 = 0.01;

%%%%%%%%%%%%%%%%%%%Generate h(:,k)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1: K
%     l(k) = rand()
%     generatePathLoss
    l(k) = 1;
end
% l = generatePathLoss(K);
if DEBUG == 1
    l(1) = 1;
    l(2) = 0.1;
    l(3) = 0.1;
    l(4) = 0.1;
end

%Generate the channel h_k
for k = 1: K
    x = (randn(1,N)+sqrt(-1)*randn(1,N))/sqrt(2);
    H(k,:) = l(k) * x;
end

Xi = diag(rand(1,K));

[U,D,V] = svd(H);
Lambda = D*D;
e = eye(K,K);
for k = 1 : k
    I(:,:,k) = e(:,k)*e(:,k)';
end

W = H'*inv(H*H'+ U*Xi*U');
%%%%%%%%%%Calculating vector x %%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    x(k) = 1/(Lambda(k,k)+Xi(k,k));
end
x = x';

%%%%%%%%%%Calculating vector b %%%%%%%%%%%%%%%%%%%%%%
b = diag(Lambda);

%%%%%%%%%%Calculating B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = b*b';

%%%%%%%%%%Calculating vectors a_kj %%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    for j = 1 : K
        a(:,k,j) = (H(k,:) * V * D * diag(U(j,:)'))';
    end
end

%%%%%%%%%%%%%%%Claculating t_k%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    if UTEST == 1
        % Teh following quantities should be the same.
        t_test(k)= H(k,:)*H'*inv(H*H'+ U*Xi*U')*e(:,k)
        t_test1(k) = H(k,:)* W(:,k)
        t(k) = H(k,:) * H'* U * inv(Lambda + Xi)*U(k,:)'
        t(k) = H(k,:) * H'* U * diag(U(k,:)') * x
        
        t(k) = H(k,:) * V * D * diag(U(k,:)')* x
        t(k) = a(:,k,k)'*x
    else
        t(k) = a(:,k,k)'*x;
    end
end

%%%%%%%%%%%%%%%%Claculating t_{\tilde{k}} %%%%%%%%%%%%%%%%%%%
for k = 1 : K
    temp_test = 0;
    temp = 0;
    temp_test1 = 0;
    temp_test2 = 0;
    for j = 1 : K
        if j ~= k
            if UTEST == 1
                temp_test = temp_test + abs(H(k,:)* W(:,j))^2;
                temp_test1 = temp_test1 + abs(H(k,:) * V * D * diag(U(j,:)')* x)^2;
                temp = temp + abs(a(:,k,j)'*x)^2;
            else
                temp = temp + abs(a(:,k,j)'*x)^2;
            end
        end
    end
    tTilde2_test(k) = temp_test
    tTilde2_test(k) = temp_test1
    tTilde2(k)=temp
end

%%%%%%%%%%%%%%%%%Calculating A_k,A_{\tilde{k}}%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    A(:,:,k) = a(:,k,k)* a(:,k,k)';
    temp = zeros(K);
    for j = 1 : K
        if j ~= k
                temp = temp + a(:,k,j)*a(:,k,j)';
        end
    end
    ATilde(:,:,k) = temp
end

%%%%%%%%%%%%%%%Calculating SINR_k%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    if UTEST == 1 
        same11 =  norm(a(:,k,k)'*x)^2
        same12  = x'*A(:,:,k)*x
        
        same21 = tTilde2(k) + sigma2 * trace(Lambda * x * x')
        same22 =  x' * (ATilde(:,:,k) + sigma2 * Lambda) * x
        
    else
        numerator = real(x'*A(:,:,k)*x);
        denomenator = real(x' * (ATilde(:,:,k) + sigma2 * Lambda) * x);
    end
end


% t = 2;
% cvx_begin
% variable Q(K,K) symmetric
% maximize( trace(Q) )
% for k = 1 : K
%     trace(t * (ATilde(:,:,k) + sigma2 * Lambda)*Q) - trace(A(:,:,k) * Q)  <= 0
% end
% subject to
% for k = 1 : K
%     trace(I(:,:,k) * Q) <= Lambda(k,k)^(-2)
% end
% cvx_end


t = 2;
cvx_begin
variable x(K)
maximize(sum(x))
for k = 1 : K
    norm(t*( ATilde(:,:,k) + sigma2 * Lambda)^(1/2) *x) <=  real(a(:,k,k))'* x
end
subject to
for k = 1 : K
  x(k) <= Lambda(k,k)^(-1)
  x(k) >= 0.001
end
cvx_end







    




