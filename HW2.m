close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A_h=1.008;
A_l=(1-0.7629*1.008)/0.2371
A=[A_h,A_l]
TPM=[0.977,1-0.977;1-0.926,0.926]
%LR_TPM=TPM^1000
%Y=rand(4,5001)


%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 10; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons_h = A(1)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 
cons_l = A(2)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 

ret_h = cons_h .^ (1 - sigma) / (1 - sigma); % return function
ret_l = cons_l .^ (1 - sigma) / (1 - sigma); % return function

% it very large negative utility
ret_h(cons_h < 0) = -Inf;
ret_l(cons_l < 0) = -Inf;

%%%% Iterat% negative consumption is not possible -> make it irrelevant by assigning

dis = [1;1]; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(2, num_k);

while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_mat_h = ret_h + beta *(TPM(1,1)* repmat(v_guess(1,:), [num_k 1])+TPM(1,2)* repmat(v_guess(2,:), [num_k 1]));
    value_mat_l = ret_l + beta *(TPM(2,1)* repmat(v_guess(1,:), [num_k 1])+TPM(2,2)* repmat(v_guess(2,:), [num_k 1]));
    
    % find the optimal k' for every k:
    [vfn_h, pol_indx_h] = max(value_mat_h, [], 2);
    vfn_h = vfn_h';
    
    [vfn_l, pol_indx_l] = max(value_mat_l, [], 2);
    vfn_l = vfn_l';
    
    % what is the distance between current guess and value function
    dis = [max(abs(vfn_h - v_guess(1,:))) ; max(abs(vfn_l - v_guess(2,:)))]
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess=[vfn_h ; vfn_l]
   
end

figure
plot(k,vfn_h,k,vfn_l)
title('Value Functions')
xlabel('k')
ylabel('V(k)')
legend('A=1.008','A=0.9743')

g_h = k(pol_indx_h); % policy function
g_l = k(pol_indx_l);

figure
plot(k,g_h,k,g_l)
title('Policy Functions')
xlabel('k')
ylabel('g(k)')
legend('A=1.008','A=0.9743')


%%Generate RAndom A sequence
rng(1234);
A_Sim=zeros(1,5001);
A_Sim(1)=A_h;
TP=rand(1,5001);
k_seq=zeros(0,5001);
k_seq(1)=25;


for i=2:5001
    [c,indx]=min(abs(k_seq(i-1)-k));
   if A_Sim(i-1)==A_h
       if TP<TPM(1,1)
           A(i)=A_h;
           k_seq(i)=g_h(indx);
       else
           A_Sim(i)=A_l;
           k_seq(i)=g_l(indx);
       end
   else
       if TP<TPM(2,2)
           A_Sim(i)=A_l;
           k_seq(i)=g_l(indx);
       else
           A_Sim(i)=A_h
           k_seq(i)=g_h(indx)
       end
   end
       
end

Y=A_Sim.*k_seq.^alpha;
std(Y)

figure
plot(Y)
title('Output over Time')
xlabel('index')
ylabel('saving')