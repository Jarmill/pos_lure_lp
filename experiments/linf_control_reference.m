load('control_plant.mat');

i = 20;

%theorem 2.4 of https://ieeexplore.ieee.org/document/751839
Delta = delta_list(i)*eye(d);


cons = [];
signs = 2*ff2n(d)-1;

Ns = size(signs, 2);
AD = cell(Ns, 1);
BF = cell(Ns, 1);

for j = 1:Ns
    scurr = signs(j, :);
    AD{j} = A+B1*diag(scurr)*Delta*C1;
    BF{j} = B2 + B1*diag(scurr)*Delta*F1;
end
    
%% declare the variables
Q = sdpvar(n, n);
S = sdpvar(n, m);
gam = sdpvar(1, 1);
mu = sdpvar(1, 1);
alpha = 1; %the parameter

%form the blocks
for j = 1:Ns
    block_1 = [P, AF{j}*Q + B3*S, BF{j}]
end