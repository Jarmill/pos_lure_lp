%Plotting a Lure trajectory

rng(30,'twister')
s = [0, 0.1, 0.6, 0.2, 0.1];
beta = [0.95; 0.9; 0.7; 0.5];

A = [s; diag(beta), zeros(4, 1)];

n = 5;
m = 2;

B11 = [eye(m); zeros(n-m, m)];
C11 = [0, 0, 1, 1, 0; 0, 0, 0, 1, 1];
delta = 0.0995;
Delta = delta*eye(m);
D = zeros(2, 2);

%% certify stability
v = sdpvar(n, 1);

sys = ss(A, B11, C11, D, 1);

ep = 1e-3;
dv = (A+B11*Delta*C11)'*v - v;
% cons = [v >= ep; dv <= -ep*v; sum(v) == 1];
cons_stab = [v >= ep; dv <= -ep*v; sum(v) == 1];
opts = sdpsettings('solver', 'mosek');
optimize(cons_stab, max(v), opts);
% optimize(cons_inf, gaminf, opts);
% 
v_rec = value(v);


%% simulate the system
% ustar = 0*ones(m, 1);
% ustar = [0.1; 0];
% ustar = [0; 0];
ustar = [0.2; 0.1];

T = 1000;
Nsim = 10;
trange = 0:T;

% x0c = cell(Nsim, 1);
Level = 3;
x0c = [ones(n, 1), Level*rand(n, Nsim-1)];
Xc = cell(Nsim, 1);
Uc = cell(Nsim, 1);
Xd = zeros(n, Nsim);

% u = @(zcurr) delta/2*(zcurr - sin(zcurr)) + ustar;
u = @(zcurr) delta*sin(zcurr) + ustar;
phases = (2*pi/n)*(0:(n-1))';
w = @(t) cos(t*2*pi/100 + phases);

for i = 1:Nsim
    x0 = x0c(:, i);
    X = zeros(n, T+1);
    U = zeros(m, T);
    X(:, 1) = x0;
    for t=1:T
        xcurr = X(:, t);
        zcurr = C11*xcurr;       
        ucurr = u(zcurr);
        xnext = A*xcurr + B11*ucurr + w(t);
        X(:, t+1) = xnext;
        U(:, t) = ucurr;
    end
    Xc{i} = X;
    Uc{i} = U;
    Xd(:, i) = xnext;
end

disp(Xd)
figure(3)
clf
subplot(1, 2, 1)
plot(trange, Xc{1})
subplot(1, 2, 2)
plot(trange, Xc{2})