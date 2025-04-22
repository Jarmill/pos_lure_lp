%Example 5.2

rng(30,'twister')
% s = [0, 0.1, 0.6, 0.2, 0.1];
% beta = [0.95; 0.9; 0.7; 0.5];

% A = [s; diag(beta), zeros(4, 1)];



n = 10;
d = 2;
e = 2;
m=4;
p=2;

A = rand(n, n)*0.25;

B1 = [eye(d); zeros(n-d, d)];
B2 = [eye(2); zeros(n-2, 2)];
B3 = [zeros(n-6, m); eye(m); zeros(2, m)];

% C1 = [0, 0, 1, 1, 0, 0; 0, 0, 0, 1, 1, 0];
C1 = [eye(d), zeros(d, n-d)];
C2 = [zeros(2, n-2), eye(2)];

% D1 = zeros(d, m);
D1 = 0.1*[eye(2), zeros(2, 2)];
D2 = zeros(p, m);

% F1 = [0.1];
F1 = 0.2*eye(e);
% F1 = 0;
F2 = zeros(1, e);

% C12 = [1, zeros(1, n-1)];

%[delta, l1 gain, linf gain]
delta_list = 0:0.005:0.4;
% delta = 0.1
% delta = 0.075
% delta = 0.0995;
% delta = 0.05;
% delta = 0.025;
% delta = 0; %1.2347    2.4693
% delta_list = 0.05;
% delta_list = 0.4;
% delta_list = 0;
% delta_list = 0.31373;
delta_crit  = 0.31373;

ld = length(delta_list);
l1_gain =zeros(ld, 1);
linf_gain =zeros(ld, 1);
l1_status =zeros(ld, 1);
linf_status =zeros(ld, 1);
linf_K =cell(ld, 1);
%% find linf gain
dindex = 1;


for i = 1:ld
% for i =dindex:dindex
    vinf = sdpvar(n, 1);
    gaminf = sdpvar(1, 1);
    S = sdpvar(m, n);
    Delta = delta_list(i)*eye(d);
    ep = 1e-3;
    AD = A+B1*Delta*C1;
    BF = B2 + B1*Delta*F1;
    dv = (AD)*vinf + sum(B1*Delta*D1*S, 2) + sum(B3*S, 2) - vinf;

    V = diag(vinf);
    cons_nonneg = [vec(A*V+B3*S) >= 0;
        % vec(AD*V+B3*S + B1*Delta*D1*S) >= 0;
        vec(C1*V + D1*S) >= 0; 
        vec(C2*V + D2*S) >= 0];
    cons_inf = [vinf >= ep; dv + sum(BF, 2) <= 0; C2*vinf - gaminf + sum(F2, 2) + sum(D2*S, 2)<= 0; gaminf >= 0];
            % vec(C1*V + D1*S) >= 0; 
    % cons = [v >= ep; dv <= -ep*v; sum(v) == 1];
    % cons_stab = [v >= ep; dv <= -ep*v; sum(v) == 1];
    opts = sdpsettings('solver', 'mosek');
    solinf = optimize([cons_inf, cons_nonneg], gaminf, opts);
    % optimize(cons_inf, gaminf, opts);
    % 
    vinf_rec = value(vinf);
    gaminf_rec = value(gaminf);
    S_rec = value(S);
    K_rec = S_rec * diag(1./vinf_rec);
    linf_gain(i) = gaminf_rec;
    linf_status(i) = solinf.problem;
    linf_K{i} = K_rec;
    
   
end


%% simulate the system
% ustar = 0*ones(m, 1);
ustar = [0.5; 0];


T = 600;
Nsim = 10;
trange = 0:T;

% x0c = cell(Nsim, 1);
x0c = 2*[ones(n, 2), rand(n, Nsim-2)];
Xc = cell(Nsim, 1);
Uc = cell(Nsim, 1);
Yc = cell(Nsim, 1);
Xd = zeros(n, Nsim);
Yd = zeros(p, Nsim);

z = @(zcurr) delta_list(dindex)/2*(zcurr + sin(1*zcurr)) + ustar;
% u = @(zcurr) delta_list(dindex)*zcurr + ustar;

wlevel = 0.05;
Kcurr = linf_K{dindex};
% wlevel = 0;
% u = @(zcurr) -delta*zcurr + ustar;
% w = @(t) wlevel*(2*rand(e, 1) - 1);
w = @(t) wlevel*(randi(2, e, 1) - 1.5)*2;
for i = 1:Nsim
    x0 = x0c(:, i);
    X = zeros(n, T+1);
    Z = zeros(d, T);
    U = zeros(m, T);
    W = zeros(e, T);
    Y = zeros(p, T);
    X(:, 1) = x0;    
    for t=1:T
        xcurr = X(:, t);
        if i<=3
            wcurr = zeros(e, 1);
        else
            wcurr = w(t);
        end
        ucurr = Kcurr*xcurr;
        % ucurr = zeros;
        zetacurr = C1*xcurr + F1*wcurr + D1*ucurr;       
        zcurr = z(zetacurr);        
        xnext = A*xcurr + B1*zcurr + B2*wcurr + B3*ucurr;
        xnext = max(xnext, zeros(n, 1));
        X(:, t+1) = xnext;
        Z(:, t) = zcurr;
        W(:, t) = wcurr;
        U(:, t) = ucurr;
        ycurr = C2*xcurr + F2*wcurr + D2*ucurr;
        Y(:, t) = ycurr;
    end
    Xc{i} = X;
    Uc{i} = Z;
    Yc{i} = Y;
    Uc{i} = U;
    Xd(:, i) = xnext;
    Yd(:, i) = ycurr;
end

% disp([delta, gam1_rec, gaminf_rec])
disp(gaminf_rec)
disp(Xd)
figure(3)
clf
subplot(1, 2, 1)
plot(trange, Xc{1})
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$x_t$', 'interpreter', 'latex', 'fontsize', 14)
title('$w_t = 0$', 'interpreter', 'latex', 'fontsize', 16)
subplot(1, 2, 2)
plot(trange, Xc{2})
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$x_t$', 'interpreter', 'latex', 'fontsize', 14)
title('$|w_t| < 0.1$', 'interpreter', 'latex', 'fontsize', 16)

ystar = Yd(end, 1);
figure(4)
clf
hold on
ybound = linf_gain(dindex)*wlevel;
for i =1:Nsim
    plot(trange(1:end-1), Yc{i});    
end

plot(xlim, [1, 1]*(ystar+ybound), 'k')
plot(xlim, [1, 1]*(ystar-ybound), 'k')
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$y_t$', 'interpreter', 'latex', 'fontsize', 14)
title('$\ell_\infty$ bounds on output', 'interpreter', 'latex', 'fontsize', 16)


%% gain feasiblity
figure(50)
clf
hold on
feas = linf_status==0;
semilogy(delta_list(feas), linf_gain(feas), 'linewidth', 3)
yl = ylim;
yl(1) = 1;
ylim(yl);
plot([1, 1]*delta_crit, yl, '--k')
xlabel('$\tau$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$\ell_\infty$ gain bound', 'interpreter', 'latex', 'fontsize', 14)
title('Bounds on controlled $\ell_\infty$ norm', 'interpreter', 'latex', 'fontsize', 16)
set(gca, 'YScale', 'log')
title('')
