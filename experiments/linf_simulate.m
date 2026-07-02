rng(30,'twister')

BOUND = 1;
SIMULATE =0;

s = [0, 0.1, 0.6, 0.2, 0.1];
beta = [0.95; 0.9; 0.7; 0.5];

A = [s; diag(beta), zeros(4, 1)];

n = 5;
m = 2;

e = 2;
p=1;
B1 = [eye(m); zeros(n-m, m)];
B2 = [eye(e); zeros(n-e, e)];
C1 = [0, 0, 1, 1, 0; 0, 0, 0, 1, 1];
C2 = [zeros(p, n-p), eye(p)];
F1 = 0.1*ones(m, e);
F2 = zeros(p, e);


% C12 = [1, zeros(1, n-1)];

%[delta, l1 gain, linf gain]
% delta_list = 0:0.025:0.25;
% delta_list = 0:0.01:0.25;

%main list
% delta_list = 0:0.0025:0.105;


delta_list = 0:0.025:0.1;

% delta = 0.1
% delta = 0.075
% delta = 0.0995;
% delta = 0.05;
% delta = 0.025;
% delta = 0; %1.2347    2.4693



%% find linf gain
if BOUND
    ld = length(delta_list);
l1_gain =zeros(ld, 1);
l1_d_gain =zeros(ld, 1);
% l2_gain =zeros(ld, 1);
linf_gain =zeros(ld, 1);
l1_status =zeros(ld, 1);
l2_status =zeros(ld, 1);
linf_status =zeros(ld, 1);
for i = 1:ld
% dindex = 3;
% for i =dindex:dindex
    vinf = sdpvar(n, 1);
    gaminf = sdpvar(1, 1);
    Delta = delta_list(i)*eye(m);
    ep = 1e-3;
    AD = A+B1*Delta*C1;
    BF = B2 + B1*Delta*F1;
    dv = (AD)*vinf - vinf;
    cons_inf = [vinf >= ep; dv + sum(BF, 2) <= 0; C2*vinf - gaminf + sum(F2, 2) <= 0];
    % cons = [v >= ep; dv <= -ep*v; sum(v) == 1];
    % cons_stab = [v >= ep; dv <= -ep*v; sum(v) == 1];
    opts = sdpsettings('solver', 'mosek', 'verbose', 0);
    solinf = optimize(cons_inf, gaminf, opts);
    % optimize(cons_inf, gaminf, opts);
    % 
    vinf_rec = value(vinf);
    gaminf_rec = value(gaminf);
    linf_gain(i) = gaminf_rec;
    linf_status(i) = solinf.problem;
    
    %% find l1 gain
    v1 = sdpvar(n, 1);
    gam1 = sdpvar(1, 1);
    dv1 = AD'*v1 - v1;
    
    cons_1 = [v1 >= ep; dv1 + sum(C2', 2) <= 0; BF'*v1 - gam1 + sum(F2', 2) <= 0];
    % cons = [v >= ep; dv <= -ep*v; sum(v) == 1];
    % cons_stab = [v >= ep; dv <= -ep*v; sum(v) == 1];
    opts = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol1 = optimize(cons_1, gam1, opts);
    % optimize(cons_inf, gaminf, opts);
    % 
    v1_rec = value(v1);
    gam1_rec = value(gam1);
    l1_gain(i) = gam1_rec;
    l1_status(i) = sol1.problem;

    %% find l1 gain through the alternative formulation
    v1_d = sdpvar(n, 1);
    gam1_d = sdpvar(1, 1);
    dv1_d = AD'*v1_d - v1_d;
    
    cons_1_d = [v1_d >= ep];
    
    for j = 1:size(BF, 2)
        dterm = dv1_d + BF(:, j);
        wterm = sum(C2, 1)*v1_d - gam1_d + sum(F2(:, j));
        cons_1_d = [cons_1_d;  dterm <= 0; wterm <= 0];
    end


    % cons_1_d = [v1_d >= ep; dv + sum(BF, 2) <= 0; C2*vinf - gaminf + sum(F2, 2) <= 0];
 % cons = [v >= ep; dv <= -ep*v; sum(v) == 1];
    % cons_stab = [v >= ep; dv <= -ep*v; sum(v) == 1];
    opts = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol1_d = optimize(cons_1_d, gam1_d, opts);
    % optimize(cons_inf, gaminf, opts);
    % 
    v1_d_rec = value(v1_d);
    gam1_d_rec = value(gam1_d);
    l1_d_gain(i) = gam1_d_rec;
    l1_d_status(i) = sol1_d.problem;
    
    
    % [l2_gain(i), l2_status(i)] = l2_pos_gain(AD, BF, C2, F2);
   

    [l1_gain, l1_d_gain]
end
end

%% simulate the system
% ustar = 0*ones(m, 1);
if SIMULATE
ustar = [0.5; 0];


T = 600;
Nsim = 20;
trange = 0:T;

% x0c = cell(Nsim, 1);
x0c = 2*[ones(n, 2), rand(n, Nsim-2)];
Xc = cell(Nsim, 1);
Uc = cell(Nsim, 1);
Yc = cell(Nsim, 1);
Xd = zeros(n, Nsim);
Yd = zeros(p, Nsim);

dindex = 21;

u = @(zcurr) delta_list(dindex)/2*(zcurr + sin(1*zcurr)) + ustar;
% u = @(zcurr) delta_list(dindex)*zcurr + ustar;

wlevel = 0.05;
% wlevel = 0;
% wlevel = 0;
% u = @(zcurr) -delta*zcurr + ustar;
% w = @(t) wlevel*(2*rand(e, 1) - 1);
% w = @(t) wlevel*(randi(2, e, 1) - 1.5)*2;

% u = @(zcurr) delta*sin(zcurr) + ustar;
% phases = (2*pi/T)*(0:(T-1))';
% phases = zeros(trange)

% w0 = @(t) 0.15*cos(t*2*pi/100 + 0.0001*t.^2);
w0 = @(t) 0*t;
w = @(t) w0(t) + wlevel*(randi(2, e, 1) - 1.5)*2;

% gain_bound = linf_gain(dindex)*2*wlevel;
for i = 1:Nsim
    x0 = x0c(:, i);
    X = zeros(n, T+1);
    U = zeros(m, T);
    W = zeros(e, T);
    Y = zeros(p, T);
    X(:, 1) = x0;    
    for t=1:T
        xcurr = X(:, t);
        if i< 3
            wcurr = zeros(e, 1);
        else
            wcurr = w(t);
        end
        zcurr = C1*xcurr + F1*wcurr;       
        ucurr = u(zcurr);        
        xnext = A*xcurr + B1*ucurr + B2*wcurr;
        xnext = max(xnext, zeros(n, 1));
        X(:, t+1) = xnext;
        U(:, t) = ucurr;
        W(:, t) = wcurr;
        ycurr = C2*xcurr + F2*wcurr;
        Y(:, t) = ycurr;
    end
    Xc{i} = X;
    Uc{i} = U;
    Yc{i} = Y;
    Xd(:, i) = xnext;
    Yd(:, i) = ycurr;
end

% disp([delta, gam1_rec, gaminf_rec])

%% plots

% disp(Xd)
figure(3)

clf
hold on
for i = 10:15
% plot(trange, Xc{5})
    xi = Xc{i}
    plot(trange(1:300), xi(3, 1:300))
end
xfix = Xc{1};
W0list = w0((1:300)');
% plot(trange(1:300), gain_bound + W0list + xfix(3, end), 'k');
% plot(trange(1:300), -gain_bound + W0list + xfix(3, end), 'k');
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$x^3_t$', 'interpreter', 'latex', 'fontsize', 14)
title('$|w_t^1 - w_t^2| < 0.1$', 'interpreter', 'latex', 'fontsize', 16)

% clf
% subplot(2, 1, 1)
% plot(trange, Xc{1})
% xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
% ylabel('$x_t$', 'interpreter', 'latex', 'fontsize', 14)
% title('$w_t = 0$', 'interpreter', 'latex', 'fontsize', 16)
% subplot(2, 1, 2)
% plot(trange, Xc{5})
% xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
% ylabel('$x_t$', 'interpreter', 'latex', 'fontsize', 14)
% title('$|w_t| < 0.1$', 'interpreter', 'latex', 'fontsize', 16)



% 
% disp(Xd)
% figure(3)
% clf
% subplot(2, 1, 1)
% plot(trange, Xc{1})
% xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
% ylabel('$x_t$', 'interpreter', 'latex', 'fontsize', 14)
% title('$w_t = 0$', 'interpreter', 'latex', 'fontsize', 16)
% subplot(2, 1, 2)
% plot(trange, Xc{5})
% xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
% ylabel('$x_t$', 'interpreter', 'latex', 'fontsize', 14)
% title('$|w_t| < 0.1$', 'interpreter', 'latex', 'fontsize', 16)
% 
ystar = Yd(end, 1);
figure(4)
clf
hold on
ybound = linf_gain(dindex)*wlevel*2;
for i =1:Nsim
    plot(trange(1:end-1), Yc{i});    
end
xlim([0, 300]);
ylim([0, 1.5])
plot(xlim, [1, 1]*(ystar+ybound), 'k')
plot(xlim, [1, 1]*(ystar-ybound), 'k')
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$y_t$', 'interpreter', 'latex', 'fontsize', 14)
% title('$\ell_\infty$ bounds on output', 'interpreter', 'latex', 'fontsize', 16)
end
% 
figure(3)
clf
hold on
plot(delta_list(~l1_status), l1_gain(~l1_status), 'LineWidth', 3)
plot(delta_list(~linf_status), linf_gain(~linf_status), '--', 'linewidth', 3)
set(gca, 'YScale', 'log')
xlabel('$\tau$', 'Interpreter','latex','FontSize',16)
ylabel('Incremental Gain', 'Interpreter','latex', 'FontSize',16)
legend({'$\ell_1$', '$\ell_\infty$'}, 'Interpreter','latex', 'FontSize',14, 'location', 'northwest')
% end