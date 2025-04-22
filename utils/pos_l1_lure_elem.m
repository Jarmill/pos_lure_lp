function [gain, out] = pos_l1_lure_elem(sys, Delta, ep)
%POS_L1_LURE_ELEM Positive L1 norm estimation under elementwise-bounded lure
%uncertainty

if nargin < 3
    ep = 1e-3;
end

[n, e] = size(sys.B1);

AF = sys.A + sys.B1 * Delta * sys.C1;

%declare variables
v1 = sdpvar(n, 1);
gam1 = sdpvar(1, 1);

cons_1 = [v1 >= ep; AF'*v1 - v1 + sum(sys.C2', 2) <= 0; sys.B2'*v1 - gam1 + sum(sys.F2', 2) <= 0];


%form the problem 
opts = sdpsettings('solver', 'mosek');

%solve the problem
sol = optimize(cons_1, gam1, opts);

%recover the solution

status = sol.problem;
if status
    out = [];
    gain = inf;
else
    v_rec = value(v1);
    gam_rec = value(gam1);
    gain = gam_rec;
    out = struct('v', v_rec, 'gam', gam_rec);
end

end

