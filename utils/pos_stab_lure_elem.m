function [status, out] = pos_stab_lure_elem(sys, Delta, ep)
%POS_STAB_LURE_ELEM positive stabilization under elementwise-bounded lure
%uncertainty

if nargin < 3
    ep = 1e-3;
end

[n, e] = size(sys.B1);

AF = sys.A + sys.B1 * Delta * sys.C1;

v = sdpvar(n, 1);
dv = AF'*v - v;

%form the problem 
cons_stab = [v >= ep; dv <= -ep*v; sum(v) == 1];
opts = sdpsettings('solver', 'mosek');

%solve the problem
sol = optimize(cons_stab, max(v), opts);

%recover the solution
status = sol.problem;
if status
    out = [];
else
    v_rec = value(v);
    out = struct('v', v_rec);
end

end

