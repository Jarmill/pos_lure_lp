function [gain, out] = pos_linf_lure_elem(sys, Delta, ep)
%POS_LINF_LURE_ELEM Positive Linfinity norm estimation under elementwise-bounded lure
%uncertainty

if nargin < 3
    ep = 1e-3;
end

[n, e] = size(sys.B1);

AF = sys.A + sys.B1 * Delta * sys.C1;

%declare variables
vinf = sdpvar(n, 1);
gaminf = sdpvar(1, 1);

cons_inf = [vinf >= ep; AF*vinf - vinf + sum(sys.B2, 2) <= 0; sys.C2*vinf - gaminf + sum(sys.F2, 2) <= 0];


%form the problem 
opts = sdpsettings('solver', 'mosek');

%solve the problem
sol = optimize(cons_inf, gaminf, opts);

%recover the solution

status = sol.problem;
if status
    out = [];
    gain = inf;
else
    v_rec = value(vinf);
    gam_rec = value(gaminf);
    gain = gam_rec;
    out = struct('v', v_rec, 'gam', gam_rec);
end

end

