function [L2err,H1err,err,pf,tf] = Compute_L2error(a, c, f, g, p, e, t, u, u_ex)
% Compute L2- and H1-error from between u (grid function, i.e., vector of length size(p,2)) and
% the exact solution u_ex, using a quadrature rule

% refine mesh to obtain accurate quadrature rule
[pf, ef, tf, Uf] = refinemesh(g, p, e, t, u);
% assemble matrices for norm computations
[Kf, Mf, ~] = assema(pf, tf, a, c, 0);
err = Uf - u_ex(pf).'; % linear interpolant of error on refined grid
L2err = sqrt(err' * Mf * err);  % L2-norm of (interpolant of) error
H1err = sqrt(err' * (Mf+Kf) * err); % % H1-norm of (interpolant of) error
end
