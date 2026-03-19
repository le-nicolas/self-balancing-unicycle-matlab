function [Ad, Bd] = unicycle_discretize(A, B, dt)
%UNICYCLE_DISCRETIZE  Exact zero-order-hold discretization.
%  Uses the Van Loan block exponential, which remains valid when A is singular.

nx = size(A, 1);
nu = size(B, 2);

M = expm([A, B; zeros(nu, nx + nu)] * dt);
Ad = M(1:nx, 1:nx);
Bd = M(1:nx, nx+1:nx+nu);
end
