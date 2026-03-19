function [A, B, C, D] = unicycle_plant(cfg)
%UNICYCLE_PLANT  Linearised state-space model of the reaction-wheel unicycle.
%
%  Physical model — two coupled inverted-pendulum axes:
%
%  PITCH axis (base wheel drives pitch recovery):
%    (I_body + m_body*L²) * θ̈_p  =  m_body*g*L*θ_p  -  τ_base
%
%  ROLL axis (reaction wheel provides roll stabilisation):
%    I_body_r * θ̈_r  =  m_body*g*L*θ_r  -  τ_rw
%    I_rw     * ω̇_rw =  τ_rw
%
%  State:   x = [θ_p, θ̇_p, θ_r, θ̇_r, ω_rw]ᵀ   (5 states)
%  Control: u = [τ_rw, τ_base]ᵀ                  (2 inputs)
%  Output:  y = [θ_p, θ_r, θ̇_p, θ̇_r]ᵀ          (4 outputs)

m  = cfg.m_body + cfg.payload_mass;
L  = cfg.L_body;
g  = cfg.g;
Ib_p = cfg.I_body   + m * L^2;   % effective pitch inertia
Ib_r = cfg.I_body_r;              % roll body inertia
Irw  = cfg.I_rw;                  % reaction wheel inertia

% A matrix (5x5)
A = zeros(5);
% pitch subsystem  [θ_p, θ̇_p]  → rows 1-2
A(1,2) =  1;
A(2,1) =  m*g*L / Ib_p;          % unstable pole

% roll subsystem   [θ_r, θ̇_r, ω_rw] → rows 3-5
A(3,4) =  1;
A(4,3) =  m*g*L / Ib_r;          % unstable pole
% ω_rw has no coupling to θ_r in A (coupling through B)

% B matrix (5x2)  columns: [τ_rw, τ_base]
B = zeros(5,2);
B(2,2) = -1 / Ib_p;              % τ_base drives pitch
B(4,1) = -1 / Ib_r;              % τ_rw drives roll
B(5,1) =  1 / Irw;               % τ_rw spins reaction wheel

% Output: observe pitch, roll, their rates, and reaction-wheel speed
C = [1 0 0 0 0;   % θ_p
     0 0 1 0 0;   % θ_r
     0 1 0 0 0;   % θ̇_p
     0 0 0 1 0;   % θ̇_r
     0 0 0 0 1];  % ω_rw encoder

D = zeros(size(C,1),2);

% Verification
fprintf('[Plant] Controllability rank: %d / %d\n', rank(ctrb(A,B)), size(A,1));
fprintf('[Plant] Observability rank:   %d / %d\n', rank(obsv(A,C)), size(A,1));
end


function xdot = unicycle_nonlinear(t, x, u, cfg)
%UNICYCLE_NONLINEAR  Full nonlinear equations of motion.
%  Used internally for validation against the linearised model.
%
%  State: x = [θ_p, θ̇_p, θ_r, θ̇_r, ω_rw]

theta_p  = x(1);
dtheta_p = x(2);
theta_r  = x(3);
dtheta_r = x(4);
omega_rw = x(5);

tau_rw   = u(1);
tau_base = u(2);

m  = cfg.m_body + cfg.payload_mass;
L  = cfg.L_body;
g  = cfg.g;
Ib_p = cfg.I_body   + m * L^2;
Ib_r = cfg.I_body_r;
Irw  = cfg.I_rw;

% Nonlinear: sin(θ) instead of θ
ddtheta_p = (m*g*L*sin(theta_p) - tau_base) / Ib_p;
ddtheta_r = (m*g*L*sin(theta_r) - tau_rw)   / Ib_r;
domega_rw = tau_rw / Irw;

xdot = [dtheta_p; ddtheta_p; dtheta_r; ddtheta_r; domega_rw];
end
