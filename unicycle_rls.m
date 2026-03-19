function varargout = unicycle_rls(action, varargin)
%UNICYCLE_RLS  Single-file RLS module for online pitch-axis identification.
switch lower(action)
    case 'init'
        varargout{1} = init_rls(varargin{:});
    case 'update'
        varargout{1} = update_rls(varargin{:});
    case 'plant'
        [varargout{1}, varargout{2}] = plant_rls(varargin{:});
    otherwise
        error('unicycle_rls:UnknownAction', ...
            'Unknown action "%s". Use init, update, or plant.', action);
end
end


function rls = init_rls(A_nom, B_nom, cfg)
%INIT_RLS  Initialise the scalar-forgetting-factor RLS estimator.
theta0 = [A_nom(2,1); -B_nom(2,2)];

rls.theta = theta0;
rls.theta_init = theta0;
rls.P = 100 * eye(2);
rls.y_filt = 0;
rls.lambda = cfg.rls_lambda;
if isfield(cfg, 'rls_P_max')
    rls.P_max = cfg.rls_P_max;
else
    rls.P_max = 50;
end
rls.update_count = 0;
rls.gain_update_count = 0;
end


function rls = update_rls(rls, phi, y, cfg)
%UPDATE_RLS  One-step RLS update for pitch stiffness/control authority.
lambda = cfg.rls_lambda;
phi = phi(:);
if isfield(cfg, 'rls_accel_alpha')
    alpha = cfg.rls_accel_alpha;
else
    alpha = 0.3;
end
rls.y_filt = alpha * y + (1 - alpha) * rls.y_filt;
y_use = rls.y_filt;

den = lambda + phi' * rls.P * phi;
if den <= eps
    return;
end

K_rls = (rls.P * phi) / den;
innovation = y_use - phi' * rls.theta;
rls.theta = rls.theta + K_rls * innovation;
rls.P = (rls.P - K_rls * phi' * rls.P) / lambda;
scale = min(1, rls.P_max ./ max(diag(rls.P), 1e-12));
if any(scale < 1)
    rls.P = diag(scale) * rls.P * diag(scale);
end
rls.update_count = rls.update_count + 1;

% Keep the identified parameters physically meaningful and numerically tame.
rls.theta(1) = min(max(rls.theta(1), 1e-6), 100);
rls.theta(2) = min(max(rls.theta(2), 1e-6), 100);
end


function [A_est, B_est] = plant_rls(rls, A_nom, B_nom)
%PLANT_RLS  Inject the identified pitch coefficients into the nominal model.
A_est = A_nom;
B_est = B_nom;

A_est(2,1) = rls.theta(1);
B_est(2,2) = -rls.theta(2);
end
