%% ============================================================
%  SIM-TO-REAL SENSITIVITY ANALYSIS
%  Mirrors: final/sim2real_sensitivity.py
%
%  Sweeps key physical parameters ±X% and reports how much
%  the composite score changes — identifies fragile assumptions.
% ============================================================
function unicycle_sim2real_sensitivity()
fprintf('\n=== SIM-TO-REAL SENSITIVITY ANALYSIS ===\n\n');

cfg0 = unicycle_config();
[A0, B0, C0, ~] = unicycle_plant(cfg0);
ctrl0 = unicycle_design_controllers(A0, B0, C0, cfg0);
% Use the best single controller for sensitivity
ctrl_ref = ctrl0.baseline_mpc;

% Baseline score
res0 = unicycle_simulate(A0, B0, C0, ctrl_ref, cfg0);
score0 = res0.composite;
fprintf('Baseline composite score: %.3f\n\n', score0);

% Parameter sweep definition:  {field, ±percentage, label}
params = {
    'm_body',   15, 'Body mass ±15%';
    'L_body',   10, 'CoM height ±10%';
    'I_body',   20, 'Body pitch inertia ±20%';
    'I_body_r', 20, 'Body roll inertia ±20%';
    'I_rw',     25, 'RW inertia ±25%';
    'dt',        0, 'Control delay +1 step';  % special case
};

results_s2r = {};
for pi = 1:size(params,1)
    pname = params{pi,1};
    pct   = params{pi,2};
    plabel = params{pi,3};

    if strcmp(pname, 'dt')
        % Special: add 1 sample of extra delay
        cfg_hi = cfg0;
        cfg_hi.imu_delay_samples = cfg0.imu_delay_samples + 1;
        cfg_lo = cfg_hi;  % symmetric
        score_hi = run_with_cfg(cfg_hi, A0, B0, C0, ctrl_ref);
        score_lo = score_hi;
    else
        val0 = cfg0.(pname);
        cfg_hi = cfg0;  cfg_hi.(pname) = val0 * (1 + pct/100);
        cfg_lo = cfg0;  cfg_lo.(pname) = val0 * (1 - pct/100);

        % Recompute plant for changed params
        [Ah, Bh, ~, ~] = unicycle_plant(cfg_hi);
        [Al, Bl, ~, ~] = unicycle_plant(cfg_lo);
        score_hi = run_with_cfg(cfg_hi, Ah, Bh, C0, ctrl_ref);
        score_lo = run_with_cfg(cfg_lo, Al, Bl, C0, ctrl_ref);
    end

    delta_hi = score_hi - score0;
    delta_lo = score_lo - score0;
    sensitivity = max(abs(delta_hi), abs(delta_lo));

    results_s2r{end+1} = {plabel, score_hi, score_lo, delta_hi, delta_lo, sensitivity}; %#ok<AGROW>

    fprintf('  %-30s  +: %+7.2f  -: %+7.2f  |Δ|_max: %.2f\n', ...
        plabel, delta_hi, delta_lo, sensitivity);
end

%% ── Plot ──────────────────────────────────────────────────────
n = numel(results_s2r);
sens_vals = zeros(n,1);
labels    = {};
for i = 1:n
    sens_vals(i) = results_s2r{i}{6};
    labels{i}    = results_s2r{i}{1}; %#ok<AGROW>
end
[sv, si] = sort(sens_vals, 'descend');
figure('Name','Sim-to-Real Sensitivity','Position',[200 200 700 400]);
barh(sv, 'FaceColor', [0.9, 0.4, 0.2]);
set(gca, 'YTick', 1:n, 'YTickLabel', labels(si));
xlabel('Max |ΔScore|');
title('Sim-to-Real Sensitivity (higher = more fragile)');
grid on;
cfg_r = unicycle_config();
saveas(gcf, fullfile(cfg_r.results_dir, 'sim2real_sensitivity.png'));
fprintf('\n[S2R] Plot saved.\n');
end


function score = run_with_cfg(cfg, A, B, C, ctrl_ref)
% Run sim with modified cfg but same controller gains
try
    res = unicycle_simulate(A, B, C, ctrl_ref, cfg);
    score = res.composite;
catch
    score = -999;
end
end
