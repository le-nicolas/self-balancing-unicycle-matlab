function unicycle_plot_all(results, cfg)
%UNICYCLE_PLOT_ALL  Generate all simulation plots.
%
%  Figures:
%   1. Per-controller time histories (pitch, roll, control)
%   2. Phase portrait (pitch vs dpitch)
%   3. Composite benchmark bar chart
%   4. Paper-style comparison (pitch envelope)

families = fieldnames(results);
colors   = lines(numel(families));
res_dir  = cfg.results_dir;

%% ── Figure 1: Time histories (one subplot group per controller) ──
fig1 = figure('Name','Time Histories','Position',[50 50 1200 800]);
n = numel(families);
cols = 3;
rows = ceil(n/cols);

for i = 1:n
    fn  = families{i};
    res = results.(fn);
    t   = res.t;
    N1  = res.crash_step;   % valid data length

    subplot(rows*3, cols, [i, i+cols, i+2*cols]);   % tall column
    % Try a proper 3-row layout instead:
end
close(fig1);

% Proper layout: one figure per controller
for i = 1:n
    fn  = families{i};
    res = results.(fn);
    t   = res.t;
    Nv  = min(res.crash_step, size(res.x_true,2));

    fig = figure('Name', res.name, 'Position', [100+i*20, 80, 780, 580]);

    subplot(3,1,1);
    plot(t(1:Nv), rad2deg(res.x_true(1,1:Nv)), 'b-', 'LineWidth', 1.5); hold on;
    plot(t(1:Nv), rad2deg(res.x_est(1,1:Nv)),  'b--', 'LineWidth', 1);
    yline(rad2deg(cfg.pitch_crash_rad),  'r--', 'Crash limit');
    yline(-rad2deg(cfg.pitch_crash_rad), 'r--');
    if res.crashed
        xline(res.crash_t, 'r:', 'CRASH');
    end
    xlabel('Time [s]'); ylabel('Pitch [deg]');
    title([res.name, ' — Pitch angle']);
    legend('True','Estimated','Crash limit','Location','best');
    grid on;

    subplot(3,1,2);
    plot(t(1:Nv), rad2deg(res.x_true(3,1:Nv)), 'm-', 'LineWidth', 1.5); hold on;
    plot(t(1:Nv), rad2deg(res.x_est(3,1:Nv)),  'm--', 'LineWidth', 1);
    yline(rad2deg(cfg.roll_crash_rad),  'r--', 'Crash limit');
    yline(-rad2deg(cfg.roll_crash_rad), 'r--');
    if res.crashed
        xline(res.crash_t, 'r:', 'CRASH');
    end
    xlabel('Time [s]'); ylabel('Roll [deg]');
    title('Roll angle');
    legend('True','Estimated','Location','best');
    grid on;

    subplot(3,1,3);
    Nu = min(res.crash_step-1, size(res.u,2));
    plot(t(1:Nu), res.u(1,1:Nu), 'g-',  'LineWidth', 1.5); hold on;
    plot(t(1:Nu), res.u(2,1:Nu), 'k-',  'LineWidth', 1.5);
    yline( cfg.tau_rw_max,   'g--', 'RW limit');
    yline(-cfg.tau_rw_max,   'g--');
    yline( cfg.tau_base_max, 'k--', 'Base limit');
    yline(-cfg.tau_base_max, 'k--');
    xlabel('Time [s]'); ylabel('Torque [N·m]');
    title('Control inputs');
    legend('τ_{rw}','τ_{base}','RW limit','','Base limit','Location','best');
    grid on;

    saveas(fig, fullfile(res_dir, ['timehist_', fn, '.png']));
end

%% ── Figure 2: All controllers — pitch comparison ──────────────
fig2 = figure('Name','Pitch Comparison','Position',[200 200 900 500]);
hold on;
for i = 1:n
    fn  = families{i};
    res = results.(fn);
    t   = res.t;
    Nv  = min(res.crash_step, size(res.x_true,2));
    plot(t(1:Nv), rad2deg(res.x_true(1,1:Nv)), ...
         'Color', colors(i,:), 'LineWidth', 1.8, 'DisplayName', res.name);
end
yline(rad2deg(cfg.pitch_crash_rad),  'r--', 'LineWidth', 1.5, 'DisplayName','Crash limit');
yline(-rad2deg(cfg.pitch_crash_rad), 'r--', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('Time [s]'); ylabel('Pitch angle [deg]');
title('Controller Comparison — Pitch (all families)');
legend('Location','best'); grid on;
saveas(fig2, fullfile(res_dir, 'comparison_pitch.png'));

%% ── Figure 3: Benchmark bar chart ─────────────────────────────
fig3 = figure('Name','Benchmark','Position',[300 100 700 400]);
scores = zeros(n,1);
names  = {};
for i = 1:n
    fn = families{i};
    scores(i) = results.(fn).composite;
    names{i}  = results.(fn).name; %#ok<AGROW>
end
[scores_s, idx_s] = sort(scores, 'descend');
bar(scores_s, 'FaceColor', [0.2, 0.6, 0.9]);
set(gca, 'XTick', 1:n, 'XTickLabel', names(idx_s), ...
    'XTickLabelRotation', 20);
ylabel('Composite Score');
title('Controller Benchmark — Composite Score (higher = better)');
grid on; ylim([0, max(scores_s)*1.2]);
saveas(fig3, fullfile(res_dir, 'benchmark_bar.png'));

%% ── Figure 4: Phase portrait — pitch ─────────────────────────
fig4 = figure('Name','Phase Portrait','Position',[400 200 600 500]);
hold on;
for i = 1:n
    fn  = families{i};
    res = results.(fn);
    Nv  = min(res.crash_step, size(res.x_true,2));
    plot(rad2deg(res.x_true(1,1:Nv)), rad2deg(res.x_true(2,1:Nv)), ...
         'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', res.name);
    plot(rad2deg(res.x_true(1,1)),    rad2deg(res.x_true(2,1)), ...
         'o', 'Color', colors(i,:), 'MarkerSize', 8, 'HandleVisibility','off');
end
xlabel('Pitch angle [deg]'); ylabel('Pitch rate [deg/s]');
title('Phase Portrait — Pitch Axis'); legend('Location','best'); grid on;
saveas(fig4, fullfile(res_dir, 'phase_portrait.png'));

%% ── Figure 5: Reaction wheel speed ────────────────────────────
fig5 = figure('Name','Wheel Speed','Position',[450 250 900 400]);
hold on;
for i = 1:n
    fn  = families{i};
    res = results.(fn);
    t   = res.t;
    Nv  = min(res.crash_step, size(res.x_true,2));
    plot(t(1:Nv), res.x_true(5,1:Nv), ...
         'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', res.name);
end
yline( cfg.omega_rw_max,  'r--', 'LineWidth',1.2, 'DisplayName','Budget limit');
yline(-cfg.omega_rw_max,  'r--', 'LineWidth',1.2, 'HandleVisibility','off');
yline( cfg.omega_rw_hard, 'k--', 'LineWidth',1.2, 'DisplayName','Hard cutoff');
yline(-cfg.omega_rw_hard, 'k--', 'LineWidth',1.2, 'HandleVisibility','off');
xlabel('Time [s]'); ylabel('\omega_{rw} [rad/s]');
title('Reaction Wheel Speed'); legend('Location','best'); grid on;
saveas(fig5, fullfile(res_dir, 'wheel_speed.png'));

fprintf('[Plot] Figures saved to %s\n', res_dir);
end
