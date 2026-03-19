function unicycle_benchmark(results, cfg)
%UNICYCLE_BENCHMARK  Print and save benchmark table (mirrors benchmark.py output).
%
%  Columns match the benchmark snapshot in the project README:
%    Controller | Survival | Crash rate | Composite score
%
%  Composite score:
%    100 * survival_fraction
%   - 50 * pitch_RMS
%   - 40 * roll_RMS
%   -  5 * ctrl_RMS
%   - 20 * mean_peak_transient_angle_per_kick
%   - 15 * mean_settling_time_fraction
%
%  The deterministic benchmark disturbances are injected on the
%  reaction-wheel channel, so the transient terms are evaluated on the
%  dominant body-angle excursion within each kick window rather than on
%  pitch alone. This keeps the added score terms informative for the
%  actual benchmark loading profile.

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════╗\n');
fprintf('║  BENCHMARK RESULTS  (MATLAB simulation equivalent)  ║\n');
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  %-26s  %8s  %10s  %9s ║\n', 'Controller', 'Survival', 'Crash', 'Score');
fprintf('╠══════════════════════════════════════════════════════╣\n');

families = fieldnames(results);

% Sort by composite score descending
scores = zeros(numel(families),1);
transient_peak = zeros(numel(families),1);
transient_settle = zeros(numel(families),1);
for i = 1:numel(families)
    [transient_peak(i), transient_settle(i)] = compute_transient_metrics(results.(families{i}), cfg);
    scores(i) = compute_benchmark_score(results.(families{i}), transient_peak(i), transient_settle(i));
end
[~, idx] = sort(scores, 'descend');

lines = {};
for ii = 1:numel(idx)
    i    = idx(ii);
    fn   = families{i};
    res  = results.(fn);
    name = res.name;
    surv = res.survival_frac;
    crash_rate = 1 - surv;
    score = scores(i);
    mpk = transient_peak(i);
    mst = transient_settle(i);
    line = sprintf('║  %-26s  %8.3f  %10.3f  %9.3f ║', ...
                   name, surv, crash_rate, score);
    fprintf('%s\n', line);
    fprintf('║    [transient] peak=%6.4f rad  settling=%5.3f         ║\n', mpk, mst);
    lines{end+1} = line; %#ok<AGROW>
end
fprintf('╚══════════════════════════════════════════════════════╝\n');

%% ── Save to text file ────────────────────────────────────────
fname = fullfile(cfg.results_dir, ...
    ['benchmark_', datestr(now,'yyyymmdd_HHMMSS'), '_summary.txt']);
fid = fopen(fname, 'w');
fprintf(fid, 'MATLAB Unicycle Benchmark — %s\n', datestr(now));
fprintf(fid, 'T_sim=%.1fs  dt=%.3fs  noise=ON  disturbance=ON\n\n', cfg.T_sim, cfg.dt);
fprintf(fid, '%-28s  %8s  %10s  %9s  %9s  %9s\n', ...
    'Controller', 'Survival', 'CrashRate', 'Score', 'PeakPitch', 'SettleFrac');
fprintf(fid, '%s\n', repmat('-', 1, 62));
for i = 1:numel(families)
    fn  = families{i};
    res = results.(fn);
    [mpk, mst] = compute_transient_metrics(res, cfg);
    score = compute_benchmark_score(res, mpk, mst);
    fprintf(fid, '%-28s  %8.3f  %10.3f  %9.3f  %9.4f  %9.3f\n', ...
        res.name, res.survival_frac, 1-res.survival_frac, score, mpk, mst);
end
fclose(fid);
fprintf('[Benchmark] Saved: %s\n', fname);

%% ── Save CSV ─────────────────────────────────────────────────
csv_fname = fullfile(cfg.results_dir, ...
    ['benchmark_', datestr(now,'yyyymmdd_HHMMSS'), '.csv']);
fid = fopen(csv_fname, 'w');
fprintf(fid, 'controller,survival,crash_rate,composite,pitch_rms,roll_rms,ctrl_rms,mean_peak_pitch_per_kick,mean_settling_time_fraction\n');
for i = 1:numel(families)
    fn  = families{i};
    res = results.(fn);
    [mpk, mst] = compute_transient_metrics(res, cfg);
    score = compute_benchmark_score(res, mpk, mst);
    fprintf(fid, '"%s",%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
        res.name, res.survival_frac, 1-res.survival_frac, ...
        score, res.pitch_rms, res.roll_rms, res.ctrl_rms, mpk, mst);
end
fclose(fid);
fprintf('[Benchmark] CSV saved: %s\n', csv_fname);
end


function score = compute_benchmark_score(res, mpk, mst)
score = 100 * res.survival_frac ...
      - 50 * res.pitch_rms ...
      - 40 * res.roll_rms ...
      -  5 * res.ctrl_rms ...
      - 20 * mpk ...
      - 15 * mst;
end


function [mpk, mst] = compute_transient_metrics(res, cfg)
%COMPUTE_TRANSIENT_METRICS  Mean peak disturbance-axis angle and settling fraction.
t = res.t;
x_pitch = res.x_true(1, :);
x_roll  = res.x_true(3, :);
dt = cfg.dt;
thresh = 0.02;
win_s = 1.5;
consec = 10;

if ~isfield(cfg, 'dist_kicks') || isempty(cfg.dist_kicks)
    mpk = 0;
    mst = 0;
    return;
end

kicks = cfg.dist_kicks;
peaks = zeros(1, numel(kicks));
fracs = zeros(1, numel(kicks));

for ki = 1:numel(kicks)
    k_start = round(kicks(ki) / dt) + 1;
    k_end = min(round((kicks(ki) + win_s) / dt) + 1, numel(t));
    if k_start > numel(t)
        peaks(ki) = 0;
        fracs(ki) = 0;
        continue;
    end

    window_pitch = x_pitch(k_start:k_end);
    window_roll = x_roll(k_start:k_end);
    if max(abs(window_roll)) >= max(abs(window_pitch))
        window_angle = window_roll;
    else
        window_angle = window_pitch;
    end
    peaks(ki) = max(abs(window_angle));

    settled = false;
    if numel(window_angle) >= consec
        [~, peak_idx] = max(abs(window_angle));
        last_start = numel(window_angle) - consec + 1;
        for ks = peak_idx:last_start
            if all(abs(window_angle(ks:ks+consec-1)) < thresh)
                fracs(ki) = min(((ks - 1) * dt) / win_s, 1.0);
                settled = true;
                break;
            end
        end
    end
    if ~settled
        fracs(ki) = 1.0;
    end
end

mpk = mean(peaks);
mst = mean(fracs);
end
