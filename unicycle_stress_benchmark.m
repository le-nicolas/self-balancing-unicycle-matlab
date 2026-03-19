%% ============================================================
%  UNICYCLE STRESS BENCHMARK
%  Mirrors: final/results/throughput_ablation_seed42.json
%           (seed=42, 80 episodes, high disturbance profile)
%
%  Usage:  unicycle_stress_benchmark
%          unicycle_stress_benchmark(20)   % override episode count
% ============================================================
function unicycle_stress_benchmark(n_episodes)
if nargin < 1, n_episodes = 20; end  % default 20 (use 80 to match original)

rng(42);  % seed = 42, matches project convention
fprintf('\n=== STRESS BENCHMARK  (seed=42, %d episodes) ===\n\n', n_episodes);

cfg = unicycle_config();

% Override to the seed-42 stress profile used in the original project.
cfg.dist_amp      = 0.30;
cfg.dist_duration = 0.05;
% Randomise kick timing across episodes
kick_times = 1.0 + 6.0 * rand(n_episodes, 3);  % up to 3 kicks per episode

[A, B, C, ~] = unicycle_plant(cfg);
ctrl_all = unicycle_design_controllers(A, B, C, cfg);
families = fieldnames(ctrl_all);

results_stress = struct();
for fi = 1:numel(families)
    fn = families{fi};
    survived = 0;
    pitch_rms_all = zeros(n_episodes,1);

    for ep = 1:n_episodes
        % Randomise initial condition
        x0 = [0.02*(2*rand()-1);
               0.01*(2*rand()-1);
               0.02*(2*rand()-1);
               0.01*(2*rand()-1);
               0];
        cfg_ep = cfg;
        cfg_ep.x0_nominal   = x0;
        cfg_ep.dist_t_kick  = kick_times(ep,1);
        cfg_ep.dist_amp     = cfg.dist_amp * (0.7 + 0.6*rand());
        cfg_ep.dist_sine_amp = 0.0;

        evalc('res = unicycle_simulate(A, B, C, ctrl_all.(fn), cfg_ep);');
        if res.survival, survived = survived + 1; end
        pitch_rms_all(ep) = res.pitch_rms;
    end

    results_stress.(fn).survival_rate = survived / n_episodes;
    results_stress.(fn).crash_rate    = 1 - survived / n_episodes;
    results_stress.(fn).mean_pitch_rms = mean(pitch_rms_all);
    results_stress.(fn).name = ctrl_all.(fn).name;
end

%% ── Print table ──────────────────────────────────────────────
fprintf('%-28s  %8s  %10s  %12s\n', ...
    'Controller','Survival','Crash Rate','Mean Pitch RMS');
fprintf('%s\n', repmat('-',1,65));
for fi = 1:numel(families)
    fn = families{fi};
    s  = results_stress.(fn);
    fprintf('%-28s  %8.3f  %10.3f  %12.6f\n', ...
        s.name, s.survival_rate, s.crash_rate, s.mean_pitch_rms);
end

%% ── Save JSON-style text ─────────────────────────────────────
cfg2 = unicycle_config();
fname = fullfile(cfg2.results_dir, ...
    sprintf('throughput_ablation_seed42_%depisodes.txt', n_episodes));
fid = fopen(fname, 'w');
fprintf(fid, '{\n  "seed": 42, "episodes": %d, "results": {\n', n_episodes);
for fi = 1:numel(families)
    fn = families{fi};
    s  = results_stress.(fn);
    fprintf(fid, '    "%s": {"survival": %.4f, "crash_rate": %.4f, "pitch_rms": %.6f}', ...
        fn, s.survival_rate, s.crash_rate, s.mean_pitch_rms);
    if fi < numel(families), fprintf(fid, ','); end
    fprintf(fid, '\n');
end
fprintf(fid, '  }\n}\n');
fclose(fid);
fprintf('\n[Stress] Results saved: %s\n', fname);
end
