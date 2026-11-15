%RUN_FINAL_EVALUATION  Generalisation test on randomized engagements (Section VI).
%
%   The finalized design (nothing is retuned) is run on N_EVAL randomized
%   cases that were NOT used during the design: random initial geometry,
%   random two-axis maneuvers up to 98 % of a_max, and one or two radar-loss
%   intervals (onset >= 1 s, duration 0.2-0.5 s). Design alternatives are
%   run on a subset of the same cases for comparison.
%
%   Produces figures/fig_evaluation.pdf, results/evaluation_metrics.txt,
%   results/tab_ablation.tex and results/evaluation.mat.
%
%   To use the instructor's evaluation script instead, copy
%   turret_controller.m and turret_plant.m next to it: both are self-contained.

clear; close all;
here = fileparts(mfilename('fullpath'));  addpath(here);
resdir = fullfile(here, 'results');  if ~exist(resdir, 'dir'), mkdir(resdir); end
eval_cases = 101:300;          % held-out cases (cases 1-40 were used for design)
abl_cases  = 101:200;
R_eng = 2;

variants = { ...
  'Proposed design',                                   [], eval_cases; ...
  'Reference led by $t_p$ (track $\hat\theta(t+t_p)$)', struct('t_lead', 0.5), abl_cases; ...
  'Feedback only: no rate aiding, linear PI',          struct('use_ff', false, 'e_int', 1e3), abl_cases; ...
  'Near-constant-velocity model ($\sigma_m=1$)',       struct('sigma_m', 1), abl_cases; ...
  'No confidence gate ($\kappa=0$, $d_{aim}=R_{eng}$)', struct('kappa', 0, 'd_aim', 2.0), abl_cases};

for v = 1:size(variants, 1)
    turret_controller('override', variants{v, 2});
    cases = variants{v, 3};  n = numel(cases);
    E = struct('success', false(n, 1), 't_int', nan(n, 1), 'n_pulses', zeros(n, 1), 'n_hits', zeros(n, 1), ...
               'miss_rms', nan(n, 1), 'miss_max', nan(n, 1), 'point_rms', nan(n, 1), 'servo_rms', nan(n, 1), ...
               'est_rms', nan(n, 1), 'pred_rms', nan(n, 1), 'in_radius', nan(n, 1), 'V_peak', nan(n, 1), ...
               'i_peak', nan(n, 1), 'a_peak', nan(n, 1), 'n_loss', zeros(n, 1));
    E.pulse_miss = [];  E.pulse_t = [];  E.pulse_age = [];  E.z = [];  E.pulse_case = [];
    for k = 1:n
        L = simulate_engagement(cases(k));  M = L.M;
        E.success(k) = M.success;  E.t_int(k) = M.t_intercept;
        E.n_pulses(k) = M.n_pulses;  E.n_hits(k) = M.n_hits;
        E.miss_rms(k) = M.miss_rms;  E.miss_max(k) = M.miss_max;
        E.point_rms(k) = M.point_rms_mrad;  E.servo_rms(k) = M.servo_rms_mrad;  E.est_rms(k) = M.est_rms_mrad;
        E.pred_rms(k) = M.pred_cross_rms;  E.in_radius(k) = M.in_radius;
        E.V_peak(k) = M.V_peak;  E.i_peak(k) = M.i_peak;
        E.a_peak(k) = max(sqrt(L.ax.^2 + L.ay.^2));
        S = scenario_definition(cases(k));  E.n_loss(k) = size(S.blackout, 1);
        E.pulse_miss = [E.pulse_miss; L.pulse_miss];  E.pulse_t = [E.pulse_t; L.pulse_t];
        E.pulse_case = [E.pulse_case; cases(k)*ones(numel(L.pulse_t), 1)];
        % age of the newest radar sample when each pulse was released
        tv = L.t(L.meas_valid > 0);
        for j = 1:numel(L.pulse_t)
            E.pulse_age(end+1, 1) = L.pulse_t(j) - max(tv(tv <= L.pulse_t(j) + 1e-9));
        end
        % filter consistency: cross-range error normalised by its predicted sigma
        idx = find(L.t >= 2 & L.meas_valid > 0);
        E.z = [E.z; L.cross_err(idx)./L.sig_d(idx)];
    end
    EV{v} = E; %#ok<SAGROW>
    fprintf('%-52s success %3d/%3d, pulses hit %4d/%4d\n', variants{v, 1}, sum(E.success), n, sum(E.n_hits), sum(E.n_pulses));
end
turret_controller('override', []);
save(fullfile(resdir, 'evaluation.mat'), 'EV', 'variants');

%% text summary
E = EV{1};  n = numel(eval_cases);
fid = fopen(fullfile(resdir, 'evaluation_metrics.txt'), 'w');
out = @(varargin) tee(fid, varargin{:});
out('Final evaluation of the finalized design on %d held-out randomized cases\n', n);
out('cases with two radar-loss intervals: %d; peak target acceleration %.1f-%.1f m/s^2\n', sum(E.n_loss == 2), min(E.a_peak), max(E.a_peak));
out('interception success          : %d/%d (%.1f %%)\n', sum(E.success), n, 100*mean(E.success));
out('pulses hit / issued           : %d/%d (%.2f %%)\n', sum(E.n_hits), sum(E.n_pulses), 100*sum(E.n_hits)/sum(E.n_pulses));
out('pulses per engagement         : mean %.2f, min %d\n', mean(E.n_pulses), min(E.n_pulses));
out('interception time [s]         : mean %.3f, std %.3f, median %.3f, max %.3f\n', mean(E.t_int, 'omitnan'), std(E.t_int, 'omitnan'), median(E.t_int, 'omitnan'), max(E.t_int));
out('miss distance at pulse [m]    : mean %.3f, median %.3f, 95th pct %.3f, max %.3f\n', mean(E.pulse_miss), median(E.pulse_miss), p95(E.pulse_miss), max(E.pulse_miss));
out('pulses released during a radar loss (sample age > 60 ms): %d, their miss: mean %.3f, max %.3f m\n', ...
    sum(E.pulse_age > 0.06), mean(E.pulse_miss(E.pulse_age > 0.06)), max([E.pulse_miss(E.pulse_age > 0.06); NaN]));
out('beam-target miss, t >= 2 s [m]: RMS %.3f, worst case %.3f\n', sqrt(mean(E.miss_rms.^2)), max(E.miss_max));
out('pointing error RMS [mrad]     : total %.3f, estimation %.3f, servo %.3f\n', sqrt(mean(E.point_rms.^2)), sqrt(mean(E.est_rms.^2)), sqrt(mean(E.servo_rms.^2)));
out('t_p-ahead prediction error, cross-range RMS [m]: %.3f\n', sqrt(mean(E.pred_rms.^2)));
out('time inside R_eng, t >= 2 s   : %.2f %%\n', 100*mean(E.in_radius));
out('filter consistency: std of normalised cross-range error = %.2f (ideal 1), |z| > 3 in %.2f %% of samples\n', std(E.z), 100*mean(abs(E.z) > 3));
out('peak voltage %.0f V, peak current %.0f A\n', max(E.V_peak), max(E.i_peak));
out('\nDesign alternatives on cases %d-%d:\n', abl_cases(1), abl_cases(end));
ftx = fopen(fullfile(resdir, 'tab_ablation.tex'), 'w');
for v = 1:size(variants, 1)
    A = EV{v};
    m = ismember(variants{v, 3}, abl_cases);        % same subset for every variant
    A.success = A.success(m);  A.t_int = A.t_int(m);  A.n_pulses = A.n_pulses(m);
    A.n_hits = A.n_hits(m);  A.miss_rms = A.miss_rms(m);
    A.pulse_miss = A.pulse_miss(ismember(A.pulse_case, abl_cases));
    na = numel(abl_cases);
    hit = 100*sum(A.n_hits)/max(sum(A.n_pulses), 1);
    out('%-52s success %5.1f %%, pulse hit rate %5.1f %% (%d pulses), t_int %.2f s, miss at pulse %.2f m, tracking miss RMS %.2f m\n', ...
        regexprep(variants{v, 1}, '[\$\\{}]', ''), 100*sum(A.success)/na, hit, sum(A.n_pulses), mean(A.t_int, 'omitnan'), mean(A.pulse_miss), sqrt(mean(A.miss_rms.^2)));
    fprintf(ftx, '%s & %.0f & %.1f & %.2f & %.2f & %.2f \\\\\n', variants{v, 1}, 100*sum(A.success)/na, hit, ...
        mean(A.t_int, 'omitnan'), mean(A.pulse_miss), sqrt(mean(A.miss_rms.^2)));
end
fclose(fid);  fclose(ftx);

%% figure
plot_evaluation();
fprintf('Saved figures/fig_evaluation.pdf and results/evaluation_metrics.txt\n');
