%RUN_SCENARIOS  Standardized scenarios A, B and C (Section V of the paper).
%
%   Produces figures/fig_scenarios.pdf, figures/fig_geometry.pdf,
%   results/scenario_metrics.txt and results/tab_scenarios.tex.
%   One representative run per scenario is plotted; the table reports the
%   mean over N_MC measurement-noise realisations.

clear; close all;
here = fileparts(mfilename('fullpath'));  addpath(here);
resdir = fullfile(here, 'results');  if ~exist(resdir, 'dir'), mkdir(resdir); end
turret_controller('override', []);
N_MC = 50;
names = {'A', 'B', 'C'};
R_eng = 2;

%% representative runs
for s = 1:3
    Lr{s} = simulate_engagement(names{s}); %#ok<SAGROW>
end

%% Monte Carlo over measurement noise
fields = {'t_intercept', 'n_pulses', 'n_hits', 'pulse_miss_mean', 'pulse_miss_max', 'miss_rms', 'miss_max', ...
          'point_rms_mrad', 'servo_rms_mrad', 'est_rms_mrad', 'cross_rms', 'pred_cross_rms', 'in_radius', 'V_peak', 'i_peak', 'V_rms_track'};
for s = 1:3
    for f = 1:numel(fields), R.(fields{f}) = nan(N_MC, 1); end
    succ = 0;
    for k = 1:N_MC
        L = simulate_engagement(struct('id', names{s}, 'seed', 100*s + k));
        for f = 1:numel(fields), R.(fields{f})(k) = L.M.(fields{f}); end
        succ = succ + L.M.success;
    end
    R.success = succ;
    MC{s} = R; %#ok<SAGROW>
    fprintf('scenario %s: %d/%d runs intercepted\n', names{s}, succ, N_MC);
end

%% text and LaTeX tables
fid = fopen(fullfile(resdir, 'scenario_metrics.txt'), 'w');
ftx = fopen(fullfile(resdir, 'tab_scenarios.tex'), 'w');
rows = { ...
  'Interception success (runs)',          @(R) sprintf('%d/%d', R.success, N_MC); ...
  'Pulses hit / issued',                  @(R) sprintf('%d/%d', sum(R.n_hits), sum(R.n_pulses)); ...
  'Interception time, mean (max) [s]',    @(R) sprintf('%.2f (%.2f)', mean(R.t_intercept, 'omitnan'), max(R.t_intercept)); ...
  'Miss at pulse, mean (max) [m]',        @(R) sprintf('%.2f (%.2f)', mean(R.pulse_miss_mean), max(R.pulse_miss_max)); ...
  'Beam--target miss, RMS [m]',           @(R) sprintf('%.2f', sqrt(mean(R.miss_rms.^2))); ...
  'Beam--target miss, worst [m]',         @(R) sprintf('%.2f', max(R.miss_max)); ...
  'Pointing error, RMS [mrad]',           @(R) sprintf('%.3f', sqrt(mean(R.point_rms_mrad.^2))); ...
  '\quad estimation part [mrad]',         @(R) sprintf('%.3f', sqrt(mean(R.est_rms_mrad.^2))); ...
  '\quad servo part [mrad]',              @(R) sprintf('%.3f', sqrt(mean(R.servo_rms_mrad.^2))); ...
  '$t_p$-ahead prediction error, RMS [m]', @(R) sprintf('%.2f', sqrt(mean(R.pred_cross_rms.^2))); ...
  'Time inside $R_{eng}$ [\%]',           @(R) sprintf('%.1f', 100*mean(R.in_radius)); ...
  'Tracking voltage, RMS [V]',            @(R) sprintf('%.0f', sqrt(mean(R.V_rms_track.^2)))};
fprintf(fid, 'Mean over %d noise realisations per scenario; tracking statistics for t >= 2 s\n\n', N_MC);
fprintf(fid, '%-42s %14s %14s %14s\n', 'metric', 'A', 'B', 'C');
for r = 1:size(rows, 1)
    v = cellfun(@(R) rows{r, 2}(R), MC, 'UniformOutput', false);
    plain = regexprep(rows{r, 1}, '\\quad |\$|\\', '');
    fprintf(fid, '%-42s %14s %14s %14s\n', plain, v{:});
    fprintf(ftx, '%s & %s & %s & %s \\\\\n', rows{r, 1}, v{:});
end
fprintf(fid, '\nRepresentative runs (default seed):\n');
for s = 1:3
    M = Lr{s}.M;
    fprintf(fid, '%s: t_int = %.3f s, pulses at %s s, miss %s m, V_peak = %.0f V, i_peak = %.0f A\n', names{s}, ...
        M.t_intercept, mat2str(Lr{s}.pulse_t', 4), mat2str(round(Lr{s}.pulse_miss', 2)), M.V_peak, M.i_peak);
end
fclose(fid);  fclose(ftx);
type(fullfile(resdir, 'scenario_metrics.txt'));

%% Figure: time histories
fs = 7;
cT = [0 0 0];  cR = [0.85 0.33 0.10];  cB = [0 0.30 0.65];  cG = [0.45 0.45 0.45];
fig = figure('Visible', 'off');
x0 = [0.065 0.385 0.705];  wdt = 0.285;
y0 = [0.775 0.545 0.315 0.085];  hgt = 0.185;
ttl = {'Scenario A: nominal', 'Scenario B: maneuver', 'Scenario C: maneuver + radar loss'};
for s = 1:3
    L = Lr{s};  t = L.t;
    S = scenario_definition(names{s});
    % (1) elevation
    ax = axes('Position', [x0(s) y0(1) wdt hgt]); hold on; shade(ax, S, [0 40]);
    plot(t, L.theta*180/pi, '-', 'Color', cT, 'LineWidth', 1.4);
    plot(t, L.theta_p*180/pi, ':', 'Color', cG, 'LineWidth', 0.9);
    plot(t, L.theta_t*180/pi, '--', 'Color', cB, 'LineWidth', 1.0);
    xlim([0 10]); ylim([0 40]); box on; grid on; set(ax, 'FontSize', fs, 'XTickLabel', []);
    title(ttl{s}, 'FontSize', fs, 'FontWeight', 'bold');
    if s == 1
        ylabel('Elevation (deg)');
        legend({'target \theta', 'predicted \theta(t+t_p)', 'turret \theta_t'}, 'FontSize', fs-1.5, 'Location', 'northwest', 'Box', 'off');
    end
    % (2) pointing error split
    ax = axes('Position', [x0(s) y0(2) wdt hgt]); hold on; shade(ax, S, [-1.5 1.5]);
    plot(t, 1e3*L.e_est, '-', 'Color', cR, 'LineWidth', 0.7);
    plot(t, -1e3*L.e_ref, '-', 'Color', cB, 'LineWidth', 0.7);
    xlim([0 10]); ylim([-1.5 1.5]); box on; grid on; set(ax, 'FontSize', fs, 'XTickLabel', []);
    if s == 1
        ylabel('Error (mrad)');
        legend({'estimation \theta_{ref}-\theta', 'servo \theta_t-\theta_{ref}'}, 'FontSize', fs-1.5, 'Location', 'northeast', 'Box', 'off');
    end
    % (3) beam-target miss distance and pulses
    ax = axes('Position', [x0(s) y0(3) wdt hgt]); hold on; shade(ax, S, [-3 3]);
    plot([0 10], [R_eng R_eng], 'k--', 'LineWidth', 0.5);  plot([0 10], -[R_eng R_eng], 'k--', 'LineWidth', 0.5);
    plot(t, L.miss, '-', 'Color', cB, 'LineWidth', 0.7);
    ip = round(L.pulse_t/(t(2) - t(1))) + 1;
    plot(L.pulse_t, L.miss(ip), 'v', 'Color', cR, 'MarkerFaceColor', cR, 'MarkerSize', 3.5);
    xlim([0 10]); ylim([-3 3]); box on; grid on; set(ax, 'FontSize', fs, 'XTickLabel', []);
    if s == 1
        ylabel('Miss (m)');
        text(9.9, 2.45, '\pm R_{eng}', 'FontSize', fs-1.5, 'HorizontalAlignment', 'right');
    end
    if s == 3, text(9.9, -2.4, 'triangles: laser pulses', 'FontSize', fs-1.5, 'HorizontalAlignment', 'right', 'Color', cR); end
    % (4) cross-range estimation error and 3-sigma bound
    ax = axes('Position', [x0(s) y0(4) wdt hgt]); hold on; shade(ax, S, [-4 4]);
    plot(t, 3*L.sig_d, '-', 'Color', cG, 'LineWidth', 0.7);  plot(t, -3*L.sig_d, '-', 'Color', cG, 'LineWidth', 0.7);
    plot(t, L.cross_err, '-', 'Color', cR, 'LineWidth', 0.7);
    xlim([0 10]); ylim([-4 4]); box on; grid on; set(ax, 'FontSize', fs);
    xlabel('Time (s)');
    if s == 1
        ylabel('Cross-range est. error (m)');
        text(9.9, 3.2, 'grey: \pm3\sigma_d from the filter', 'FontSize', fs-1.5, 'HorizontalAlignment', 'right');
    end
end
save_fig(fig, 'fig_scenarios', 7.16, 4.6);
close(fig);

%% Figure: engagement geometry
fig = figure('Visible', 'off');
ax = axes('Position', [0.13 0.17 0.84 0.80]); hold on;
phi = linspace(-10, 85, 200)*pi/180;
plot(800*cos(phi), 800*sin(phi), ':', 'Color', cG, 'LineWidth', 0.7);
plot(2000*cos(phi), 2000*sin(phi), ':', 'Color', cG, 'LineWidth', 0.7);
LA = Lr{1};  LB = Lr{2};
plot(LA.x, LA.y, '-', 'Color', cB, 'LineWidth', 1.1);
plot(LB.x, LB.y, '-', 'Color', cR, 'LineWidth', 1.1);
for n = [101 102 103 104 105 106]
    [xe, ye] = target_truth(LA.t, n);
    plot(xe, ye, '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.5);
end
plot(LA.x, LA.y, '-', 'Color', cB, 'LineWidth', 1.1);
plot(LB.x, LB.y, '-', 'Color', cR, 'LineWidth', 1.1);
plot(LA.x(1), LA.y(1), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
plot(0, 0, 'ks', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
text(30, 70, 'turret', 'FontSize', fs-1);
text(560, 150, 'r = 800 m', 'FontSize', fs-1, 'Color', cG);
text(1540, 1010, 'r = 2000 m', 'FontSize', fs-1, 'Color', cG, 'Rotation', -52);
text(LA.x(1) - 20, LA.y(1) + 45, 't = 0', 'FontSize', fs-1, 'HorizontalAlignment', 'right');
text(1000, 470, 'A', 'FontSize', fs, 'Color', cB, 'FontWeight', 'bold');
text(1010, 690, 'B, C', 'FontSize', fs, 'Color', cR, 'FontWeight', 'bold');
axis equal; xlim([0 2050]); ylim([-50 1250]); box on; grid on; set(ax, 'FontSize', fs);
xlabel('x (m)'); ylabel('y (m)');
save_fig(fig, 'fig_geometry', 3.45, 2.2);
close(fig);
fprintf('Saved figures/fig_scenarios.pdf, figures/fig_geometry.pdf and results/scenario_metrics.txt\n');
