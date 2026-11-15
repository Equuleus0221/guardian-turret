function plot_evaluation()
%PLOT_EVALUATION  Fig. 5 of the paper from results/evaluation.mat.
%   Called by run_final_evaluation; can be re-run alone to restyle the figure.
here = fileparts(mfilename('fullpath'));
D = load(fullfile(here, 'results', 'evaluation.mat'));
E = D.EV{1};  n = numel(E.success);  R_eng = 2;
fs = 7;  cB = [0 0.30 0.65];  cR = [0.85 0.33 0.10];

fig = figure('Visible', 'off');
ax = axes('Position', [0.075 0.19 0.25 0.72]); hold on;
edges = 0:0.05:2.2;
cnt = histc(E.pulse_miss, edges); %#ok<HISTC>
bar(edges + 0.025, cnt, 1, 'FaceColor', cB, 'EdgeColor', 'none');
yl = [0, 1.1*max(cnt)];
plot([R_eng R_eng], yl, 'k--', 'LineWidth', 0.6);
text(R_eng - 0.05, 0.9*yl(2), 'R_{eng}', 'FontSize', fs-1, 'HorizontalAlignment', 'right');
xlim([0 2.2]); ylim(yl); box on; grid on; set(ax, 'FontSize', fs);
xlabel('Miss distance at pulse (m)'); ylabel('Number of pulses');
title(sprintf('(a) %d pulses, %d cases', numel(E.pulse_miss), n), 'FontSize', fs, 'FontWeight', 'normal');

ax = axes('Position', [0.405 0.19 0.25 0.72]); hold on;
ts = sort(E.t_int(~isnan(E.t_int)));
stairs([0; ts; 10], [0; (1:numel(ts))'/n; numel(ts)/n]*100, '-', 'Color', cB, 'LineWidth', 1.2);
xlim([0 2]); ylim([0 105]); box on; grid on; set(ax, 'FontSize', fs);
xlabel('Interception time (s)'); ylabel('Cases intercepted (%)');
title('(b) Time of first successful pulse', 'FontSize', fs, 'FontWeight', 'normal');

ax = axes('Position', [0.735 0.19 0.25 0.72]); hold on;
plot(E.a_peak(E.n_loss == 1), E.miss_max(E.n_loss == 1), 'o', 'Color', cB, 'MarkerSize', 2.5);
plot(E.a_peak(E.n_loss == 2), E.miss_max(E.n_loss == 2), 's', 'Color', cR, 'MarkerSize', 2.5);
plot([11 20.5], [R_eng R_eng], 'k--', 'LineWidth', 0.6);
xlim([11.5 20]); ylim([0, max(5.5, 1.35*max(E.miss_max))]); box on; grid on; set(ax, 'FontSize', fs);
xlabel('Peak target acceleration (m/s^2)'); ylabel('Worst beam-target miss (m)');
legend({'one radar loss', 'two radar losses'}, 'FontSize', fs-1.5, 'Location', 'northwest', 'Box', 'off');
title('(c) Worst tracking miss, t \geq 2 s', 'FontSize', fs, 'FontWeight', 'normal');
save_fig(fig, 'fig_evaluation', 7.16, 2.0);
close(fig);
end
