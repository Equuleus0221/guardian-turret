%DESIGN_CONTROLLER  Plant derivation check, compensator design and loop analysis.
%
%   Reproduces the numbers of Sections III-IV of the paper and Fig. 3
%   (Bode diagram, root locus and step response). Toolbox-free.

clear; close all;
here = fileparts(mfilename('fullpath'));  addpath(here);
resdir = fullfile(here, 'results');  if ~exist(resdir, 'dir'), mkdir(resdir); end
turret_controller('override', []);
P = turret_plant('params');
C = turret_controller('params');
fid = fopen(fullfile(resdir, 'design_summary.txt'), 'w');
out = @(varargin) tee(fid, varargin{:});

%% Plant
out('================ PLANT ================\n');
out('J_eq = %.0f kg m^2 (turret %.0f + reflected rotor %.0f)\n', P.Jeq, P.J_turret, P.N^2*P.J_m);
out('Gp(s) = %.0f / (%.0f s^3 + %.0f s^2 + %.0f s)\n', P.num, P.den(1), P.den(2), P.den(3));
out('      = %.4f / ( s (s^2 + %.3f s + %.1f) )\n', P.num/P.den(1), P.den(2)/P.den(1), P.den(3)/P.den(1));
pp = roots(P.den);
out('poles: 0, %.3f +/- j%.3f   (wn = %.2f rad/s, zeta = %.4f)\n', real(pp(2)), abs(imag(pp(2))), P.wn_p, P.zeta_p);
out('velocity gain k_v = %.4e rad/s per V  (%.0f V per rad/s)\n', P.k_v, 1/P.k_v);
out('resonant magnification 1/(2 zeta) = %.2f (%.1f dB)\n', 1/(2*P.zeta_p), 20*log10(1/(2*P.zeta_p)));

% numerical check of turret_plant.m against the transfer function:
% steady-state speed for 100 V and resonance frequency of the free response
th = 0; om = 0; cur = 0; h = 1e-4; n = 30000; omlog = zeros(n, 1);
for k = 1:n
    [a1, b1, c1] = turret_plant(0, 100, th, om, cur);
    [a2, b2, c2] = turret_plant(0, 100, th + h/2*a1, om + h/2*b1, cur + h/2*c1);
    [a3, b3, c3] = turret_plant(0, 100, th + h/2*a2, om + h/2*b2, cur + h/2*c2);
    [a4, b4, c4] = turret_plant(0, 100, th + h*a3, om + h*b3, cur + h*c3);
    th = th + h/6*(a1 + 2*a2 + 2*a3 + a4);  om = om + h/6*(b1 + 2*b2 + 2*b3 + b4);
    cur = cur + h/6*(c1 + 2*c2 + 2*c3 + c4);  omlog(k) = om;
end
osc = omlog - 100*P.k_v;  zc = find(osc(1:end-1).*osc(2:end) < 0);
wd_sim = pi/mean(diff(zc(1:10))*h);
out('check of turret_plant.m: speed at 100 V = %.5f rad/s (TF: %.5f), ringing at %.2f rad/s (TF: %.2f)\n', ...
    omlog(end), 100*P.k_v, wd_sim, abs(imag(pp(2))));

%% Requirements that drive the design
R_eng = 2; r_min = 800; r_max = 2000; v_max = 150; a_max = 20;
out('\n================ REQUIREMENTS ================\n');
out('pointing tolerance R_eng/r: %.2f mrad at %d m ... %.2f mrad at %d m\n', 1e3*R_eng/r_min, r_min, 1e3*R_eng/r_max, r_max);
out('max LOS rate v_max/r_min = %.4f rad/s  -> needs %.0f V of back-emf\n', v_max/r_min, v_max/r_min/P.k_v);
al_max = a_max/r_min + v_max^2/r_min^2;
out('max LOS acceleration a_max/r_min + v_max^2/r_min^2 = %.4f rad/s^2\n', al_max);
out('type-1 loop would need Kv >= %.0f 1/s for 1 mrad at max LOS rate -> type-2 loop chosen\n', v_max/r_min/1e-3);

%% Compensator and loop
A = loop_analysis(C, P);
out('\n================ COMPENSATOR ================\n');
out('Gc(s) = Kp (s + z_i)/s * (s^2 + 2 zeta_p wn s + wn^2)/(s^2 + 2 zeta_d wn s + wn^2)\n');
out('Kp = %.0f V/rad, z_i = %.2f rad/s (Ki = %.0f V/(rad s)), wn = %.2f rad/s, zeta_p = %.4f, zeta_d = %.2f\n', ...
    C.Kp, C.z_i, C.Ki, C.wn_p, C.zeta_p, C.zeta_d);
out('closed-loop poles:\n');
cp = sort(A.cl_poles, 'descend', 'ComparisonMethod', 'real');
for k = 1:numel(cp)
    if imag(cp(k)) >= 0
        if abs(imag(cp(k))) < 1e-6, out('   %9.3f\n', real(cp(k)));
        else, out('   %9.3f +/- j%.3f  (wn = %.2f, zeta = %.3f)\n', real(cp(k)), imag(cp(k)), abs(cp(k)), -real(cp(k))/abs(cp(k))); end
    end
end
out('gain crossover %.2f rad/s, phase margin %.1f deg, gain margin %.1f dB at %.1f rad/s\n', A.wc, A.PM, A.GM_dB, A.w_pc);
out('delay margin %.1f ms, Ms = %.2f, Mt = %.2f, -3 dB bandwidth %.1f rad/s, noise bandwidth %.1f rad/s\n', ...
    1e3*A.delay_margin, A.Ms, A.Mt, A.w_bw, A.w_neb);
out('step: overshoot %.1f %%, rise time %.3f s, 2 %% settling %.2f s\n', A.overshoot, A.t_rise, A.t_settle2);
out('error constants: Kp = inf, Kv = inf, Ka = %.1f 1/s^2\n', A.Ka);
out('steady-state error for LOS acceleration: %.3f mrad at %.3f rad/s^2 (worst case), %.3f mrad at 0.02 rad/s^2\n', ...
    1e3*al_max/A.Ka, al_max, 1e3*0.02/A.Ka);

% why the notch is needed: same PI without the notch
Cn = C;  Cn.Gc_num = C.Kp*[1, C.z_i];  Cn.Gc_den = [1, 0];
An = loop_analysis(Cn, P);
out('PI alone (no notch), same gains: stable = %d, max Re(pole) = %+.2f, |L| at resonance = %.2f\n', ...
    An.stable, max(real(An.cl_poles)), max(abs(An.L(An.w > 40 & An.w < 70))));
% largest gain crossover a PI could reach without the notch (6 dB gain margin at resonance)
out('PI alone: 6 dB margin at the resonance needs Kp <= %.0f V/rad, i.e. crossover <= %.1f rad/s\n', ...
    0.5*P.wn_p*2*P.zeta_p/P.k_v, 0.5*P.wn_p*2*P.zeta_p);

% robustness to a mismatch between notch and plant resonance (inductance error)
out('\nrobustness to resonance mismatch (notch fixed):\n');
for sc = [0.8 0.9 1.0 1.1 1.2]
    Ar = loop_analysis(C, P, sc);
    out('   wn x %.1f: stable = %d, PM = %.1f deg, GM = %.1f dB\n', sc, Ar.stable, Ar.PM, Ar.GM_dB);
end
fclose(fid);

%% Figure: Bode, root locus, step response
fs = 7;
fig = figure('Visible', 'off');
% --- Bode magnitude
ax1 = axes('Position', [0.07 0.60 0.36 0.35]);
Gp = polyval(P.num, 1i*A.w)./polyval(P.den, 1i*A.w);
semilogx(A.w, 20*log10(abs(C.Kp*Gp)), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8); hold on;
semilogx(An.w, An.mag_dB, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8);
semilogx(A.w, A.mag_dB, '-', 'Color', [0 0.30 0.65], 'LineWidth', 1.2);
plot([0.1 1000], [0 0], 'k:', 'LineWidth', 0.5);
plot(A.wc, 0, 'o', 'Color', [0 0.30 0.65], 'MarkerSize', 3.5, 'MarkerFaceColor', [0 0.30 0.65]);
ylim([-80 80]); xlim([0.3 1000]); grid on; set(gca, 'FontSize', fs, 'XTickLabel', []);
ylabel('Magnitude (dB)');
legend({'K_pG_p', 'PI\cdotG_p (no notch)', 'L = G_cG_p'}, 'FontSize', fs-1, 'Location', 'southwest', 'Box', 'off');
title('(a) Open-loop frequency response', 'FontSize', fs, 'FontWeight', 'normal');
% --- Bode phase
ax2 = axes('Position', [0.07 0.14 0.36 0.40]);
semilogx(An.w, An.phase_deg, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8); hold on;
semilogx(A.w, A.phase_deg, '-', 'Color', [0 0.30 0.65], 'LineWidth', 1.2);
plot([0.1 1000], [-180 -180], 'k:', 'LineWidth', 0.5);
plot([A.wc A.wc], [-180, -180 + A.PM], '-', 'Color', [0 0.30 0.65], 'LineWidth', 0.8);
text(A.wc*1.15, -180 + A.PM/2, sprintf('PM = %.0f%s', A.PM, char(176)), 'FontSize', fs-1);
text(A.w_pc*1.25, -150, sprintf('GM = %.1f dB', A.GM_dB), 'FontSize', fs-1);
ylim([-450 -90]); xlim([0.3 1000]); set(gca, 'YTick', -450:90:-90, 'FontSize', fs); grid on;
xlabel('Frequency (rad/s)'); ylabel('Phase (deg)');
% --- root locus in Kp
ax3 = axes('Position', [0.52 0.14 0.20 0.81]);
kk = logspace(-2, 1.2, 1500);
num1 = conv(C.Gc_num/C.Kp, A.Gp_num);
rl = zeros(numel(A.L_den) - 1, numel(kk));
for j = 1:numel(kk)
    q = A.L_den;  q(end-numel(num1)+1:end) = q(end-numel(num1)+1:end) + kk(j)*C.Kp*num1;
    rl(:, j) = roots(q);
end
plot(real(rl(:)), imag(rl(:)), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 2); hold on;
ol = roots(A.L_den);  oz = roots(num1);
plot(real(ol), imag(ol), 'kx', 'MarkerSize', 5, 'LineWidth', 0.8);
plot(real(oz), imag(oz), 'ko', 'MarkerSize', 4, 'LineWidth', 0.8);
plot(real(A.cl_poles), imag(A.cl_poles), 's', 'Color', [0 0.30 0.65], 'MarkerFaceColor', [0 0.30 0.65], 'MarkerSize', 3.5);
plot([0 0], [-90 90], 'k:', 'LineWidth', 0.5);
xlim([-60 15]); ylim([-90 90]); grid on; set(gca, 'FontSize', fs);
xlabel('Re (rad/s)'); ylabel('Im (rad/s)');
title('(b) Root locus in K_p', 'FontSize', fs, 'FontWeight', 'normal');
% --- step response
ax4 = axes('Position', [0.80 0.14 0.19 0.81]);
plot(A.step_t, A.step_y, '-', 'Color', [0 0.30 0.65], 'LineWidth', 1.2); hold on;
plot([0 1.5], [1 1], 'k:', 'LineWidth', 0.5);
xlim([0 1.5]); ylim([0 1.3]); grid on; set(gca, 'FontSize', fs);
xlabel('Time (s)'); ylabel('\theta / \theta_{ref}');
title('(c) Step response of T(s)', 'FontSize', fs, 'FontWeight', 'normal');
save_fig(fig, 'fig_design', 7.16, 2.5);
close(fig);
fprintf('\nSaved figures/fig_design.pdf and results/design_summary.txt\n');
