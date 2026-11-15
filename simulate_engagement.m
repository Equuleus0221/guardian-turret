function L = simulate_engagement(scenario, opts)
%SIMULATE_ENGAGEMENT  Closed-loop simulation of one 10 s engagement.
%
%   L = simulate_engagement(scenario)
%   L = simulate_engagement(scenario, opts)
%
%   Uses exactly the three prescribed functions
%       target_motion_generator, turret_controller, turret_plant
%   and integrates the plant with a fixed-step RK4 scheme (the motor voltage
%   is held constant over each step, i.e. the controller runs at 1/dt).
%
%   opts.dt      integration / controller step [s]      (default 1e-3)
%   opts.t_end   engagement duration [s]                (default 10)
%   opts.theta0  initial turret elevation [rad]         (default 0)
%
%   L is a log struct (time histories, pulse list and the metrics L.M).
%   Laser rule: a pulse is issued on a 0->1 transition of fire_cmd, provided
%   at least 1.0 s has elapsed since the previous issued pulse. A pulse is a
%   hit when the distance between the beam (a ray along the turret
%   elevation) and the true target position is <= R_eng = 2 m.

if nargin < 2, opts = struct(); end
dt     = getopt(opts, 'dt', 1e-3);
t_end  = getopt(opts, 't_end', 10);
theta0 = getopt(opts, 'theta0', 0);
R_eng = 2.0;  t_recharge = 1.0;

C  = turret_controller('params');
n  = round(t_end/dt) + 1;
t  = (0:n-1)'*dt;
[x, y, vx, vy, ax, ay, lost] = target_truth(t, scenario);

z = zeros(n, 1);
L.t = t;  L.x = x;  L.y = y;  L.vx = vx;  L.vy = vy;  L.ax = ax;  L.ay = ay;  L.lost = lost;
L.r = sqrt(x.^2 + y.^2);  L.theta = atan2(y, x);
L.r_meas = z;  L.theta_meas = z;  L.meas_valid = z;
L.theta_t = z;  L.omega_t = z;  L.i = z;  L.V = z;  L.fire = z;
L.xh = nan(n, 1);  L.yh = nan(n, 1);  L.vxh = nan(n, 1);  L.vyh = nan(n, 1);
L.axh = nan(n, 1); L.ayh = nan(n, 1);
L.theta_ref = nan(n, 1);  L.theta_p = nan(n, 1);  L.xp = nan(n, 1);  L.yp = nan(n, 1);
L.d_hat = nan(n, 1);  L.sig_d = nan(n, 1);
pulse_t = [];  pulse_miss = [];

th = theta0;  om = 0;  cur = 0;
fire_prev = 0;  t_last_pulse = -inf;
for k = 1:n
    tk = t(k);
    [rm, thm, valid] = target_motion_generator(tk, scenario);
    [V, fire] = turret_controller(tk, rm, thm, valid, th, om);
    S = turret_controller('state');

    L.r_meas(k) = rm;  L.theta_meas(k) = thm;  L.meas_valid(k) = valid;
    L.theta_t(k) = th;  L.omega_t(k) = om;  L.i(k) = cur;  L.V(k) = V;  L.fire(k) = fire;
    if S.init
        L.xh(k) = S.x(1);  L.vxh(k) = S.x(2);  L.axh(k) = S.x(3);
        L.yh(k) = S.x(4);  L.vyh(k) = S.x(5);  L.ayh(k) = S.x(6);
        L.theta_ref(k) = S.theta_ref;  L.theta_p(k) = S.theta_p;
        L.xp(k) = S.x_p(1);  L.yp(k) = S.x_p(2);
        L.d_hat(k) = S.d_hat;  L.sig_d(k) = S.sig_d;
    end

    % laser pulse on a rising edge, subject to the recharge interval
    if fire && ~fire_prev && (tk - t_last_pulse) >= t_recharge - 1e-9
        t_last_pulse = tk;
        pulse_t(end+1, 1) = tk; %#ok<AGROW>
        pulse_miss(end+1, 1) = beam_miss(x(k), y(k), th); %#ok<AGROW>
    end
    fire_prev = fire;

    % RK4 step of the plant with the voltage held
    if k < n
        [k1a, k1b, k1c] = turret_plant(tk,        V, th,             om,             cur);
        [k2a, k2b, k2c] = turret_plant(tk + dt/2, V, th + dt/2*k1a,  om + dt/2*k1b,  cur + dt/2*k1c);
        [k3a, k3b, k3c] = turret_plant(tk + dt/2, V, th + dt/2*k2a,  om + dt/2*k2b,  cur + dt/2*k2c);
        [k4a, k4b, k4c] = turret_plant(tk + dt,   V, th + dt*k3a,    om + dt*k3b,    cur + dt*k3c);
        th  = th  + dt/6*(k1a + 2*k2a + 2*k3a + k4a);
        om  = om  + dt/6*(k1b + 2*k2b + 2*k3b + k4b);
        cur = cur + dt/6*(k1c + 2*k2c + 2*k3c + k4c);
    end
end

% ---- derived signals -------------------------------------------------------
L.e_point = L.theta_t - L.theta;                          % pointing error [rad]
L.miss    = -x.*sin(L.theta_t) + y.*cos(L.theta_t);       % signed beam-target miss [m]
L.e_ref   = L.theta_ref - L.theta_t;                      % servo error [rad]
L.e_est   = L.theta_ref - L.theta;                        % reference (estimation) error [rad]
L.pos_err = sqrt((L.xh - x).^2 + (L.yh - y).^2);          % position estimation error [m]
L.cross_err = (-(L.xh - x).*y + (L.yh - y).*x)./L.r;      % cross-range estimation error [m]
% prediction made at t for t + t_p, compared with the truth at t + t_p
[xf, yf] = target_truth(t + C.t_p, scenario);
L.pred_err = sqrt((L.xp - xf).^2 + (L.yp - yf).^2);       % [m]
L.pred_cross = (-(L.xp - xf).*yf + (L.yp - yf).*xf)./sqrt(xf.^2 + yf.^2);
L.pulse_t = pulse_t;  L.pulse_miss = pulse_miss;  L.pulse_hit = pulse_miss <= R_eng;

% ---- metrics -----------------------------------------------------------------
M.n_pulses = numel(pulse_t);
M.n_hits   = sum(L.pulse_hit);
M.success  = M.n_hits > 0;
if M.success, M.t_intercept = pulse_t(find(L.pulse_hit, 1)); else, M.t_intercept = NaN; end
M.pulse_miss_mean = mean(pulse_miss);
M.pulse_miss_max  = max([pulse_miss; NaN]);
trk = t >= 2.0;                                           % after acquisition
M.miss_rms   = sqrt(mean(L.miss(trk).^2));
M.miss_max   = max(abs(L.miss(trk)));
M.point_rms_mrad = 1e3*sqrt(mean(L.e_point(trk).^2));
M.servo_rms_mrad = 1e3*sqrt(mean(L.e_ref(trk).^2));
M.est_rms_mrad   = 1e3*sqrt(mean(L.e_est(trk).^2));
M.cross_rms  = sqrt(mean(L.cross_err(trk).^2));
M.pos_rms    = sqrt(mean(L.pos_err(trk).^2));
pv = trk & t <= t_end - C.t_p;
M.pred_rms   = sqrt(mean(L.pred_err(pv).^2));
M.pred_cross_rms = sqrt(mean(L.pred_cross(pv).^2));
M.in_radius  = mean(abs(L.miss(trk)) <= R_eng);           % fraction of time on target
M.V_peak = max(abs(L.V));  M.i_peak = max(abs(L.i));
M.V_rms_track = sqrt(mean(L.V(trk).^2));
L.M = M;
L.name = scenario_name(scenario);
end

% ----------------------------------------------------------------------------
function d = beam_miss(x, y, th)
% minimum distance between the beam (ray from the origin) and the target
along = x*cos(th) + y*sin(th);
if along > 0
    d = abs(-x*sin(th) + y*cos(th));
else
    d = sqrt(x^2 + y^2);
end
end

function v = getopt(o, f, d)
if isfield(o, f), v = o.(f); else, v = d; end
end

function s = scenario_name(scenario)
S = scenario_definition(scenario);
s = S.name;
end
