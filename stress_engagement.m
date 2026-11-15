function R = stress_engagement(T, opts)
%STRESS_ENGAGEMENT  Independent evaluator for the stress tests.
%
%   R = stress_engagement(T, opts)
%
%   Written separately from simulate_engagement.m, in the way an external
%   evaluation script might be written. It shares no code with it: it has
%   its own radar model (MATLAB randn), its own integration loop and its own
%   hit test, and it uses only turret_controller.m and turret_plant.m.
%
%   T     trajectory from stress_trajectory.m
%   opts  .dt          step [s], a multiple of 1 ms             (default 1e-3)
%         .integrator  'rk4'      plant RK4, voltage held        (default)
%                      'euler'    forward Euler
%                      'rk4stage' controller called at every RK4 stage
%         .meas        'pulse'    meas_valid = 1 only on a new sample (default)
%                      'held'     meas_valid = 1 whenever the radar sees the
%                                 target, measurement held between samples
%                      'nan'      as 'pulse', but r_meas = theta_meas = NaN
%                                 during a radar loss
%         .theta0      initial turret elevation [rad]            (default 0)
%         .phase       time of the first radar sample [s]        (default 0)
%         .t_first     time of the first controller call [s]     (default 0)
%         .seed        noise seed                                (default T.seed)

if nargin < 2, opts = struct(); end
dt     = getopt(opts, 'dt', 1e-3);
integ  = getopt(opts, 'integrator', 'rk4');
mmode  = getopt(opts, 'meas', 'pulse');
th     = getopt(opts, 'theta0', 0);
phase  = getopt(opts, 'phase', 0);
tfirst = getopt(opts, 't_first', 0);
seed   = getopt(opts, 'seed', T.seed);

Ts = 0.05;  sig_r = 1.0;  sig_th = 0.02*pi/180;  R_eng = 2.0;  t_rech = 1.0;
rs = RandStream('mt19937ar', 'Seed', seed);
step = round(dt/1e-3);
idx = (1 + round(tfirst/1e-3)):step:10001;       % truth samples used (0..10 s)
n = numel(idx);

om = 0;  cur = 0;
rm = 0;  thm = 0;  seen = false;
next_sample = phase;
fire_prev = 0;  t_last = -inf;
pulse_t = [];  pulse_miss = [];
miss = zeros(n, 1);  tt = T.t(idx);  Vlog = zeros(n, 1);

for k = 1:n
    tk = tt(k);  i = idx(k);
    lost = any(tk >= T.blackout(:,1) - 1e-9 & tk < T.blackout(:,1) + T.blackout(:,2) - 1e-9);
    new = 0;
    if tk >= next_sample - 1e-9
        next_sample = next_sample + Ts;
        if ~lost
            rm  = hypot(T.x(i), T.y(i)) + sig_r*randn(rs);
            thm = atan2(T.y(i), T.x(i)) + sig_th*randn(rs);
            new = 1;  seen = true;
        end
    end
    switch mmode
        case 'pulse', valid = new;  r_in = rm;  th_in = thm;
        case 'held',  valid = double(seen && ~lost);  r_in = rm;  th_in = thm;
        case 'nan',   valid = new;  r_in = rm;  th_in = thm;
                      if lost, r_in = NaN;  th_in = NaN; end
        otherwise, error('unknown meas mode');
    end

    [V, fire] = turret_controller(tk, r_in, th_in, valid, th, om);
    Vlog(k) = V;

    % hit test: distance from the target to the beam ray
    along = T.x(i)*cos(th) + T.y(i)*sin(th);
    perp  = -T.x(i)*sin(th) + T.y(i)*cos(th);
    if along > 0, miss(k) = abs(perp); else, miss(k) = hypot(T.x(i), T.y(i)); end
    if fire == 1 && fire_prev == 0 && (tk - t_last) >= t_rech - 1e-9
        t_last = tk;
        pulse_t(end+1, 1) = tk; %#ok<AGROW>
        pulse_miss(end+1, 1) = miss(k); %#ok<AGROW>
    end
    fire_prev = fire;

    if k == n, break, end
    switch integ
        case 'euler'
            [a, b, c] = turret_plant(tk, V, th, om, cur);
            th = th + dt*a;  om = om + dt*b;  cur = cur + dt*c;
        case 'rk4'
            [a1, b1, c1] = turret_plant(tk,        V, th,            om,            cur);
            [a2, b2, c2] = turret_plant(tk + dt/2, V, th + dt/2*a1,  om + dt/2*b1,  cur + dt/2*c1);
            [a3, b3, c3] = turret_plant(tk + dt/2, V, th + dt/2*a2,  om + dt/2*b2,  cur + dt/2*c2);
            [a4, b4, c4] = turret_plant(tk + dt,   V, th + dt*a3,    om + dt*b3,    cur + dt*c3);
            th  = th  + dt/6*(a1 + 2*a2 + 2*a3 + a4);
            om  = om  + dt/6*(b1 + 2*b2 + 2*b3 + b4);
            cur = cur + dt/6*(c1 + 2*c2 + 2*c3 + c4);
        case 'rk4stage'
            % the controller is re-evaluated inside every stage
            vh = valid*strcmp(mmode, 'held');
            [a1, b1, c1] = turret_plant(tk, V, th, om, cur);
            V2 = turret_controller(tk + dt/2, r_in, th_in, vh, th + dt/2*a1, om + dt/2*b1);
            [a2, b2, c2] = turret_plant(tk + dt/2, V2, th + dt/2*a1, om + dt/2*b1, cur + dt/2*c1);
            V3 = turret_controller(tk + dt/2, r_in, th_in, vh, th + dt/2*a2, om + dt/2*b2);
            [a3, b3, c3] = turret_plant(tk + dt/2, V3, th + dt/2*a2, om + dt/2*b2, cur + dt/2*c2);
            V4 = turret_controller(tk + dt, r_in, th_in, vh, th + dt*a3, om + dt*b3);
            [a4, b4, c4] = turret_plant(tk + dt, V4, th + dt*a3, om + dt*b3, cur + dt*c3);
            th  = th  + dt/6*(a1 + 2*a2 + 2*a3 + a4);
            om  = om  + dt/6*(b1 + 2*b2 + 2*b3 + b4);
            cur = cur + dt/6*(c1 + 2*c2 + 2*c3 + c4);
        otherwise, error('unknown integrator');
    end
end

hit = pulse_miss <= R_eng;
R.name = T.name;
R.n_pulses = numel(pulse_t);  R.n_hits = sum(hit);
R.success = any(hit);
R.t_intercept = NaN;  if R.success, R.t_intercept = pulse_t(find(hit, 1)); end
R.pulse_t = pulse_t;  R.pulse_miss = pulse_miss;
trk = tt >= 2;
R.miss_rms = sqrt(mean(miss(trk).^2));  R.miss_max = max(miss(trk));
R.in_radius = mean(miss(trk) <= R_eng);
R.V_peak = max(abs(Vlog));
R.finite = all(isfinite(Vlog)) && isfinite(th);
el = atan2(T.y(idx), T.x(idx));
R.los_rate_max = max(abs(diff(el)))/dt;
end

function v = getopt(o, f, d)
if isfield(o, f), v = o.(f); else, v = d; end
end
