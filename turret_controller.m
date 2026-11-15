function [V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret)
%TURRET_CONTROLLER  Predictive fire-control law for the Guardian turret.
%
%   [V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, ...
%                                       theta_turret, theta_dot_turret)
%
%   Signal flow (one fixed design for every scenario, nothing is retuned):
%
%   1. TARGET-STATE ESTIMATION. A Kalman filter in Cartesian coordinates
%      with a Singer (exponentially correlated acceleration) model per axis,
%      state [x vx ax y vy ay]. Radar range/elevation is converted to (x, y)
%      and its covariance is rotated with the line of sight.
%   2. TARGET PREDICTION. The estimate is propagated open loop with the same
%      model. On every radar sample this issues a future trajectory valid
%      over the horizon t_p = 0.50 s; with no new sample (between samples or
%      during a measurement loss) the turret keeps following the trajectory
%      issued last. Because the loss interval never exceeds t_p, the turret
%      never runs past the end of an issued prediction.
%   3. PREDICTED ELEVATION REFERENCE. theta_ref(t) = atan2(y_hat, x_hat) of
%      the predicted position that is DUE at time t, limited to the turret
%      travel. (Laser time of flight is zero, so a pulse only hits if the
%      beam is on the target's present position; see the paper, Sec. II.)
%   4. FEEDBACK CONTROLLER. LTI compensator Gc(s) = Kp (s + z_i)/s * N(s),
%      N(s) a notch on the 56 rad/s electromechanical resonance. Type-2
%      loop: zero steady-state error for a constant line-of-sight rate.
%      Rate aiding adds the back-emf voltage of the predicted line-of-sight
%      rate outside the loop. Large-signal provisions: voltage limiter
%      upstream of the notch and conditional integration (anti-windup).
%   5. FIRE LOGIC. A pulse is requested when the estimated miss distance d
%      is small (aim gate), d + kappa*sigma_d is inside the engagement
%      radius (confidence gate, sigma_d from the filter covariance), the
%      beam is dwelling on the target rather than sweeping across it, and
%      the laser has recharged.
%
%   Auxiliary calls (not used by the evaluation loop):
%     C  = turret_controller('params')        controller parameter struct
%     S  = turret_controller('state')         internal state (for plotting)
%     turret_controller('override', struct)   override parameters (ablation)
%     turret_controller('override', [])       restore the default design
%
%   The controller keeps persistent state. A call with t = 0 starts a new
%   engagement. The file is self-contained and runs in MATLAB and Octave.

persistent C S OV
if isempty(C)
    C = ctrl_params(OV);
end
if ischar(t)
    switch lower(t)
        case 'params',   V_m = C;
        case 'state',    V_m = S;
        case 'override', OV = r_meas; C = ctrl_params(OV); S = []; V_m = C;
        otherwise, error('turret_controller:cmd', 'Unknown command "%s".', t);
    end
    fire_cmd = 0;
    return
end

% ---- start of a new engagement ------------------------------------------
if isempty(S) || t <= 1e-12 || t < S.t - 0.2
    S = init_state(t);
end
dt = t - S.t;
if dt < 0                       % solver stepped back slightly: hold outputs
    V_m = S.V; fire_cmd = S.fire;
    return
end
S.t = t;

% ---- 1. estimator: time update to the present instant --------------------
if S.init && dt > 0
    [F, Q] = singer_FQ(dt, C.alpha, C.q);
    S.x = [F*S.x(1:3); F*S.x(4:6)];
    Pxx = F*S.P(1:3,1:3)*F' + Q;
    Pyy = F*S.P(4:6,4:6)*F' + Q;
    Pxy = F*S.P(1:3,4:6)*F';
    S.P = [Pxx, Pxy; Pxy', Pyy];
end

% ---- 1. estimator: measurement update (only on a NEW radar sample) -------
% A sample is new only if its values differ from the last one used: two
% noisy radar samples are never bit-identical, whereas a caller that holds
% meas_valid = 1 between samples (or right after a radar loss) repeats the
% old values, and those must not be fused again as if they were current.
is_new = meas_valid && isfinite(r_meas) && isfinite(theta_meas) && r_meas > 0 && ...
         (~S.init || r_meas ~= S.z_last(1) || theta_meas ~= S.z_last(2));
if is_new
    c = cos(theta_meas);  s = sin(theta_meas);
    z = [r_meas*c; r_meas*s];
    Jc = [c, -r_meas*s; s, r_meas*c];
    R = Jc*[C.sigma_r^2, 0; 0, C.sigma_th^2]*Jc';
    if ~S.init
        S.x = [z(1); 0; 0; z(2); 0; 0];
        S.P = zeros(6);
        S.P([1 4], [1 4]) = R;
        S.P(2,2) = C.v0_var;      S.P(5,5) = C.v0_var;
        S.P(3,3) = C.sigma_m^2;   S.P(6,6) = C.sigma_m^2;
        S.init = true;
    else
        H = zeros(2, 6);  H(1,1) = 1;  H(2,4) = 1;
        Sinn = H*S.P*H' + R;
        K = (S.P*H')/Sinn;
        S.x = S.x + K*(z - H*S.x);
        IKH = eye(6) - K*H;
        S.P = IKH*S.P*IKH' + K*R*K';          % Joseph form
    end
    S.n_upd = S.n_upd + 1;
    S.t_upd = t;
    S.z_last = [r_meas, theta_meas];
end

if ~S.init                       % no target yet: hold position
    V_m = 0; fire_cmd = 0; S.V = 0; S.fire = 0;
    return
end

% ---- 2./3. prediction and elevation reference -----------------------------
% Present-time point of the issued trajectory (plus optional lead t_lead).
if C.t_lead > 0
    Fl = singer_FQ(C.t_lead, C.alpha, C.q);
    xr = Fl(1,:)*S.x(1:3);  yr = Fl(1,:)*S.x(4:6);
    vxr = Fl(2,:)*S.x(1:3); vyr = Fl(2,:)*S.x(4:6);
else
    xr = S.x(1);  yr = S.x(4);  vxr = S.x(2);  vyr = S.x(5);
end
theta_los = atan2(yr, xr);
theta_ref = min(max(theta_los, C.theta_min + C.theta_margin), C.theta_max - C.theta_margin);
clamped = (theta_ref ~= theta_los);
% End point of the issued trajectory: predicted state at the horizon t_p.
S.x_p = [C.F_tp(1,:)*S.x(1:3); C.F_tp(1,:)*S.x(4:6)];
S.theta_p = atan2(S.x_p(2), S.x_p(1));
S.theta_ref = theta_ref;

% ---- 4. feedback compensator Gc(s) = Kp (s + z_i)/s * N(s) ----------------
e = theta_ref - theta_turret;
% Rate aiding: back-emf voltage of the predicted line-of-sight rate,
% V_ff = theta_ref_dot/k_v. It enters outside the loop (stability margins
% are those of Gc*Gp) and is zero when the reference sits on a travel limit.
V_ff = 0;
if C.use_ff && ~clamped
    V_ff = C.k_ff*(xr*vyr - yr*vxr)/(xr^2 + yr^2)/C.k_v;
end
if dt > 0
    if dt ~= S.h                 % (re)build the Tustin matrices of the notch
        S.h = dt;
        M = eye(2) - C.An*dt/2;
        S.M1 = M\(eye(2) + C.An*dt/2);
        S.M2 = M\(C.Bn*dt/2);
    end
    % PI stage (trapezoidal). Anti-windup by conditional integration: the
    % integrator runs only inside the error band |e| <= e_int and while the
    % limiter is not active, so a large slew is flown without integrator
    % charge and arrives without the slow overshoot tail of the PI zero.
    % Inside the band (all of normal tracking) the loop is the LTI Gc(s).
    if abs(e) <= C.e_int && S.sat == 0
        S.xi = S.xi + C.Kp*C.z_i*dt*(e + S.e_prev)/2;
    end
    u_raw = C.Kp*e + S.xi + V_ff;
    % The voltage limiter acts on the PI output, UPSTREAM of the notch, so
    % that the notch also removes the resonance-frequency content of the
    % clipped command (a clipped voltage applied directly would ring the
    % 56 rad/s mode, which the cancellation leaves undamped).
    u_pi = min(max(u_raw, -C.V_lim), C.V_lim);
    S.sat = sign(u_raw)*(u_pi ~= u_raw);
    % notch stage (Tustin)
    S.xn = S.M1*S.xn + S.M2*(u_pi + S.u_prev);
else
    u_pi = min(max(C.Kp*e + S.xi + V_ff, -C.V_lim), C.V_lim);
end
v = C.Cn*S.xn + u_pi;            % notch has unity direct gain
V_m = min(max(v, -C.V_max), C.V_max);          % hard limit of the drive
S.e_prev = e;
S.u_prev = u_pi;

% ---- 5. fire logic --------------------------------------------------------
ct = cos(theta_turret);  st = sin(theta_turret);
xh = S.x(1);  yh = S.x(4);
d_hat = -xh*st + yh*ct;                       % predicted miss (perpendicular to beam)
along = xh*ct + yh*st;                        % target must be in front of the turret
hrow = [-st, ct];
sig_d = sqrt(max(hrow*S.P([1 4], [1 4])*hrow', 0));
% rate of change of the miss: target cross-velocity minus beam sweep velocity
d_dot = (-S.x(2)*st + S.x(5)*ct) - along*theta_dot_turret;
% worst miss over a short look-ahead tau_fire: the beam must be dwelling on
% the target, not sweeping across it, when the pulse is released
d_gate = max(abs(d_hat), abs(d_hat + d_dot*C.tau_fire));
S.d_hat = d_hat;  S.sig_d = sig_d;  S.d_gate = d_gate;

if S.fire && (t - S.t_fire) >= C.fire_hold
    S.fire = 0;                               % release, so that a new 0->1 edge is possible
end
ready = (S.n_upd >= C.fire_min_updates) && ((t - S.t_fire) >= C.t_recharge + 1e-9) && along > 0;
if ~S.fire && ready && (d_gate <= C.d_aim) && (d_gate + C.kappa*sig_d <= C.R_fire)
    S.fire = 1;
    S.t_fire = t;
end
fire_cmd = S.fire;
S.V = V_m;
end

% =========================================================================
function S = init_state(t)
S.t = t;
S.init = false;  S.n_upd = 0;  S.t_upd = -inf;  S.z_last = [NaN, NaN];
S.x = zeros(6, 1);  S.P = zeros(6);
S.x_p = [0; 0];  S.theta_p = 0;  S.theta_ref = 0;
S.xi = 0;  S.xn = [0; 0];  S.e_prev = 0;  S.u_prev = 0;  S.sat = 0;
S.h = -1;  S.M1 = eye(2);  S.M2 = [0; 0];
S.fire = 0;  S.t_fire = -inf;  S.V = 0;
S.d_hat = NaN;  S.sig_d = NaN;  S.d_gate = NaN;
end

% =========================================================================
function C = ctrl_params(OV)
% ---- radar and engagement data (from the problem statement) --------------
C.Ts        = 0.05;                 % radar sampling interval [s]
C.sigma_r   = 1.0;                  % range noise [m]
C.sigma_th  = 0.02*pi/180;          % elevation noise [rad]
C.t_p       = 0.50;                 % prediction horizon [s]
C.R_eng     = 2.0;                  % engagement radius [m]
C.t_recharge = 1.0;                 % minimum inter-pulse interval [s]
C.theta_min = -10*pi/180;           % elevation stops [rad]
C.theta_max =  85*pi/180;
a_max = 20;  v_max = 150;           % target envelope

% ---- estimator (Singer model) --------------------------------------------
% sigma_m: per-axis acceleration std of a vector uniformly distributed in a
%          disc of radius a_max  ->  a_max/2
% tau_m  : maneuver correlation time; the weave periods that fit inside the
%          envelope are a few seconds
C.sigma_m = a_max/2;                % [m/s^2]
C.tau_m   = 2.0;                    % [s]
C.v0_var  = v_max^2/2;              % prior velocity variance per axis [m^2/s^2]

% ---- reference generation -------------------------------------------------
C.t_lead       = 0;                 % extra lead of the reference [s] (0: due-time reference)
C.theta_margin = 0.5*pi/180;        % keep the reference clear of the stops [rad]

% ---- plant model used for the compensator --------------------------------
N = 100; Kt = 10; Ke = 10; Rm = 1.0; Lm = 0.10; b = 1500;
Jeq = 1200 + N^2*0.20;
den = [Jeq*Lm, Jeq*Rm + b*Lm, b*Rm + N^2*Kt*Ke];
C.wn_p   = sqrt(den(3)/den(1));               % 55.94 rad/s
C.zeta_p = den(2)/(2*sqrt(den(1)*den(3)));    % 0.0936
C.k_v    = N*Kt/den(3);                       % 9.985e-4 rad/s per V

% ---- feedback compensator  Gc(s) = Kp (s + z_i)/s * N(s) ------------------
C.w_c    = 15;                      % design crossover [rad/s]
C.z_i    = 3;                       % PI zero [rad/s] (= w_c/5)
C.zeta_d = 0.70;                    % damping of the notch poles
C.Kp     = [];                      % set below for unity loop gain at w_c
C.V_max  = 1000;                    % drive voltage limit [V] (design assumption)
C.V_lim  = 900;                     % limiter on the PI output, upstream of the notch [V]
C.e_int  = 2e-3;                    % integrator active for |e| <= e_int [rad]
C.use_ff = true;                    % rate aiding V_ff = theta_ref_dot/k_v
C.k_ff   = 1.0;                     % rate-aiding gain (1 = exact back-emf compensation)

% ---- fire logic -----------------------------------------------------------
C.kappa  = 3.0;                     % confidence factor on the miss estimate
C.R_fire = C.R_eng;                 % confidence gate: d_gate + kappa*sigma_d <= R_fire
C.d_aim  = 0.5;                     % aim gate: estimated miss must be below d_aim [m]
C.tau_fire = 0.02;                  % look-ahead over which the aim must hold [s]
C.fire_min_updates = 4;             % radar samples needed before firing
C.fire_hold = 0.10;                 % duration of the fire_cmd = 1 pulse [s]

% ---- overrides (ablation studies only) ------------------------------------
if nargin > 0 && isstruct(OV)
    fn = fieldnames(OV);
    for k = 1:numel(fn), C.(fn{k}) = OV.(fn{k}); end
end

% ---- derived quantities ----------------------------------------------------
C.alpha = 1/C.tau_m;
C.q     = 2*C.alpha*C.sigma_m^2;              % PSD of the acceleration driving noise
C.F_tp  = singer_FQ(C.t_p, C.alpha, C.q);
% notch N(s) = (s^2 + 2 zeta_p wn s + wn^2)/(s^2 + 2 zeta_d wn s + wn^2)
C.notch_num = [1, 2*C.zeta_p*C.wn_p, C.wn_p^2];
C.notch_den = [1, 2*C.zeta_d*C.wn_p, C.wn_p^2];
C.An = [0, 1; -C.notch_den(3), -C.notch_den(2)];
C.Bn = [0; 1];
C.Cn = [0, C.notch_num(2) - C.notch_den(2)];
if isempty(C.Kp)
    % |L(j w_c)| = 1 with L = Kp (s+z_i)/s * k_v wn^2 / (s (s^2 + 2 zeta_d wn s + wn^2))
    w = C.w_c;  s = 1i*w;
    Lw = (s + C.z_i)/s * C.k_v*C.wn_p^2/(s*(s^2 + C.notch_den(2)*s + C.notch_den(3)));
    C.Kp = 1/abs(Lw);
end
C.Ki = C.Kp*C.z_i;
C.Ka = C.Ki*C.k_v;                            % acceleration error constant [1/s^2]
% compensator polynomials, Gc(s) = Gc_num/Gc_den
C.Gc_num = C.Kp*conv([1, C.z_i], C.notch_num);
C.Gc_den = conv([1, 0], C.notch_den);
end

% =========================================================================
function [F, Q] = singer_FQ(T, al, q)
% Transition matrix and process-noise covariance of the Singer model
%   d/dt [p; v; a] = [0 1 0; 0 0 1; 0 0 -al] [p; v; a] + [0; 0; 1] w,  E[w w'] = q
aT = al*T;  e1 = exp(-aT);
F = [1, T, (aT - 1 + e1)/al^2;
     0, 1, (1 - e1)/al;
     0, 0, e1];
if nargout < 2, return, end
if aT < 0.02
    % limit al*T -> 0 (white-noise jerk); avoids cancellation in the closed form
    Q = q*[T^5/20, T^4/8, T^3/6;
           T^4/8,  T^3/3, T^2/2;
           T^3/6,  T^2/2, T];
else
    e2 = exp(-2*aT);
    q11 = (1 - e2 + 2*aT + 2*aT^3/3 - 2*aT^2 - 4*aT*e1)/(2*al^5);
    q12 = (e2 + 1 - 2*e1 + 2*aT*e1 - 2*aT + aT^2)/(2*al^4);
    q13 = (1 - e2 - 2*aT*e1)/(2*al^3);
    q22 = (4*e1 - 3 - e2 + 2*aT)/(2*al^3);
    q23 = (e2 + 1 - 2*e1)/(2*al^2);
    q33 = (1 - e2)/(2*al);
    Q = q*[q11, q12, q13; q12, q22, q23; q13, q23, q33];
end
end
