function S = scenario_definition(scenario)
%SCENARIO_DEFINITION  Parameters of a standardized or randomized engagement.
%
%   S = scenario_definition(scenario)
%
%   scenario = 'A'  nominal, non-maneuvering target
%            = 'B'  smooth bounded maneuver in x and y starting at t = 1 s
%            = 'C'  scenario B motion + radar loss of 0.50 s starting at 1.00 s
%            = n    (numeric, n >= 1) randomized evaluation case n: random
%                   initial geometry, random bounded maneuver in x and y and
%                   one or two random measurement-loss intervals
%            = struct('id', <any of the above>, 'seed', k)  same motion, but
%                   measurement-noise realisation k (for Monte-Carlo runs)
%
%   Target motion (per axis, tau = t - t_m >= 0):
%       a(t) = sum_j A_j sin(w_j tau)                 (zero and continuous at onset)
%       v(t) = v0 + sum_j A_j/w_j (1 - cos(w_j tau))
%       p(t) = p0 + v0 t + sum_j A_j/w_j (tau - sin(w_j tau)/w_j)
%
%   Every scenario is verified against the prescribed envelope
%       800 m <= r <= 2000 m,  |v| <= 150 m/s,  |a| <= 20 m/s^2   for 0 <= t <= 10 s.
%   Randomized cases are drawn by rejection sampling until the envelope holds.
%
%   A portable Park-Miller generator is used, so scenarios are bit-identical
%   in MATLAB and Octave.

persistent cache_key cache_val
noise_seed = [];
id = scenario;
if isstruct(scenario)
    id = scenario.id;
    if isfield(scenario, 'seed'), noise_seed = scenario.seed; end
end
if ischar(id), key = upper(id); else, key = sprintf('N%d', id); end
fullkey = sprintf('%s#%d', key, sum([noise_seed, -1]));
if ~isempty(cache_key) && strcmp(cache_key, fullkey)
    S = cache_val;
    return
end

L = scenario_limits();

if ischar(id)
    S.x0 = 1850;  S.y0 = 600;       % initial position [m]
    S.vx0 = -100; S.vy0 = -5;       % initial velocity [m/s]
    S.t_m = 1.0;                    % maneuver onset [s]
    S.Ax = []; S.wx = []; S.Ay = []; S.wy = [];
    S.blackout = zeros(0, 2);       % rows: [start, duration]
    switch upper(id)
        case 'A'
            S.name = 'A'; S.seed = 11;
        case 'B'
            S.name = 'B'; S.seed = 12;
            S.Ax = 12;  S.wx = 1.3;
            S.Ay = 15;  S.wy = 2.0;
        case 'C'
            S.name = 'C'; S.seed = 13;
            S.Ax = 12;  S.wx = 1.3;
            S.Ay = 15;  S.wy = 2.0;
            S.blackout = [1.00, 0.50];
        otherwise
            error('scenario_definition:unknown', 'Unknown scenario "%s".', id);
    end
    [ok, msg] = check_envelope(S, L);
    if ~ok, error('scenario_definition:envelope', 'Scenario %s: %s', S.name, msg); end
else
    S = random_case(id, L);
end
if ~isempty(noise_seed), S.seed = noise_seed; end

cache_key = fullkey;
cache_val = S;
end

% -------------------------------------------------------------------------
function L = scenario_limits()
L.r_min = 800;  L.r_max = 2000;     % operating range [m]
L.v_max = 150;                      % speed limit [m/s]
L.a_max = 20;                       % acceleration limit [m/s^2]
L.t_end = 10;                       % engagement duration [s]
L.el_min = -8*pi/180;               % keep the LOS inside the turret travel
L.el_max = 83*pi/180;
end

% -------------------------------------------------------------------------
function S = random_case(n, L)
s = 1000 + 7919*n;                  % generator state
for attempt = 1:5000
    [u, s] = lcg(s, 20);
    r0  = 1450 + 500*u(1);
    el0 = (8 + 37*u(2))*pi/180;
    spd = 70 + 65*u(3);
    hdg = el0 + pi + (70*u(4) - 35)*pi/180;      % roughly inbound
    S.x0 = r0*cos(el0);  S.y0 = r0*sin(el0);
    S.vx0 = spd*cos(hdg); S.vy0 = spd*sin(hdg);
    S.t_m = 2*u(5);
    S.Ax = 2*u(6:7)' - 1;   S.wx = 0.6 + 2.4*u(8:9)';
    S.Ay = 2*u(10:11)' - 1; S.wy = 0.6 + 2.4*u(12:13)';
    % scale the maneuver so that its peak acceleration is 60..98 % of a_max
    tt = linspace(0, L.t_end, 2001);
    [~, ~, ~, ~, ax, ay] = motion(S, tt);
    apk = max(sqrt(ax.^2 + ay.^2));
    k = (0.60 + 0.38*u(14))*L.a_max/apk;
    S.Ax = k*S.Ax;  S.Ay = k*S.Ay;
    % one or two measurement-loss intervals (onset >= 1 s, duration <= 0.5 s)
    b1 = 1.0 + 3.5*u(15);  d1 = 0.20 + 0.30*u(16);
    S.blackout = [b1, d1];
    if u(17) > 0.5
        b2 = b1 + d1 + 1.5 + 2.5*u(18);  d2 = 0.20 + 0.30*u(19);
        if b2 + d2 < L.t_end, S.blackout = [S.blackout; b2, d2]; end
    end
    S.name = sprintf('E%02d', n);
    S.seed = 5000 + n;
    if check_envelope(S, L), return, end
end
error('scenario_definition:sampling', 'No admissible trajectory found for case %d.', n);
end

% -------------------------------------------------------------------------
function [ok, msg] = check_envelope(S, L)
tt = linspace(0, L.t_end, 5001);
[x, y, vx, vy, ax, ay] = motion(S, tt);
r = sqrt(x.^2 + y.^2);  el = atan2(y, x);
ok = true; msg = '';
if min(r) < L.r_min || max(r) > L.r_max, ok = false; msg = 'range limits violated'; return, end
if max(sqrt(vx.^2 + vy.^2)) > L.v_max, ok = false; msg = 'speed limit violated'; return, end
if max(sqrt(ax.^2 + ay.^2)) > L.a_max, ok = false; msg = 'acceleration limit violated'; return, end
if min(el) < L.el_min || max(el) > L.el_max, ok = false; msg = 'elevation outside turret travel'; return, end
end

% -------------------------------------------------------------------------
function [x, y, vx, vy, ax, ay] = motion(S, t)
tau = max(t - S.t_m, 0);
[dx, dvx, dax] = axis_maneuver(S.Ax, S.wx, tau);
[dy, dvy, day] = axis_maneuver(S.Ay, S.wy, tau);
x = S.x0 + S.vx0*t + dx;   vx = S.vx0 + dvx;   ax = dax;
y = S.y0 + S.vy0*t + dy;   vy = S.vy0 + dvy;   ay = day;
end

function [p, v, a] = axis_maneuver(A, w, tau)
p = zeros(size(tau)); v = p; a = p;
for j = 1:numel(A)
    a = a + A(j)*sin(w(j)*tau);
    v = v + A(j)/w(j)*(1 - cos(w(j)*tau));
    p = p + A(j)/w(j)*(tau - sin(w(j)*tau)/w(j));
end
end

% -------------------------------------------------------------------------
function [u, s] = lcg(s, n)
% Park-Miller minimal-standard generator, exact in double precision.
u = zeros(n, 1);
for k = 1:n
    s = mod(16807*s, 2147483647);
    u(k) = s/2147483647;
end
end
