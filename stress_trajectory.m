function T = stress_trajectory(family, k)
%STRESS_TRAJECTORY  Target trajectories for the stress tests.
%
%   T = stress_trajectory(family, k)   k-th random member of a family
%
%   These trajectories are deliberately NOT of the sinusoidal type used in
%   scenario_definition.m. The target is flown in its own path frame with a
%   normal (turn) and a tangential (speed) acceleration command that is
%   piecewise constant and passed through a first-order lag, so that the
%   acceleration is continuous.
%
%   family = 'straight'  non-maneuvering, any heading (also receding/diving)
%            'jink'      bang-bang turns of 12-19.5 m/s^2 with speed changes
%            'turn'      sustained 17-19.5 m/s^2 turn with one or two reversals
%            'reversal'  full turn that reverses exactly when the radar is
%                        lost for 0.5 s (worst case for prediction)
%            'crosser'   150 m/s pass at 810-900 m: highest line-of-sight rate
%            'highel'    engagement at 60-84 deg elevation, near the upper stop
%            'multiloss' jinking target with four or five radar losses
%
%   Every trajectory satisfies 800 m <= r <= 2000 m, |v| <= 150 m/s and
%   |a| <= 20 m/s^2 on 0 <= t <= 10 s (rejection sampling), and stays inside
%   the turret travel.
%
%   T: struct with t (1 ms grid, 0..10.6 s), x, y, vx, vy, ax, ay,
%      blackout (rows [start, duration]) and name.

fams = {'straight', 'jink', 'turn', 'reversal', 'crosser', 'highel', 'multiloss'};
fi = find(strcmp(fams, family));
if isempty(fi), error('stress_trajectory:family', 'Unknown family "%s".', family); end
rs = RandStream('mt19937ar', 'Seed', 10000*fi + k);
U = @(a, b) a + (b - a)*rand(rs);

dt = 1e-3;  t = (0:dt:10.6)';  n = numel(t);
tau_a = 0.10;                       % lag of the acceleration command [s]
a_cap = 19.5;  v_cap = 147;
el_lim = [-8, 83]*pi/180;
if strcmp(family, 'highel'), el_lim = [-8, 84.2]*pi/180; end

for attempt = 1:5000
    % ---- initial state -----------------------------------------------------
    r0 = U(1300, 1950);  el0 = U(5, 60)*pi/180;  spd = U(60, 145);
    hdg = el0 + pi + U(-60, 60)*pi/180;
    if strcmp(family, 'straight'), hdg = U(0, 2*pi); end
    if strcmp(family, 'highel'),   el0 = U(60, 80)*pi/180;  r0 = U(1100, 1900); end
    p = r0*[cos(el0), sin(el0)];  v = spd*[cos(hdg), sin(hdg)];

    % ---- acceleration schedule: rows [t_start, a_normal, a_tangential] -----
    sched = [0, 0, 0];
    blackout = zeros(0, 2);
    switch family
        case 'straight'
        case {'jink', 'multiloss', 'highel'}
            ts = U(0, 1.5);  sgn = sign(U(-1, 1));
            while ts < 10.6
                amp = U(12, 19.5);
                if strcmp(family, 'highel'), amp = U(6, 15); end
                sched(end+1, :) = [ts, sgn*amp, U(-5, 5)]; %#ok<AGROW>
                ts = ts + U(0.6, 2.0);  sgn = -sgn;
            end
        case 'turn'
            sgn = sign(U(-1, 1));  amp = U(17, 19.5);
            sched(end+1, :) = [U(0, 1), sgn*amp, 0];
            trev = sort([U(2, 5), U(5, 9)]);
            nrev = 1 + (rand(rs) > 0.5);
            for j = 1:nrev
                sgn = -sgn;  sched(end+1, :) = [trev(j), sgn*amp, 0]; %#ok<AGROW>
            end
        case 'reversal'
            b0 = U(1.5, 8);  sgn = sign(U(-1, 1));
            sched(end+1, :) = [b0 - U(1.0, 1.5), sgn*a_cap, 0];
            sched(end+1, :) = [b0, -sgn*a_cap, 0];
            sched(end+1, :) = [b0 + 1.5, 0, 0];
            blackout = [b0, 0.5];
        case 'crosser'
            rc = U(810, 900);  tc = U(4, 6);  phic = U(34, 42)*pi/180;
            spd = U(135, 149);  sgn = sign(U(-1, 1));
            v = sgn*spd*[-sin(phic), cos(phic)];
            p = rc*[cos(phic), sin(phic)] - v*tc;
    end
    % ---- radar-loss intervals ----------------------------------------------
    if strcmp(family, 'multiloss')
        b = 1.0 + U(0, 0.5);
        for j = 1:(4 + (rand(rs) > 0.5))
            d = U(0.3, 0.5);
            if b + d < 9.8, blackout(end+1, :) = [b, d]; end %#ok<AGROW>
            b = b + d + U(0.25, 0.8);
        end
    elseif isempty(blackout)
        b = U(1.0, 7.5);
        if rand(rs) < 0.3, b = 1.0; end          % loss at the earliest allowed instant
        blackout = [b, U(0.3, 0.5)];
        if rand(rs) < 0.5
            b2 = b + blackout(1, 2) + U(0.3, 2.0);
            if b2 < 9.3, blackout(2, :) = [b2, U(0.3, 0.5)]; end
        end
    end

    % ---- fly the trajectory --------------------------------------------------
    X = zeros(n, 6);  an = 0;  at = 0;  js = 1;
    for i = 1:n
        while js < size(sched, 1) && t(i) >= sched(js+1, 1), js = js + 1; end
        an = an + dt/tau_a*(sched(js, 2) - an);
        at = at + dt/tau_a*(sched(js, 3) - at);
        s = norm(v);  vh = v/s;  nh = [-vh(2), vh(1)];
        att = at;
        if (s > v_cap && att > 0) || (s < 45 && att < 0), att = 0; end
        a = att*vh + an*nh;
        na = norm(a);  if na > a_cap, a = a*a_cap/na; end
        X(i, :) = [p, v, a];
        p = p + v*dt + 0.5*a*dt^2;
        v = v + a*dt;
    end
    r = sqrt(X(:,1).^2 + X(:,2).^2);  el = atan2(X(:,2), X(:,1));
    in = t <= 10.0;
    ok = min(r(in)) >= 800 && max(r(in)) <= 2000 && ...
         max(sqrt(X(in,3).^2 + X(in,4).^2)) <= 150 && ...
         max(sqrt(X(in,5).^2 + X(in,6).^2)) <= 20 && ...
         min(el(in)) >= el_lim(1) && max(el(in)) <= el_lim(2);
    if strcmp(family, 'highel'), ok = ok && max(el(in)) >= 70*pi/180; end
    if strcmp(family, 'crosser'), ok = ok && max(abs(diff(el(in)))/dt) >= 0.15; end
    if ok, break, end
end
if ~ok, error('stress_trajectory:sampling', 'No admissible %s trajectory for k = %d.', family, k); end

T.name = sprintf('%s-%02d', family, k);
T.t = t;  T.x = X(:,1);  T.y = X(:,2);  T.vx = X(:,3);  T.vy = X(:,4);  T.ax = X(:,5);  T.ay = X(:,6);
T.blackout = blackout;
T.seed = 10000*fi + k;
end
