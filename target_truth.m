function [x, y, vx, vy, ax, ay, lost] = target_truth(t, scenario)
%TARGET_TRUTH  True target kinematics in the elevation plane.
%
%   [x, y, vx, vy, ax, ay, lost] = target_truth(t, scenario)
%
%   t may be a vector. x is horizontal distance from the turret [m], y is
%   height above the trunnion [m]. 'lost' is true while the radar
%   measurement is unavailable. See SCENARIO_DEFINITION for the scenarios.

S = scenario_definition(scenario);
tau = max(t - S.t_m, 0);
[dx, dvx, ax] = axis_maneuver(S.Ax, S.wx, tau);
[dy, dvy, ay] = axis_maneuver(S.Ay, S.wy, tau);
x = S.x0 + S.vx0*t + dx;   vx = S.vx0 + dvx;
y = S.y0 + S.vy0*t + dy;   vy = S.vy0 + dvy;

lost = false(size(t));
for k = 1:size(S.blackout, 1)
    t0 = S.blackout(k, 1);  t1 = t0 + S.blackout(k, 2);
    lost = lost | (t >= t0 - 1e-9 & t < t1 - 1e-9);
end
end

function [p, v, a] = axis_maneuver(A, w, tau)
p = zeros(size(tau)); v = p; a = p;
for j = 1:numel(A)
    a = a + A(j)*sin(w(j)*tau);
    v = v + A(j)/w(j)*(1 - cos(w(j)*tau));
    p = p + A(j)/w(j)*(tau - sin(w(j)*tau)/w(j));
end
end
