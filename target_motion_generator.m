function [r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
%TARGET_MOTION_GENERATOR  Radar measurement of the target (range, elevation).
%
%   [r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
%
%   The radar delivers a new measurement every Ts = 50 ms (20 Hz):
%       r_meas     = r(t_k)     + n_r,      n_r     ~ N(0, (1 m)^2)
%       theta_meas = theta(t_k) + n_theta,  n_theta ~ N(0, (0.02 deg)^2)
%   meas_valid = 1 only on the call at which a NEW measurement becomes
%   available. Between radar samples, and during measurement-loss intervals,
%   meas_valid = 0 and the outputs hold the last delivered measurement.
%
%   The function must be called with non-decreasing t within one run; a call
%   with t = 0 (or with a different scenario) starts a new run.
%   Noise comes from a portable generator seeded by the scenario, so results
%   are reproducible and identical in MATLAB and Octave.

persistent key k_last r_hold th_hold s
Ts = 0.05;  sigma_r = 1.0;  sigma_th = 0.02*pi/180;

S = scenario_definition(scenario);
newkey = sprintf('%s#%d', S.name, S.seed);
k = floor(t/Ts + 1e-6);
if isempty(key) || ~strcmp(key, newkey) || t <= 0 || k < k_last
    key = newkey;
    k_last = -1;
    r_hold = 0; th_hold = 0;
    s = 48271 + 104729*S.seed;          % noise generator state
    for j = 1:10, s = mod(16807*s, 2147483647); end
end

meas_valid = 0;
if k > k_last
    k_last = k;
    tk = k*Ts;
    [x, y, ~, ~, ~, ~, lost] = target_truth(tk, scenario);
    % draw the noise even while the radar is blind, so that the noise
    % sequence after a measurement loss does not depend on its duration
    [n, s] = gauss2(s);
    if ~lost
        r_hold  = sqrt(x^2 + y^2) + sigma_r*n(1);
        th_hold = atan2(y, x) + sigma_th*n(2);
        meas_valid = 1;
    end
end
r_meas = r_hold;
theta_meas = th_hold;
end

function [n, s] = gauss2(s)
% two independent N(0,1) samples (Box-Muller on a Park-Miller generator)
s = mod(16807*s, 2147483647);  u1 = s/2147483647;
s = mod(16807*s, 2147483647);  u2 = s/2147483647;
rho = sqrt(-2*log(u1));
n = [rho*cos(2*pi*u2); rho*sin(2*pi*u2)];
end
