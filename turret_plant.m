function [theta_dot_turret, theta_ddot_turret, i_dot] = turret_plant(t, V_m, theta_turret, theta_dot_turret, i)
%TURRET_PLANT  Elevation-drive dynamics of the Guardian laser turret.
%
%   [theta_dot, theta_ddot, i_dot] = turret_plant(t, V_m, theta, theta_dot, i)
%
%   Inputs
%     t            time [s] (unused, the plant is time-invariant)
%     V_m          motor armature voltage [V]
%     theta_turret turret elevation [rad]
%     theta_dot_turret  turret angular velocity [rad/s]
%     i            motor armature current [A]
%   Outputs (state derivatives)
%     theta_dot_turret  d(theta)/dt   [rad/s]
%     theta_ddot_turret d2(theta)/dt2 [rad/s^2]
%     i_dot             di/dt         [A/s]
%
%   Model (all quantities referred to the turret/trunnion side, N = 100):
%     motor speed      w_m = N*theta_dot
%     electrical       L di/dt = V_m - R i - Ke*N*theta_dot
%     mechanical       (J_turret + N^2 J_m) theta_ddot = N*Kt*i - b*theta_dot
%   Assumptions: rigid, lossless, backlash-free gear train; trunnion-balanced
%   cradle (no gravity torque about the elevation axis); viscous damping only.
%
%   Transfer function:  theta(s)/V(s) = N Kt / ( s [ (Jeq s + b)(L s + R) + N^2 Kt Ke ] )
%
%   P = turret_plant('params') returns the parameter struct.
%
%   This file is self-contained (no dependencies).

persistent P
if isempty(P)
    P = plant_params();
end
if nargin == 1 && ischar(t)
    theta_dot_turret = P;
    return
end

w = theta_dot_turret;

% Mechanical and electrical dynamics
theta_ddot_turret = (P.N*P.Kt*i - P.b*w) / P.Jeq;
i_dot             = (V_m - P.Rm*i - P.N*P.Ke*w) / P.Lm;

% Mechanical elevation stops (-10 deg ... +85 deg). An inelastic stop:
% outward motion is arrested, inward motion is unaffected.
at_hi = (theta_turret >= P.theta_max) && (w > 0 || (w == 0 && theta_ddot_turret > 0));
at_lo = (theta_turret <= P.theta_min) && (w < 0 || (w == 0 && theta_ddot_turret < 0));
if at_hi || at_lo
    theta_dot_turret  = 0;
    theta_ddot_turret = -w / P.tau_stop;   % velocity is driven to zero within ~tau_stop
end
end

% -------------------------------------------------------------------------
function P = plant_params()
P.J_turret = 1200;      % turret + cradle inertia about trunnion [kg m^2]
P.J_m      = 0.20;      % motor rotor inertia [kg m^2]
P.N        = 100;       % gear reduction (motor speed / turret speed)
P.b        = 1500;      % turret-side viscous damping [N m s/rad]
P.Kt       = 10;        % torque constant [N m/A]
P.Ke       = 10;        % back-emf constant [V s/rad]
P.Rm       = 1.0;       % armature resistance [ohm]
P.Lm       = 0.10;      % armature inductance [H]
P.theta_min = -10*pi/180;   % lower elevation stop [rad]
P.theta_max =  85*pi/180;   % upper elevation stop [rad]
P.tau_stop  = 5e-3;         % stop arrest time constant [s]

% Derived quantities
P.Jeq = P.J_turret + P.N^2*P.J_m;                 % 3200 kg m^2
% Gp(s) = num / (den(1) s^3 + den(2) s^2 + den(3) s)
P.num = P.N*P.Kt;
P.den = [P.Jeq*P.Lm, P.Jeq*P.Rm + P.b*P.Lm, P.b*P.Rm + P.N^2*P.Kt*P.Ke, 0];
P.wn_p   = sqrt(P.den(3)/P.den(1));               % resonance [rad/s]
P.zeta_p = P.den(2)/(2*sqrt(P.den(1)*P.den(3)));  % damping ratio
P.k_v    = P.num/P.den(3);                        % velocity gain [rad/s per V]
end
