function A = loop_analysis(C, P, scale_wn)
%LOOP_ANALYSIS  Linear analysis of the turret elevation loop (toolbox-free).
%
%   A = loop_analysis()            uses the design in turret_controller.m
%   A = loop_analysis(C, P)        uses the given controller/plant structs
%   A = loop_analysis(C, P, k)     plant resonance scaled by k (robustness)
%
%   Returns loop polynomials, closed-loop poles, gain/phase margins,
%   crossover, bandwidth, sensitivity peaks, error constants and the unit
%   step response of T(s) = L/(1+L).

if nargin < 1 || isempty(C), C = turret_controller('params'); end
if nargin < 2 || isempty(P), P = turret_plant('params'); end
if nargin < 3, scale_wn = 1; end

% plant Gp = num/den, optionally with a shifted resonance (inductance error)
den = P.den;
if scale_wn ~= 1
    Lm = P.Lm/scale_wn^2;                      % wn ~ 1/sqrt(L)
    den = [P.Jeq*Lm, P.Jeq*P.Rm + P.b*Lm, P.b*P.Rm + P.N^2*P.Kt*P.Ke, 0];
end
A.Gp_num = P.num/den(1);  A.Gp_den = den/den(1);
A.L_num = conv(C.Gc_num, A.Gp_num);
A.L_den = conv(C.Gc_den, A.Gp_den);
A.cl_den = polyadd(A.L_den, A.L_num);
A.cl_poles = roots(A.cl_den);
A.stable = all(real(A.cl_poles) < 0);

w = logspace(-1, 3, 20000)';
s = 1i*w;
Lw = polyval(A.L_num, s)./polyval(A.L_den, s);
A.w = w;  A.L = Lw;
A.mag_dB = 20*log10(abs(Lw));
A.phase_deg = unwrap(angle(Lw))*180/pi;
% the loop has two integrators: phase starts at -180 deg
A.phase_deg = A.phase_deg - 360*round((A.phase_deg(1) + 180)/360);
A.S = 1./(1 + Lw);  A.T = Lw./(1 + Lw);

% gain crossover(s) and phase margin (worst case)
ix = find(diff(sign(A.mag_dB)) ~= 0);
wc = zeros(numel(ix), 1);  pm = wc;
for k = 1:numel(ix)
    j = ix(k);
    f = A.mag_dB(j)/(A.mag_dB(j) - A.mag_dB(j+1));
    wc(k) = w(j)*(w(j+1)/w(j))^f;
    pm(k) = 180 + A.phase_deg(j) + f*(A.phase_deg(j+1) - A.phase_deg(j));
end
A.wc_all = wc;  A.pm_all = pm;
[A.PM, kk] = min(pm);  A.wc = wc(kk);
% phase crossover(s) and gain margin
ph = A.phase_deg + 180;
ix = find(diff(sign(ph)) ~= 0);
gm = inf;  wpc = NaN;
for k = 1:numel(ix)
    j = ix(k);
    f = ph(j)/(ph(j) - ph(j+1));
    g = -(A.mag_dB(j) + f*(A.mag_dB(j+1) - A.mag_dB(j)));
    if abs(g) < abs(gm), gm = g;  wpc = w(j)*(w(j+1)/w(j))^f; end
end
A.GM_dB = gm;  A.w_pc = wpc;
A.Ms = max(abs(A.S));  A.Mt = max(abs(A.T));
A.delay_margin = (A.PM*pi/180)/A.wc;
jb = find(abs(A.T) < 1/sqrt(2), 1);
A.w_bw = w(jb);
A.Ka = C.Ki*P.num/den(3);                      % lim s^2 L(s)
% noise-equivalent bandwidth of T(s) [rad/s] (one-sided integral of |T|^2)
wl = linspace(0, 2000, 400001)';
Tl = polyval(A.L_num, 1i*wl)./polyval(A.cl_den, 1i*wl);
A.w_neb = trapz(wl, abs(Tl).^2);

% unit step response of T(s) via its state-space form
nd = numel(A.cl_den) - 1;
a = A.cl_den/A.cl_den(1);
bnum = [zeros(1, nd + 1 - numel(A.L_num)), A.L_num]/A.cl_den(1);
Ac = [zeros(nd-1, 1), eye(nd-1); -fliplr(a(2:end))];
Bc = [zeros(nd-1, 1); 1];
Cc = fliplr(bnum(2:end)) - bnum(1)*fliplr(a(2:end));
h = 1e-3;  tt = (0:h:3)';
Md = expm([Ac, Bc; zeros(1, nd+1)]*h);
Ad = Md(1:nd, 1:nd);  Bd = Md(1:nd, end);
xs = zeros(nd, 1);  ys = zeros(size(tt));
for k = 1:numel(tt)
    ys(k) = Cc*xs + bnum(1);
    xs = Ad*xs + Bd;
end
A.step_t = tt;  A.step_y = ys;
A.overshoot = 100*(max(ys) - 1);
A.t_rise = tt(find(ys >= 0.9, 1)) - tt(find(ys >= 0.1, 1));
out2 = find(abs(ys - 1) > 0.02, 1, 'last');
A.t_settle2 = tt(min(out2 + 1, numel(tt)));
end

function p = polyadd(a, b)
n = max(numel(a), numel(b));
p = [zeros(1, n - numel(a)), a] + [zeros(1, n - numel(b)), b];
end
