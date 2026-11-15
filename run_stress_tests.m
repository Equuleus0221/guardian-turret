%RUN_STRESS_TESTS  Stress tests of the finalized design with an independent evaluator.
%
%   Part 1: seven target-motion families that differ from the sinusoidal
%           maneuvers of scenario_definition.m (see stress_trajectory.m).
%   Part 2: the same controller driven through different calling conventions
%           (step size, integrator, meas_valid convention, initial state),
%           as an external evaluation script might use.
%
%   Nothing is retuned. Output: results/stress_metrics.txt, results/stress.mat

clear; close all;
here = fileparts(mfilename('fullpath'));  addpath(here);
resdir = fullfile(here, 'results');  if ~exist(resdir, 'dir'), mkdir(resdir); end
turret_controller('override', []);
N_FAM = 40;  N_IF = 30;
fid = fopen(fullfile(resdir, 'stress_metrics.txt'), 'w');
out = @(varargin) tee(fid, varargin{:});
hdr = '%-34s %9s %11s %7s %7s %7s %7s %7s %7s\n';
row = '%-34s %4d/%-4d %5d/%-5d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n';

%% Part 1: target-motion families
families = {'straight', 'jink', 'turn', 'reversal', 'crosser', 'highel', 'multiloss'};
out('PART 1: target-motion families (%d cases each, default interface)\n', N_FAM);
out(hdr, 'family', 'intercept', 'pulses hit', 't_int', 'max', 'miss', 'max', 'trkRMS', 'trkMax');
ALL = struct('label', {}, 'R', {});
tot = zeros(1, 4);
for f = 1:numel(families)
    clear Rs
    for k = 1:N_FAM
        T = stress_trajectory(families{f}, k);
        Rs(k) = stress_engagement(T); %#ok<SAGROW>
    end
    ALL(end+1) = struct('label', families{f}, 'R', Rs); %#ok<SAGROW>
    tot = tot + report(out, row, families{f}, Rs);
end
out('%-34s %4d/%-4d %5d/%-5d\n', 'TOTAL', tot(1), tot(2), tot(3), tot(4));

%% Part 2: calling conventions (jink family)
variants = { ...
  'baseline: rk4, 1 ms, pulse',        struct(); ...
  'euler, 1 ms',                       struct('integrator', 'euler'); ...
  'rk4, 2 ms',                         struct('dt', 2e-3); ...
  'rk4, 5 ms',                         struct('dt', 5e-3); ...
  'rk4, 10 ms',                        struct('dt', 1e-2); ...
  'controller inside RK4 stages',      struct('integrator', 'rk4stage'); ...
  'meas_valid held high',              struct('meas', 'held'); ...
  'held + controller in RK4 stages',   struct('meas', 'held', 'integrator', 'rk4stage'); ...
  'NaN measurements during loss',      struct('meas', 'nan'); ...
  'turret starts at +85 deg',          struct('theta0', 85*pi/180); ...
  'turret starts at -10 deg',          struct('theta0', -10*pi/180); ...
  'radar phase 20 ms',                 struct('phase', 0.02); ...
  'first call at t = 1 ms',            struct('t_first', 1e-3)};
out('\nPART 2: calling conventions (%d jink cases each)\n', N_IF);
out(hdr, 'variant', 'intercept', 'pulses hit', 't_int', 'max', 'miss', 'max', 'trkRMS', 'trkMax');
tot2 = zeros(1, 4);
for v = 1:size(variants, 1)
    clear Rs
    for k = 1:N_IF
        T = stress_trajectory('jink', k);
        Rs(k) = stress_engagement(T, variants{v, 2}); %#ok<SAGROW>
    end
    ALL(end+1) = struct('label', variants{v, 1}, 'R', Rs); %#ok<SAGROW>
    tot2 = tot2 + report(out, row, variants{v, 1}, Rs);
end
out('%-34s %4d/%-4d %5d/%-5d\n', 'TOTAL', tot2(1), tot2(2), tot2(3), tot2(4));

%% misses, if any
out('\nPulses outside R_eng = 2 m:\n');
nmiss = 0;
for a = 1:numel(ALL)
    for k = 1:numel(ALL(a).R)
        R = ALL(a).R(k);
        for j = find(R.pulse_miss > 2)'
            out('  %-34s %-14s t = %.3f s, miss = %.2f m\n', ALL(a).label, R.name, R.pulse_t(j), R.pulse_miss(j));
            nmiss = nmiss + 1;
        end
        if ~R.finite, out('  %-34s %-14s NON-FINITE OUTPUT\n', ALL(a).label, R.name); end
    end
end
if nmiss == 0, out('  none\n'); end
fclose(fid);
save(fullfile(resdir, 'stress.mat'), 'ALL');

function c = report(out, row, label, Rs)
ti = [Rs.t_intercept];  pm = vertcat(Rs.pulse_miss);
c = [sum([Rs.success]), numel(Rs), sum([Rs.n_hits]), sum([Rs.n_pulses])];
out(row, label, c(1), c(2), c(3), c(4), mean(ti, 'omitnan'), max(ti), mean(pm), max([pm; NaN]), ...
    sqrt(mean([Rs.miss_rms].^2)), max([Rs.miss_max]));
end
