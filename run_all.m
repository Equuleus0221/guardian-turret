%RUN_ALL  Reproduce every number, table and figure of the paper.
%
%   Run from the repository root:   >> run_all
%   Runtime: about 30 minutes in MATLAB R2025a on a laptop.
%
%   1. design_controller      plant check, compensator, margins   (Fig. 2)
%   2. run_scenarios          scenarios A, B, C                   (Figs. 3-4, Table II)
%   3. run_final_evaluation   200 held-out randomized cases       (Fig. 5, Table III)
%   4. run_stress_tests       independent evaluator, 670 engagements (Section VI)
%
%   Outputs go to ./figures and ./results.

design_controller;
run_scenarios;
run_final_evaluation;
run_stress_tests;
fprintf('\nAll results reproduced. See ./results and ./figures.\n');
