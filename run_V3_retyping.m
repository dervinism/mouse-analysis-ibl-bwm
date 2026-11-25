% run_V3_retyping.m - rerun what depends on the functional area types after the V3 retyping of params.m
% (2026-11-08: RT and OT unclassified; LP sensory; VISa, VISrl association; CL, CM, PCN, PF, PIL, PR
% association; motor = MOp, MOs, VAL, VM). Coherence results are per area acronym and are not recomputed.
% Backups: params_before_V3_20261108.m.txt, dependencies\bwmAnalysisResults_before_V3_20261108.mat,
% tables\*_before_V3_20261108.*
cd(fileparts(mfilename('fullpath')));
diary('run_V3_retyping.log'); diary on
qualityFilter = false;
t = tic; fprintf('%s summarise_areas\n', datetime); summarise_areas
fprintf('%s correlate_pupil (%.1f min so far)\n', datetime, toc(t)/60); correlate_pupil
close all; clearvars -except t
qualityFilter = false;
fprintf('%s phase_ordering_stats (%.1f min so far)\n', datetime, toc(t)/60); phase_ordering_stats(false);
fprintf('%s firing_rate_figures (%.1f min so far)\n', datetime, toc(t)/60); firing_rate_figures
close all; clearvars -except t
fprintf('%s manuscript_figures 3/5/8 (%.1f min so far)\n', datetime, toc(t)/60);
manuscript_figures('C:\Users\m329786\Work\Leicester\Paper_1\Figures_v5', ["3" "5" "8"], 'bright', 'unit');
fprintf('%s manuscript_numbers (%.1f min so far)\n', datetime, toc(t)/60); manuscript_numbers(false);
fprintf('%s ALL DONE in %.1f min\n', datetime, toc(t)/60);
diary off
