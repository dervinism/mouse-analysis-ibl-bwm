% run_medians_fig5.m - Figure 5 with medians / Wilcoxon / Brown-Forsythe and sourced UoL/Allen data (2026-11-08),
% after run_rhoGE0.m (IBL rho >= 0). Outputs: figures_qc0?\firing_rate_distros, Paper_1\Figures_v6\Figure5, tables.
cd(fileparts(mfilename('fullpath')));
diary('run_medians_fig5.log'); diary on
t = tic;
qualityFilter = false; useMedians = true;
fprintf('%s firing_rate_figures\n', datetime); firing_rate_figures
close all; clearvars -except t
fprintf('%s manuscript_figures 5 (%.1f min so far)\n', datetime, toc(t)/60);
manuscript_figures('C:\Users\m329786\Work\Leicester\Paper_1\Figures_v6', "5", 'bright', 'unit');
fprintf('%s manuscript_numbers (%.1f min so far)\n', datetime, toc(t)/60); manuscript_numbers(false);
fprintf('%s ALL DONE in %.1f min\n', datetime, toc(t)/60);
diary off
