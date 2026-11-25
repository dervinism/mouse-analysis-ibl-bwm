% run_rhoGE0.m - IBL positive subpopulation switched from rho > 0 to rho >= 0, as in the UoL/Allen code
% (2026-11-08, manuscript comment on the Unit classification paragraph). correlate_pupil.m and
% coherence_pupil.m were edited; spikingPupilCohPositive is rebuilt from spikingPupilCoh (no coherence rerun).
% Backups: dependencies\bwmAnalysisResults_before_rhoGE0_20261108.mat, tables\*_before_rhoGE0_20261108.*,
% correlate_pupil_before_rhoGE0_20261108.m.txt, coherence_pupil_before_rhoGE0_20261108.m.txt
cd(fileparts(mfilename('fullpath')));
diary('run_rhoGE0.log'); diary on
t = tic;
qualityFilter = false;
params
f = fullfile(processedDataFolder, ['bwmAnalysisResults' resultsTag '.mat']);
b = fullfile(processedDataFolder, ['bwmAnalysisResults' resultsTag '_before_rhoGE0_20261108.mat']);
if ~exist(b, 'file'), copyfile(f, b); fprintf('%s backup %s\n', datetime, b); end
tb = dir(fullfile('tables', '*_qc0*.*'));
for k = 1:numel(tb)
  if contains(tb(k).name, '_before_'), continue; end
  [~, nm, ext] = fileparts(tb(k).name);
  dst = fullfile('tables', [nm '_before_rhoGE0_20261108' ext]);
  if ~exist(dst, 'file'), copyfile(fullfile('tables', tb(k).name), dst); end
end
fprintf('%s rebuild positive population\n', datetime);
load(f, 'infraslowAnalyses');
infraslowAnalyses = rebuild_positive_rhoGE0(infraslowAnalyses);
saveAnalysisResults(f, infraslowAnalyses);
clearvars -except t
qualityFilter = false;
fprintf('%s correlate_pupil (%.1f min so far)\n', datetime, toc(t)/60); correlate_pupil
close all; clearvars -except t
qualityFilter = false;
fprintf('%s phase_ordering_stats (%.1f min so far)\n', datetime, toc(t)/60); phase_ordering_stats(false);
fprintf('%s firing_rate_figures (%.1f min so far)\n', datetime, toc(t)/60); firing_rate_figures
close all; clearvars -except t
fprintf('%s manuscript_figures 3/5/8 (%.1f min so far)\n', datetime, toc(t)/60);
manuscript_figures('C:\Users\m329786\Work\Leicester\Paper_1\Figures_v6', ["3" "5" "8"], 'bright', 'unit');
fprintf('%s manuscript_numbers (%.1f min so far)\n', datetime, toc(t)/60); manuscript_numbers(false);
fprintf('%s ALL DONE in %.1f min\n', datetime, toc(t)/60);
diary off
