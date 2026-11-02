function run_figures_probe1_fix(qualityFilter)
% run_figures_probe1_fix(qualityFilter)
%
% Regenerates the IBL figures and statistics after the probe01 fix
% (2026-10-30) for one unit selection: qualityFilter = false (all Kilosort
% clusters, figures in ./figures) or true (units passing IBL quality
% control, figures in ./figures_qc). Loads the preprocessed data and the
% analysis results once and runs firing_rate_figures,
% phase_distribution_figures, phase_profiles_pupil and
% phase_coherence_correlations, logging the printed statistics and the
% duration of each script to run_figures_probe1_fix_qc<0|1>.log.

cd(fileparts(mfilename('fullpath')));
diary(sprintf('run_figures_probe1_fix_qc%d.log', qualityFilter)); diary on
tAll = tic;
fprintf('%s qualityFilter = %d\n', datetime, qualityFilter);
params
load(fullfile(processedDataFolder, 'bwmPreprocessedData2.mat')); %#ok<*LOAD>
load(fullfile(processedDataFolder, ['bwmAnalysisResults' resultsTag '.mat']));
for script = {'firing_rate_figures', 'phase_distribution_figures', ...
    'phase_profiles_pupil', 'phase_coherence_correlations'}
  t = tic;
  try
    eval(script{1});
    fprintf('%s %s done in %.1f min\n', datetime, script{1}, toc(t)/60);
  catch err
    fprintf(2, '%s %s FAILED: %s\n', datetime, script{1}, getReport(err, 'extended'));
  end
  close all
end
fprintf('%s all done in %.1f min\n', datetime, toc(tAll)/60);
diary off
end
