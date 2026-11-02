function run_one_figure_script(qualityFilter, scriptName)
% run_one_figure_script(qualityFilter, scriptName)
%
% Runs one IBL figure script (e.g. 'firing_rate_figures') for one unit
% selection (qualityFilter false/true), loading the preprocessed data and
% the matching analysis results first, and logs its output to
% run_<scriptName>_qc<0|1>.log.

cd(fileparts(mfilename('fullpath')));
diary(sprintf('run_%s_qc%d.log', scriptName, qualityFilter)); diary on
t = tic;
fprintf('%s %s, qualityFilter = %d\n', datetime, scriptName, qualityFilter);
params
load(fullfile(processedDataFolder, 'bwmPreprocessedData2.mat')); %#ok<*LOAD>
load(fullfile(processedDataFolder, ['bwmAnalysisResults' resultsTag '.mat']));
eval(scriptName);
close all
fprintf('%s %s done in %.1f min\n', datetime, scriptName, toc(t)/60);
diary off
end
