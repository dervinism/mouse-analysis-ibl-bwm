function T = phase_ordering_stats(qualityFilter, freqOI)
% T = phase_ordering_stats(qualityFilter, freqOI)
%
% Table 7 of the manuscript: phase ordering of IBL units at freqOI (default
% 0.3 Hz) referenced to the pupil signal, for the joint population and the
% positive and negative subpopulations. Unit-level phases are pooled across
% recordings within each area group (Hp, nCx, Th, sensory, association,
% motor, association-motor) and two groups are compared as described in
% the Methods: the joint population with the Mardia-Watson-Wheeler test
% (MWW, W referred to chi^2 with 2 df); the subpopulations with the
% Watson-Williams test (WW, F with 1 and n-2 df), unless the
% sample-size-weighted mean resultant length rw is below 0.45, in which
% case MWW is used instead (marked with a dagger).
%
% Args:
%   qualityFilter (logical, optional, positional): false (default) for
%     all Kilosort clusters, true for units passing IBL quality control;
%     selects the analysis results file via params.m.
%   freqOI (numeric, optional, positional): frequency in Hz (default 0.3).
%
% Returns:
%   T (table): one row per comparison and population; also written to
%     tables/phase_ordering_stats_qc<0|1>.csv.

if nargin < 1 || isempty(qualityFilter)
  qualityFilter = false; %#ok<NASGU> read by params
end
if nargin < 2 || isempty(freqOI)
  freqOI = 0.3;
end
cd(fileparts(mfilename('fullpath')));
params
load(fullfile(processedDataFolder, ['bwmAnalysisResults' resultsTag '.mat']), 'infraslowAnalyses');
areaTable = infraslowAnalyses.areaSummaries.areaTable;

% pooled version: association-motor vs sensory (+ association vs motor, which decides the pooling);
% non-pooled version: association vs sensory and motor vs sensory as well
comparisons = {'Hp', 'nCx'; 'nCx', 'Th'; 'association-motor', 'sensory'; 'association', 'motor'; ...
  'association', 'sensory'; 'motor', 'sensory'};
populations = {'all', 'spikingPupilCoh'; 'positive', 'spikingPupilCohPositive'; ...
  'negative', 'spikingPupilCohNegative'};
rwThreshold = 0.45;

rows = {};
for iPop = 1:size(populations,1)
  cohStruct = infraslowAnalyses.(populations{iPop,2});
  for iComp = 1:size(comparisons,1)
    phase = cell(1,2);
    for iGroup = 1:2
      phase{iGroup} = local_group_phases(cohStruct, areaTable, comparisons{iComp,iGroup}, freqOI);
    end
    n = cellfun(@numel, phase);
    if any(n < 2)
      rows(end+1,:) = {comparisons{iComp,1}, comparisons{iComp,2}, populations{iPop,1}, n(1), n(2), ...
        NaN, NaN, NaN, NaN, '', NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
      continue
    end
    mu = cellfun(@circ_mean, phase);
    med = cellfun(@local_circ_median, phase);
    r = cellfun(@circ_r, phase);
    rw = sum(n .* r) / sum(n);
    if strcmpi(populations{iPop,1}, 'all') || rw < rwThreshold
      [~, p, W] = mardiatestn_circ_equal(phase);
      test = 'MWW';
      if ~strcmpi(populations{iPop,1}, 'all')
        test = 'MWW (rw < 0.45)';
      end
      stat = W; df1 = 2; df2 = NaN;
    else
      [p, wwTable] = circ_wwtest(phase{1}, phase{2});
      test = 'WW';
      stat = wwTable{2,5}; df1 = wwTable{2,2}; df2 = wwTable{3,2};
    end
    rows(end+1,:) = {comparisons{iComp,1}, comparisons{iComp,2}, populations{iPop,1}, n(1), n(2), ...
      mu(1), mu(2), med(1), med(2), test, stat, df1, df2, p, rw}; %#ok<AGROW>
  end
end
T = cell2table(rows, 'VariableNames', {'group1', 'group2', 'population', 'n1', 'n2', ...
  'mean1', 'mean2', 'median1', 'median2', 'test', 'statistic', 'df1', 'df2', 'p', 'rw'});
if ~isfolder('tables')
  mkdir('tables');
end
writetable(T, fullfile('tables', sprintf('phase_ordering_stats_qc%d.csv', qualityFilter)));
disp(T);
end


function phase = local_group_phases(cohStruct, areaTable, group, freqOI)
% Unit phases at freqOI of all units in the areas making up the group.
areaInds = find(getAreaInds(group, areaTable));
phase = [];
for iArea = 1:numel(areaInds)
  areaData = getAreaCohData(cohStruct, areaTable, areaInds(iArea));
  for iRec = 1:numel(areaData)
    if isstruct(areaData{iRec}) && isfield(areaData{iRec}, 'fullInterpCoherence')
      recData = areaData{iRec}.fullInterpCoherence;
      [dist, freqInd] = min(abs(recData.frequency(1,:) - freqOI));
      if dist < 1e-6
        phase = [phase; recData.phase(:,freqInd)]; %#ok<AGROW>
      end
    end
  end
end
phase = phase(~isnan(phase));
end


function med = local_circ_median(alpha)
% Circular median: the angle minimising the mean circular distance to the
% sample, searched on a 0.001 rad grid (circ_median is O(n^2) in memory).
grid = (-pi:0.001:pi)';
cost = zeros(size(grid));
for i = 1:numel(grid)
  cost(i) = mean(abs(angle(exp(1i*(alpha - grid(i))))));
end
[~, iMin] = min(cost);
med = grid(iMin);
end
