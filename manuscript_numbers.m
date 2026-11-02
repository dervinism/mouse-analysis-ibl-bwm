function manuscript_numbers(qualityFilter)
% manuscript_numbers(qualityFilter)
%
% IBL numbers quoted in the manuscript text after the probe01 fix, from the
% saved analysis results (same computations as the figure scripts, which do
% not print them):
%   - Table 4: per-recording positive fraction per area (mean, SD, 95% CI
%     half-width from the t distribution, recording count);
%   - Results t-tests (two-sample, pooled variance, recordings as units):
%     Th vs nCx, nCx vs Hp, sensory vs motor + association (recordings of
%     the two concatenated), motor vs association;
%   - Fig. 5C: DG log10 firing rate in the most constricted vs most dilated
%     pupil (bottom/top 12.5%), meanTest 'ANOVARM' and varTest as in
%     firing_rate_figures.m;
%   - Fig. 5E: Spearman correlation of log firing rate change with positive
%     fraction across the non-overlapping area groups (incMask and the
%     hard-coded UoL/Allen values of firing_rate_figures.m);
%   - Fig. 8D-F: circular-linear correlation of area phase at 0.3 Hz with
%     positive fraction for all/positive/negative units (incMask of
%     phase_profiles_pupil.m);
%   - dataset counts (sessions, insertions, mice, units).
% Writes tables/table4_qc<0|1>.csv and tables/manuscript_numbers_qc<0|1>.txt.

if nargin < 1 || isempty(qualityFilter)
  qualityFilter = false; %#ok<NASGU> read by params
end
cd(fileparts(mfilename('fullpath')));
params
load(fullfile(processedDataFolder, ['bwmAnalysisResults' resultsTag '.mat']), 'infraslowAnalyses');
areaTable = infraslowAnalyses.areaSummaries.areaTable;
corrOI = infraslowAnalyses.spikingPupilCorr.areasOI;
outFile = fullfile('tables', sprintf('manuscript_numbers_qc%d.txt', qualityFilter));
fid = fopen(outFile, 'w');

%% Table 4
rows = {'VPL-VPM','VPL-VPM'; 'Sensory thalamus','sensory-Th'; 'Brainstem','BS'; ...
  'Association thalamus','association-Th'; 'Thalamus','Th'; 'Epithalamus','EPT'; ...
  'Cerebellum','Cereb'; 'Motor thalamus','motor-Th'; 'Subthalamus','STh'; ...
  'Sensory areas','sensory'; 'V2','V2'; 'Midbrain','MBr'; 'Visual cortex','VIS'; ...
  'Sensory neocortex','sensory-nCx'; 'Basal ganglia','BG'; 'LG','LG'; 'SSp','SSp'; ...
  'Paleocortex','pCx'; 'Claustrum','CLA'; 'RSC','RSP'; 'Neocortex','nCx'; ...
  'Motor areas','motor'; 'Periallocortex','paCx'; 'Proisocortex','piCx'; ...
  'Hypothalamus','Hyp'; 'Association areas','association'; 'V1','VISp'; ...
  'Amygdalar areas','AA'; 'Basal forebrain','BF'; 'Association neocortex','association-nCx'; ...
  'Motor neocortex','motor-nCx'; 'CA','CA'; 'Hippocampus','Hp'; 'DG','DG'};
perRec = @(acr) local_per_rec(corrOI, areasOI, acr);
T4 = cell(size(rows,1), 6);
for i = 1:size(rows,1)
  x = perRec(rows{i,2}); n = numel(x);
  T4(i,:) = {rows{i,1}, rows{i,2}, mean(x), std(x), tinv(0.975, n-1)*std(x)/sqrt(n), n};
end
T4 = cell2table(T4, 'VariableNames', {'area', 'acronym', 'mean', 'sd', 'ci95', 'recordings'});
T4 = sortrows(T4, 'mean', 'descend');
writetable(T4, fullfile('tables', sprintf('table4_qc%d.csv', qualityFilter)));
say(fid, 'Table 4 (per-recording positive fraction)\n');
for i = 1:height(T4)
  say(fid, '  %-22s %.3f  SD %.3f  CI %.3f  n %d   (stored mean %.3f)\n', T4.area{i}, T4.mean(i), T4.sd(i), ...
    T4.ci95(i), T4.recordings(i), corrOI.positiveSpearmanFractionsMeans(strcmp(areasOI, T4.acronym{i}), 1));
end

%% Results t-tests (paragraph on the IBL gradient)
tests = {'Th', {'Th'}, 'nCx', {'nCx'}; 'nCx', {'nCx'}, 'Hp', {'Hp'}; ...
  'sensory', {'sensory'}, 'motor + association', {'motor', 'association'}; ...
  'motor', {'motor'}, 'association', {'association'}};
say(fid, '\nt-tests (recordings as units, pooled variance)\n');
for i = 1:size(tests,1)
  a = cell2mat(cellfun(perRec, tests{i,2}, 'UniformOutput', false)');
  b = cell2mat(cellfun(perRec, tests{i,4}, 'UniformOutput', false)');
  [~, p, ~, st] = ttest2(a, b);
  say(fid, '  %s vs %s: %.3f±%.3f (n %d) vs %.3f±%.3f (n %d), t(%d) = %.2f, p = %.2g\n', tests{i,1}, tests{i,3}, ...
    mean(a), tinv(0.975, numel(a)-1)*std(a)/sqrt(numel(a)), numel(a), mean(b), ...
    tinv(0.975, numel(b)-1)*std(b)/sqrt(numel(b)), numel(b), st.df, st.tstat, p);
end

%% Fig. 8D-F: phase vs positive fraction (phase_profiles_pupil.m)
incMask8 = [false false false false false true true false true true true ...
  true true true false false false false false false false false false ...
  false false false false false false false false false false false false ...
  true true true true true true true true true true true true false false ...
  false true false];
pops = {'all', 'spikingPupilCoh'; 'positive', 'spikingPupilCohPositive'; 'negative', 'spikingPupilCohNegative'};
say(fid, '\nFig. 8D-F: circular-linear correlation, area phase at 0.3 Hz vs positive fraction (%d area groups: %s)\n', ...
  sum(incMask8), strjoin(areasOI(incMask8), ', '));
for iPop = 1:3
  cohStruct = infraslowAnalyses.(pops{iPop,2});
  ph = nan(numel(areasOI), 1);
  for iA = 1:numel(areasOI)
    if ~incMask8(iA), continue, end
    inds = find(getAreaInds(areasOI{iA}, areaTable)); phase = [];
    for k = 1:numel(inds)
      d = getAreaCohData(cohStruct, areaTable, inds(k));
      for r = 1:numel(d)
        if isstruct(d{r}) && isfield(d{r}, 'fullInterpCoherence')
          fi = abs(d{r}.fullInterpCoherence.frequency(1,:) - 0.3) < 1e-6;
          phase = [phase; d{r}.fullInterpCoherence.phase(:, fi)]; %#ok<AGROW>
        end
      end
    end
    if ~isempty(phase), ph(iA) = datamean(phase, 'circularNP'); end
  end
  [rho, pv] = corrMulti(ph(incMask8)', corrOI.positiveSpearmanFractionsMeans(incMask8, 1)', 'circlinearnp');
  say(fid, '  %s units: rho = %.2f, p = %.3g\n', pops{iPop,1}, rho, pv);
end

%% Fig. 5C and 5E: firing rates per pupil state (firing_rate_figures.m)
load(fullfile(processedDataFolder, 'bwmPreprocessedData2.mat'), 'infraslowData');
firingRateThr = 50/3600;
nAreas = size(areaLabels,1);
firingRates = cell(nAreas,1);
for iArea = 1:nAreas
  areaInds = infraslowAnalyses.areaSummaries.groupedUnitInds{iArea};
  fr = nan(size(areaInds,1), 2);
  for iExp = unique(areaInds(:,1))'
    unitMask = areaInds(:,1) == iExp;
    expData = infraslowData.experimentData{iExp};
    if ~isempty(expData.leftPupilAreaSize) && ~isempty(expData.rightPupilAreaSize)
      pupil = mean([expData.leftPupilAreaSize; expData.rightPupilAreaSize]);
    elseif ~isempty(expData.leftPupilAreaSize), pupil = expData.leftPupilAreaSize;
    elseif ~isempty(expData.rightPupilAreaSize), pupil = expData.rightPupilAreaSize;
    else, continue
    end
    counts = expData.spikeCounts(areaInds(unitMask,2),:);
    dur = size(counts,2) / expData.samplingRate;
    rate = sum(counts,2) / dur; active = rate > firingRateThr;
    pv = prctile(pupil, [12.5 87.5]);
    lo = zeros(sum(unitMask),1); hi = lo;
    lo(active) = sum(counts(active, pupil <= pv(1)), 2) / dur;   % whole-session duration, as in the figure script
    hi(active) = sum(counts(active, pupil > pv(2)), 2) / dur;
    lo(~lo) = nan; hi(~hi) = nan;
    fr(unitMask,:) = [log10(lo) log10(hi)];
  end
  firingRates{iArea} = fr;
end
frOI = @(acr) concatenateCells(firingRates(getAreaInds(acr, areaTable)));
logRates = frOI('DG');
sM = meanTest(logRates, 'ANOVARM'); sV = varTest(logRates);
paired = all(isfinite(logRates), 2);
[~, pT, ~, stT] = ttest(logRates(paired,2), logRates(paired,1));
say(fid, '\nFig. 5C (IBL DG): units %d, with both states %d\n', size(logRates,1), sum(paired));
say(fid, '  mean constricted %.3f, dilated %.3f (omitnan); paired t(%d) = %.1f, p = %.3g\n', ...
  mean(logRates(:,1), 'omitnan'), mean(logRates(:,2), 'omitnan'), stT.df, stT.tstat, pT);
say(fid, '  SD constricted %.2f, dilated %.2f\n', std(logRates(:,1), 'omitnan'), std(logRates(:,2), 'omitnan'));
say(fid, '  meanTest ANOVARM: %s\n', strtrim(evalc('disp(sM)')));
say(fid, '  varTest (Levene quadratic): %s\n', strtrim(evalc('disp(sV)')));

incMask5 = [false false false false false true true false true true true ...
  true true true false false false false false false false false false ...
  false false false false false false false false false false false false ...
  true true true true true true true true true true true true false false ...
  false true false];
uolChange = [0.145666 - 0.298024; 0.341518 - 0.259733; 0.289202 - 0.458551; 0.249652 - 0.424188; ...
  0.434660 - 0.418259; 0.300893 - 0.190590; 0.770454 - 0.415658];
uolFrac = [0.475467; 0.674299; 0.400523; 0.461996; 0.607605; 0.69995; 0.874251];
uolInc = [true false true false true true true];
allenChange = [0.675289 - 0.635800; 0.539034 - 0.615952; 0.647850 - 0.631803; 0.933029 - 0.693471; ...
  0.452149 - 0.481036; 0.476366 - 0.500610; 0.471062 - 0.496322];
allenFrac = [0.400404; 0.235446; 0.384781; 0.781899; 0.446529; 0.473333; 0.467418];
allenInc = [true true false true true true false];
rateChange = nan(numel(areasOI),1);
for iA = 1:numel(areasOI)
  x = frOI(areasOI{iA});
  if ~isempty(x), rateChange(iA) = mean(x(:,2), 'omitnan') - mean(x(:,1), 'omitnan'); end
end
fracAll = [corrOI.positiveSpearmanFractionsMeans(:,1); uolFrac; allenFrac];
changeAll = [rateChange; uolChange; allenChange];
inc = [incMask5 uolInc allenInc];
[rho5, p5] = corrMulti(fracAll(inc)', changeAll(inc)', 'Spearman');
[rho5i, p5i] = corrMulti(corrOI.positiveSpearmanFractionsMeans(incMask5,1)', rateChange(incMask5)', 'Spearman');
say(fid, '\nFig. 5E: Spearman rho = %.2f, p = %.3g (%d area groups: %d IBL + %d UoL + %d Allen); IBL only rho = %.2f, p = %.3g\n', ...
  rho5, p5, sum(inc), sum(incMask5), sum(uolInc), sum(allenInc), rho5i, p5i);

%% Dataset counts
nExp = numel(infraslowData.experimentData);
hasPupil = cellfun(@(e) ~isempty(e.leftPupilAreaSize) || ~isempty(e.rightPupilAreaSize), infraslowData.experimentData);
nUnits = sum(cellfun(@(e) size(e.spikeCounts,1), infraslowData.experimentData));
nUnitsPupil = sum(cellfun(@(e) size(e.spikeCounts,1), infraslowData.experimentData(hasPupil)));
r = infraslowAnalyses.spikingPupilCorr.recordings.rSpearman;
if iscell(r), nClassified = sum(cellfun(@(x) sum(isfinite(x(:))), r)); else, nClassified = sum(isfinite(r(:))); end
say(fid, '\nCounts: %d sessions (%d with pupil), %d units (%d in sessions with pupil), %d units with a pupil correlation\n', ...
  nExp, sum(hasPupil), nUnits, nUnitsPupil, nClassified);
f = fieldnames(infraslowData.experimentData{1});
say(fid, '  experimentData fields: %s\n', strjoin(f', ', '));
fclose(fid);
end


function x = local_per_rec(corrOI, areasOI, acr)
% per-recording positive fractions of one area group (column 1 = Spearman)
i = find(strcmp(areasOI, acr));
x = corrOI.positiveSpearmanFractionsPerRec{i};
if isempty(x), x = []; else, x = x(:,1); x = x(isfinite(x)); end
end


function say(fid, varargin)
% print to the command window and to the output file
fprintf(varargin{:}); fprintf(fid, varargin{:});
end
