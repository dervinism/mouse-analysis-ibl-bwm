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
say(fid, '\nt-tests (recordings as units, pooled variance); pHolm = Holm-Bonferroni across these %d tests\n', size(tests,1));
pRaw = nan(size(tests,1), 1);
for i = 1:size(tests,1)
  a = cell2mat(cellfun(perRec, tests{i,2}, 'UniformOutput', false)');
  b = cell2mat(cellfun(perRec, tests{i,4}, 'UniformOutput', false)');
  [~, pRaw(i)] = ttest2(a, b);
end
[~, o] = sort(pRaw); m = numel(pRaw); pHolm = nan(m, 1); pHolm(o) = min(1, cummax((m:-1:1)' .* pRaw(o)));
for i = 1:size(tests,1)
  a = cell2mat(cellfun(perRec, tests{i,2}, 'UniformOutput', false)');
  b = cell2mat(cellfun(perRec, tests{i,4}, 'UniformOutput', false)');
  [~, p, ~, st] = ttest2(a, b);
  say(fid, '  %s vs %s: %.3f±%.3f (n %d) vs %.3f±%.3f (n %d), t(%d) = %.2f, p = %.2g, pHolm = %.2g\n', tests{i,1}, tests{i,3}, ...
    mean(a), tinv(0.975, numel(a)-1)*std(a)/sqrt(numel(a)), numel(a), mean(b), ...
    tinv(0.975, numel(b)-1)*std(b)/sqrt(numel(b)), numel(b), st.df, st.tstat, p, pHolm(i));
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
    [~, lr] = stateRates(counts, pupil, expData.samplingRate, 12.5);   % time-in-state denominator, zero -> -3 (as UoL/Allen)
    lr = reshape(lr, [], 2); lr(~active,:) = nan;
    fr(unitMask,:) = lr;
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

% Median-based statistics (2026-11-08): medians, IQR, paired Wilcoxon signed-rank, Brown-Forsythe
load(fullfile(processedDataFolder, 'uolAllen_fig5.mat'), 'uolAllenFig5');
ua = @(lbl) uolAllenFig5(strcmp({uolAllenFig5.label}, lbl));
say(fid, '\nFig. 5 median statistics (log10 rate, silent in a state = -3; constricted vs dilated)\n');
fig5med(fid, 'IBL DG', logRates);
fig5med(fid, 'UoL Th', ua('Th UoL').logRates);
fig5med(fid, 'Allen Th', ua('Th Allen').logRates);
fig5med(fid, 'Allen VIS (Cx)', ua('Cx Allen').logRates);

incMask5 = [false false false false false true true false true true true ...
  true true true false false false false false false false false false ...
  false false false false false false false false false false false false ...
  true true true true true true true true true true true true false false ...
  false true false];
uolLbl = {'CA UoL', 'neoCx UoL', 'DG UoL', 'CA-DG UoL', 'RSP UoL', 'S1 UoL', 'Th UoL'};
allenLbl = {'CA Allen', 'DG Allen', 'CA-DG Allen', 'Th Allen', 'V1 Allen', 'V2 Allen', 'Cx Allen'};
medChange = @(L) median(L(:,2), 'omitnan') - median(L(:,1), 'omitnan');
uolChange = cellfun(@(l) medChange(ua(l).logRates), uolLbl)';
uolFrac = cellfun(@(l) ua(l).positiveFractionTable, uolLbl)';
uolInc = [true false true false true true true];
allenChange = cellfun(@(l) medChange(ua(l).logRates), allenLbl)';
allenFrac = cellfun(@(l) ua(l).positiveFractionTable, allenLbl)';
allenInc = [true true false true true true false];
rateChange = nan(numel(areasOI),1);
for iA = 1:numel(areasOI)
  x = frOI(areasOI{iA});
  if ~isempty(x), rateChange(iA) = medChange(x); end   % difference of medians (2026-11-08)
end
fracAll = [corrOI.positiveSpearmanFractionsMeans(:,1); uolFrac; allenFrac];
changeAll = [rateChange; uolChange; allenChange];
inc = [incMask5 uolInc allenInc];
[rho5, p5] = corrMulti(fracAll(inc)', changeAll(inc)', 'Spearman');
[rho5i, p5i] = corrMulti(corrOI.positiveSpearmanFractionsMeans(incMask5,1)', rateChange(incMask5)', 'Spearman');
for k = 1:numel(uolLbl), say(fid, '  %-12s median change %+.3f, positive fraction %.3f\n', uolLbl{k}, uolChange(k), uolFrac(k)); end
for k = 1:numel(allenLbl), say(fid, '  %-12s median change %+.3f, positive fraction %.3f\n', allenLbl{k}, allenChange(k), allenFrac(k)); end
say(fid, '\nFig. 5E: Spearman rho = %.2f, p = %.3g (%d area groups: %d IBL + %d UoL + %d Allen); IBL only rho = %.2f, p = %.3g\n', ...
  rho5, p5, sum(inc), sum(incMask5), sum(uolInc), sum(allenInc), rho5i, p5i);

%% Dataset counts
nExp = numel(infraslowData.experimentData);
hasPupil = cellfun(@(e) (isfield(e, 'leftPupilAreaSize') && ~isempty(e.leftPupilAreaSize)) || ...
  (isfield(e, 'rightPupilAreaSize') && ~isempty(e.rightPupilAreaSize)), infraslowData.experimentData);   % some sessions lack the fields
nUnits = sum(cellfun(@(e) local_n(e), infraslowData.experimentData));
nUnitsPupil = sum(cellfun(@(e) local_n(e), infraslowData.experimentData(hasPupil)));
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


function n = local_n(e)
% unit count of a session (0 for sessions without spike data)
if isfield(e, 'spikeCounts'), n = size(e.spikeCounts, 1); else, n = 0; end
end


function fig5med(fid, name, L)
% medians, IQR, paired Wilcoxon signed-rank (normal approximation) and Brown-Forsythe test for one area
ok = all(~isnan(L), 2); c = L(ok,1); d = L(ok,2); n = sum(ok);
[pW, ~, sW] = signrank(c, d, 'method','approximate');
[pBF, sBF] = vartestn([c d], 'TestType','BrownForsythe', 'Display','off');
say(fid, '  %-15s n %d | median %.3f -> %.3f (IQR %.3f-%.3f -> %.3f-%.3f; widths %.2f -> %.2f) | median paired change %+.3f\n', ...
  name, n, median(c), median(d), prctile(c,25), prctile(c,75), prctile(d,25), prctile(d,75), iqr(c), iqr(d), median(d - c));
say(fid, '  %-15s Wilcoxon z = %.2f, p = %.2g | Brown-Forsythe F(%d, %d) = %.1f, p = %.2g | mean %.3f -> %.3f\n', '', ...
  sW.zval, pW, sBF.df(1), sBF.df(2), sBF.fstat, pBF, mean(c), mean(d));
end
