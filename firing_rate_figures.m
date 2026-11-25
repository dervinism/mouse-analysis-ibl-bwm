% Figures quantifying firing rates and firing rate changes in correlated
% and anticorrelated cells in different brain areas

clear dataAugmented

% Load parameters
params
alpha = 0.05; % Significance level
excludeMovement = false;
includeAllAreas = true;
percentiles = [12.5 25 100/3 37.5 50 62.5 200/3 75 87.5];
firingRateThr = 50/3600;
if ~exist('useMedians', 'var') || isempty(useMedians)
  useMedians = true;   % 2026-11-08 (co-author request): medians, paired Wilcoxon, Brown-Forsythe instead of means, t, Levene
end
% UoL and Allen per-unit data for Figure 5 (extract_uol_allen_fig5.m; Allen = Aug/Sep 2025 run)
load(fullfile(processedDataFolder, 'uolAllen_fig5.mat'), 'uolAllenFig5');

% Load preprocessed data
if ~exist('infraslowData', 'var')
  preprocessedDataFile = fullfile(processedDataFolder, 'bwmPreprocessedData2.mat');
  load(preprocessedDataFile);
end

% Load data analysis results
if ~exist('infraslowAnalyses', 'var')
  analysisResultsFile = fullfile(processedDataFolder, ['bwmAnalysisResults' resultsTag '.mat']);
  load(analysisResultsFile);
end

% Compile firing rates
nAreas = size(areaLabels,1);
nCentiles = ceil(numel(percentiles)/2);
firingRates = cell(nAreas,1);
for iArea = 1:nAreas
  areaInds = infraslowAnalyses.areaSummaries.groupedUnitInds{iArea};
  firingRates_area = nan(size(areaInds,1),34);
  expInds = unique(areaInds(:,1));
  nExps = numel(expInds);
  for iExp = 1:nExps
    unitMask = areaInds(:,1) == expInds(iExp);
    nUnits = sum(unitMask);

    % Get experiment data
    expData = infraslowData.experimentData{expInds(iExp)};

    % Get pupil correlation data
    if isempty(infraslowAnalyses.spikingPupilCorr.singleAreas.rSpearman{iArea})
      rSpearman = nan(nUnits,1);
      pvalSpearman = nan(nUnits,1);
    else
      rSpearman = infraslowAnalyses.spikingPupilCorr.singleAreas.rSpearman{iArea}(unitMask);
      pvalSpearman = infraslowAnalyses.spikingPupilCorr.singleAreas.pvalSpearman{iArea}(unitMask);
      if isempty(rSpearman)
        rSpearman = nan(nUnits,1);
        pvalSpearman = nan(nUnits,1);
      end
    end

    % Get the pupil area size
    if ~isempty(expData.leftPupilAreaSize) && ~isempty(expData.rightPupilAreaSize)
      pupilAreaSize = mean([expData.leftPupilAreaSize; expData.rightPupilAreaSize]);
    elseif ~isempty(expData.leftPupilAreaSize)
      pupilAreaSize = expData.leftPupilAreaSize;
    elseif ~isempty(expData.rightPupilAreaSize)
      pupilAreaSize = expData.rightPupilAreaSize;
    else
      pupilAreaSize = [];
    end

    % Calculate firing rates
    spikeCounts = expData.spikeCounts(areaInds(unitMask,2),:);
    nSamples = size(spikeCounts,2);
    samplingRate = expData.samplingRate;
    expDuration = nSamples/samplingRate;
    firingRates_exp = sum(spikeCounts,2)/expDuration;
    firingRates_exp(~firingRates_exp) = nan;
    logFiringRates_exp = log10(firingRates_exp);

    activeUnits = firingRates_exp > firingRateThr;
    if ~isempty(pupilAreaSize)
      % State rates as in the UoL/Allen analysis (2026-11-08 fix): spike count in the state / time spent
      % in it (was / whole-recording duration, ~8x too low), silent units kept at a zero rate with
      % log10 = -3 (were dropped from that state only, unpairing the comparison). Inactive units: NaN.
      [stRates, stLogRates] = stateRates(spikeCounts, pupilAreaSize, samplingRate, percentiles(1:nCentiles));
      bottomFiringRate_exp = stRates(:,:,1); topFiringRate_exp = stRates(:,:,2);
      logBottomFiringRate_exp = stLogRates(:,:,1); logTopFiringRate_exp = stLogRates(:,:,2);
      bottomFiringRate_exp(~activeUnits,:) = nan; topFiringRate_exp(~activeUnits,:) = nan;
      logBottomFiringRate_exp(~activeUnits,:) = nan; logTopFiringRate_exp(~activeUnits,:) = nan;
      firingRateChange_exp = topFiringRate_exp - bottomFiringRate_exp;
      logFiringRateChange_exp = logTopFiringRate_exp - logBottomFiringRate_exp;
    else
      bottomFiringRate_exp = nan(nUnits,nCentiles);   % no pupil: no state rates (was 0)
      topFiringRate_exp = nan(nUnits,nCentiles);
      firingRateChange_exp = nan(nUnits,nCentiles);
      logBottomFiringRate_exp = nan(nUnits,nCentiles);
      logTopFiringRate_exp = nan(nUnits,nCentiles);
      logFiringRateChange_exp = nan(nUnits,nCentiles);
    end
    firingRates_area(unitMask,:) = [rSpearman pvalSpearman ...               % 1 2
      firingRates_exp logFiringRates_exp ...                                 % 3 4
      bottomFiringRate_exp topFiringRate_exp firingRateChange_exp ...        % 5:9  10:14  15:19
      logBottomFiringRate_exp logTopFiringRate_exp logFiringRateChange_exp]; % 20:24  25:29  30:34
  end
  firingRates{iArea} = firingRates_area;
end

% Sort areas large to small positive fraction
[sortedAreasAll, areaOrderAll] = sort( ...
  infraslowAnalyses.spikingPupilCorr.singleAreas.positiveSpearmanFractionsMeans(:,1), 'descend');
areaOrderAll = areaOrderAll(~isnan(sortedAreasAll));
%areaOrderAll = intersect(intersect(areaOrderAll, find(contains(areaLabels(:,1), {'Th','Cx','Hp'})), 'stable'), nonemptyEntries, 'stable');

% Compare overall firing rates
figure;
nAreasToVisualise = numel(areaOrderAll);
for iArea = 1:nAreasToVisualise
  areaFiringRates = firingRates{areaOrderAll(iArea)};
  positiveUnitMask = areaFiringRates(:,1) > 0;
  meanPositiveFiringRate = mean(areaFiringRates(positiveUnitMask,4), 'omitnan');
  meanNegativeFiringRate = mean(areaFiringRates(~positiveUnitMask,4), 'omitnan');
  meanFiringRate = mean(areaFiringRates(:,4), 'omitnan');
  %plot(iArea, meanPositiveFiringRate, '.r', 'MarkerSize',10);
  if iArea == 1
    hold on
  end
  plot(iArea, meanNegativeFiringRate, '.b', 'MarkerSize',10);
  %plot(iArea, meanFiringRate, '.k', 'MarkerSize',10);
end
hold off

% Compare firing rate changes
figure;
nAreasToVisualise = numel(areaOrderAll);
for iArea = 1:nAreasToVisualise
  areaFiringRates = firingRates{areaOrderAll(iArea)};
  positiveUnitMask = areaFiringRates(:,1) > 0;
  meanPositiveFiringRateChange = mean(abs(areaFiringRates(positiveUnitMask,30)), 'omitnan');
  meanNegativeFiringRateChange = mean(abs(areaFiringRates(~positiveUnitMask,30)), 'omitnan');
  meanFiringRateChange = mean(abs(areaFiringRates(:,30)), 'omitnan');
  plot(iArea, meanPositiveFiringRateChange, '.r', 'MarkerSize',10);
  if iArea == 1
    hold on
  end
  %plot(iArea, meanNegativeFiringRateChange, '.b', 'MarkerSize',10);
  %plot(iArea, meanFiringRateChange, '.k', 'MarkerSize',10);
end
hold off

% Correlate firing rates and firing rate changes
figure;
nAreasToVisualise = numel(areaOrderAll);
for iArea = 1:nAreasToVisualise
  areaFiringRates = firingRates{areaOrderAll(iArea)};
  firingRate = areaFiringRates(:,4);
  firingRateChange = abs(areaFiringRates(:,30));
  plot(firingRate, firingRateChange, '.', 'MarkerSize',1);
  if iArea == 1
    hold on
  end
end
hold off

% Correlate positive fraction and firing rate changes
nAreasToVisualise = numel(areaOrderAll);
for type = 1:3
  figure;
  for iArea = 1:nAreasToVisualise
    areaFiringRates = firingRates{areaOrderAll(iArea)};
    positiveUnitMask = areaFiringRates(:,1)>=0;
    positiveFraction = sum(positiveUnitMask)/size(areaFiringRates,1);
    if type == 1
      unitsOI = true(size(positiveUnitMask));
    elseif type == 2
      unitsOI = positiveUnitMask;
    elseif type == 3
      unitsOI = ~positiveUnitMask;
    end
    firingRateChange = mean(abs(areaFiringRates(unitsOI,30)), 'omitnan');
    plot(positiveFraction, firingRateChange, '.', 'MarkerSize',5);
    if iArea == 1
      hold on
    end
  end
  hold off
end
%close all

% Create firing rate distribution figures: All types of units combined
nAreas = numel(areasOI);
firingRatesOI = cell(nAreas,1);
binSize = 0.35;
binEdges = -3.5:binSize:2.5;
bins = binEdges(1:end-1)+binSize/2;
distroFolder = fullfile(figFolder, 'firing_rate_distros');
if ~isfolder(distroFolder)
  mkdir(distroFolder);
end
for iArea = 1:nAreas
  areaInds = getAreaInds(areasOI{iArea}, infraslowAnalyses.areaSummaries.areaTable);
  firingRatesOI{iArea} = concatenateCells(firingRates(areaInds));
  rSpearman = firingRatesOI{iArea}(:,1);
  bottomFiringRate = firingRatesOI{iArea}(:,20);
  topFiringRate = firingRatesOI{iArea}(:,25);
  plotRateDistros(bottomFiringRate, topFiringRate, rSpearman, areasOI{iArea}, distroFolder, binEdges, useMedians);
end
% UoL and Allen areas shown in Figure 5 (A: UoL thalamus or Allen thalamus; B: Allen visual cortex)
for lbl = {'Th UoL', 'Th Allen', 'Cx Allen'}
  e = uolAllenFig5(strcmp({uolAllenFig5.label}, lbl{1}));
  plotRateDistros(e.logRates(:,1), e.logRates(:,2), e.rho, strrep(lbl{1}, ' ', '_'), distroFolder, binEdges, useMedians);
end


% Additional data
includeExtraData = true;
uolData_areas = {'CA UoL', 'neoCx UoL', 'DG UoL', 'CA-DG UoL', 'RSP UoL', ...
  'S1 UoL', 'Th UoL'};
uolDdata_incMask = [true false true false true true true];
allenData_areas = {'CA Allen', 'DG Allen', 'CA-DG Allen', 'Th Allen', ...
  'V1 Allen', 'V2 Allen', 'Cx Allen'};
allenDdata_incMask = [true true false true true true false];
% Rate changes and positive fractions of the UoL and Allen areas (2026-11-08): rate changes computed from the
% per-unit data in uolAllen_fig5.mat (previously hard-coded, see firing_rate_figures_before_medians_20261108.m.txt);
% positive fractions = the per-recording means behind Tables 2 and 3 (see extract_uol_allen_fig5.m).
uolData_firingRateChange = nan(numel(uolData_areas), 1); uolData_positiveFraction = nan(numel(uolData_areas), 1);
for k = 1:numel(uolData_areas)
  e = uolAllenFig5(strcmp({uolAllenFig5.label}, uolData_areas{k}));
  uolData_firingRateChange(k) = stateChange(e.logRates(:,1), e.logRates(:,2), useMedians);
  uolData_positiveFraction(k) = e.positiveFractionTable;
end
allenData_firingRateChange = nan(numel(allenData_areas), 1); allenData_positiveFraction = nan(numel(allenData_areas), 1);
for k = 1:numel(allenData_areas)
  e = uolAllenFig5(strcmp({uolAllenFig5.label}, allenData_areas{k}));
  allenData_firingRateChange(k) = stateChange(e.logRates(:,1), e.logRates(:,2), useMedians);
  allenData_positiveFraction(k) = e.positiveFractionTable;
end

% Firing rate summary figures: Firing rate change
incMask = [false false false false false true true false true true true ...
  true true true false false false false false false false false false ...
  false false false false false false false false false false false false ...
  true true true true true true true true true true true true false false ...
  false true false];
positiveSpearmanFractionsMeans = ...
  infraslowAnalyses.spikingPupilCorr.areasOI.positiveSpearmanFractionsMeans(:,1);
if includeExtraData
  if ~exist('dataAugmented', 'var')
    nAreasInit = numel(areasOI);
    areasOI_full = [areasOI_full; uolData_areas'; allenData_areas'];
    areasOI = [areasOI; uolData_areas'; allenData_areas'];
  end
  positiveSpearmanFractionsMeans = [positiveSpearmanFractionsMeans;
    uolData_positiveFraction; allenData_positiveFraction];
  incMask = [incMask uolDdata_incMask allenDdata_incMask];
  dataAugmented = true;
end
[sortedAreasOI, areaOrderOI] = sort( ...
  positiveSpearmanFractionsMeans, 'descend');
areaOrderOI = areaOrderOI(~isnan(sortedAreasOI));

%includeAreaInds = [1:14 17:21 24 29:numel(incMask)];
%includeAreaInds = [1:14 17:18, 24 29:32 36:numel(incMask)];
includeAreaInds = [1:14 17:18, 24 29:32 36:47 51 53:numel(incMask)];
nAreas = numel(areaOrderOI);
fontSize = 18;
fH = figure;
lineCount = 0;
for iArea = 1:nAreas
  if ismember(areaOrderOI(iArea), includeAreaInds)
    if ismember(areasOI_full{areaOrderOI(iArea)}, uolData_areas)
      rateChange = uolData_firingRateChange( ...
        ismember(uolData_areas, areasOI_full{areaOrderOI(iArea)}));
    elseif ismember(areasOI_full{areaOrderOI(iArea)}, allenData_areas)
      rateChange = allenData_firingRateChange( ...
        ismember(allenData_areas, areasOI_full{areaOrderOI(iArea)}));
    else
      areaInds = getAreaInds(areasOI{areaOrderOI(iArea)}, infraslowAnalyses.areaSummaries.areaTable);
      firingRatesOI{iArea} = concatenateCells(firingRates(areaInds)); %#ok<*SAGROW>
      bottomFiringRate = firingRatesOI{iArea}(:,20);
      topFiringRate = firingRatesOI{iArea}(:,25);
      rateChange = stateChange(bottomFiringRate, topFiringRate, useMedians);
    end

    if contains(areasOI(areaOrderOI(iArea)), {'Th', 'VP', 'PO', 'LG', 'LP'}) && ~strcmpi(areasOI(areaOrderOI(iArea)), 'STh')
      lineColour = [0.620 0.855 0.957]; %'b';
    elseif contains(areasOI(areaOrderOI(iArea)), {'Cx', 'SS', 'RSP', 'VIS'})
      lineColour = [0.843 0.788 0.361]; %'g';
    elseif contains(areasOI(areaOrderOI(iArea)), {'Hp', 'CA', 'DG'})
      lineColour = [0.867 0.561 0.891]; %'m';
    else
      lightGrey = 211;
      lineColour = [lightGrey lightGrey lightGrey]./256;
    end
    lineCount = lineCount + 1;
    x = [lineCount lineCount];
    y = [0 rateChange];
    %y = [0 rateChange_alt];
    plot(x, y, '-', 'color',lineColour, 'LineWidth',1.5);
    if lineCount == 1
      hold on
    end
  end
end
plot([0.5 lineCount+0.5], [0 0], ':');
hold off
xticks(1:lineCount);
areasOI_full{56} = 'Hp UoL';
areasOI_full{62} = 'Hp Allen';
xticklabels(areasOI_full(areaOrderOI(ismember(areaOrderOI, includeAreaInds))));
xtickangle(45)
ylabel('Log_{10}(Firing rate) change')

% Tidy the figure
set(fH, 'Color', 'white');
ax = gca;
set(ax, 'box', 'off');
set(ax, 'TickDir', 'out');
yTicks = get(ax, 'YTick');
if numel(yTicks) > 8
  set(ax, 'YTick', yTicks(1:2:end));
end
ax.FontSize = fontSize - 4;
%set(get(ax, 'YAxis'), 'FontWeight', 'bold');
xlim([0.5 lineCount+0.5])
%ylim([0.35 0.85])


% % Firing rate summary figures: Firing rate variance
% positiveSpearmanFractionsMeans = ...
%   infraslowAnalyses.spikingPupilCorr.areasOI.positiveSpearmanFractionsMeans(:,1);
% [sortedAreasOI, areaOrderOI] = sort( ...
%   positiveSpearmanFractionsMeans, 'descend');
% areaOrderOI = areaOrderOI(~isnan(sortedAreasOI));
% 
% includeAreaInds = [1:14 17:21 24 29:numel(areaOrderOI)];
% nAreas = numel(areaOrderOI);
% fontSize = 18;
% fH = figure;
% lineCount = 0;
% for iArea = 1:nAreas
%   if ismember(areaOrderOI(iArea), includeAreaInds)
%     areaInds = getAreaInds(areasOI{areaOrderOI(iArea)}, infraslowAnalyses.areaSummaries.areaTable);
%     firingRatesOI{iArea} = concatenateCells(firingRates(areaInds)); %#ok<*SAGROW>
%     firingRateVar = var(firingRatesOI{iArea}(:,4));
% 
%     if contains(areasOI(areaOrderOI(iArea)), {'Th', 'VP', 'PO', 'LG', 'LP'}) && ~strcmpi(areasOI(areaOrderOI(iArea)), 'STh')
%       lineColour = 'b';
%     elseif contains(areasOI(areaOrderOI(iArea)), {'Cx', 'SS', 'RSP', 'VIS'})
%       lineColour = 'g';
%     elseif contains(areasOI(areaOrderOI(iArea)), {'Hp', 'CA', 'DG'})
%       lineColour = 'm';
%     else
%       lightGrey = 211;
%       lineColour = [lightGrey lightGrey lightGrey]./256;
%     end
%     lineCount = lineCount + 1;
%     x = [lineCount lineCount];
%     y = [0 firingRateVar];
%     plot(x, y, '-', 'color',lineColour, 'LineWidth',1.5);
%     if lineCount == 1
%       hold on
%     end
%   end
% end
% plot([0.5 lineCount+0.5], [0 0], ':');
% hold off
% xticks(1:lineCount);
% xticklabels(areasOI_full(areaOrderOI(ismember(areaOrderOI, includeAreaInds))));
% ylabel('Fraction')
% 
% % Tidy the figure
% set(fH, 'Color', 'white');
% ax = gca;
% set(ax, 'box', 'off');
% set(ax, 'TickDir', 'out');
% yTicks = get(ax, 'YTick');
% if numel(yTicks) > 8
%   set(ax, 'YTick', yTicks(1:2:end));
% end
% ax.FontSize = fontSize - 4;
% set(get(ax, 'YAxis'), 'FontWeight', 'bold');
% xlim([0.5 lineCount+0.5])
% %ylim([0.35 0.85])


% Firing rate summary figures: Firing rate change correlation
% Aggregate firing rates
rateChange = nan(nAreasInit,1);
for iArea = 1:nAreasInit
  areaInds = getAreaInds(areasOI{iArea}, infraslowAnalyses.areaSummaries.areaTable);
  firingRatesOI{iArea} = concatenateCells(firingRates(areaInds)); %#ok<*SAGROW>
  bottomFiringRate = firingRatesOI{iArea}(:,20);
  topFiringRate = firingRatesOI{iArea}(:,25);
  rateChange(iArea) = stateChange(bottomFiringRate, topFiringRate, useMedians);
end
if includeExtraData
  rateChange = [rateChange; uolData_firingRateChange; allenData_firingRateChange];
end

% Correlate firing rate change and the positive fraction
[rSpearman, pvalSpearman] = corrMulti( ...
  positiveSpearmanFractionsMeans(incMask)', rateChange(incMask)', 'Spearman');

% Fit the data
[~, slope, coefficients] = fitLine( ...
  positiveSpearmanFractionsMeans(incMask), rateChange(incMask), type='linear-linear');

% Draw the figure
fontSize = 18;
fH = figure;
plot(positiveSpearmanFractionsMeans(incMask), rateChange(incMask), ...
  'o', 'MarkerSize',5, 'MarkerEdgeColor','b', 'MarkerFaceColor','b');
xlabel('Positive fraction');
ylabel('Log_{10}(Firing rate) change')

% Draw the line fit
xLim = xlim;
yFit = coefficients(1)*xLim + coefficients(2);
hold on; plot(xLim, yFit, ':k'); hold off

% Add text
yLim = ylim;
xAxisSize = diff(xLim);
yAxisSize = diff(yLim);
xPos = xLim(1) + 0.05*xAxisSize;
yPos = yLim(1) + 0.95*yAxisSize;
txtStr = ['rho = ' num2str(round(rSpearman, 1,'decimals')) ...
  '  p = ' num2str(round(pvalSpearman, 1,'significant'))];
text(xPos, yPos, txtStr);
xlim(xLim);
ylim(yLim);

% Tidy the figure
set(fH, 'Color', 'white');
ax = gca;
set(ax, 'box', 'off');
set(ax, 'TickDir', 'out');
ax.FontSize = fontSize - 4;
%set(get(ax, 'XAxis'), 'FontWeight', 'bold');
%set(get(ax, 'YAxis'), 'FontWeight', 'bold');


%% Local functions
function d = stateChange(bottom, top, useMedians)
% Change of log10 firing rate from the constricted to the dilated state: difference of medians (useMedians) or means.
if useMedians
  d = median(top, 'omitnan') - median(bottom, 'omitnan');
else
  d = mean(top, 'omitnan') - mean(bottom, 'omitnan');
end
end


function plotRateDistros(bottomFiringRate, topFiringRate, rSpearman, areaName, distroFolder, binEdges, useMedians)
% Constricted vs dilated log firing rate distributions (all units; and split by the sign of rho), with the
% location test (paired t / Wilcoxon signed-rank) and the spread test (Levene quadratic / Brown-Forsythe).
bottomFiringRateDistro = histcounts(bottomFiringRate, binEdges);
topFiringRateDistro = histcounts(topFiringRate, binEdges);
logRates = [bottomFiringRate topFiringRate];
options.legendLabels = {'constricted','dilated'};
options.lineStyles = {'-','-'};
options.lineWidths = [2 2];
options.markerStyles = {'.','o'};
options.displayMeans = true;
options.displayVars = true;
if useMedians, options.centre = 'median'; end
options.markerVPos = 1.15*max([max(bottomFiringRateDistro) max(topFiringRateDistro)]);
options.xLim = [-3 2];   % widened from 1.5 (user, 2026-11-08): the fastest units were cut off
options.xLabel = 'Log_{10}(firing rate)';
options.yLabel = 'Unit count';
options.saveFig = false;
options.figSize = 18;
options.figName = ['Firing_rates_' areaName];
options.fH = [];
options.fH = histPlotFR(binEdges, logRates, options);

% Display stats
digit = 3;
ok = all(~isnan(logRates), 2);
if useMedians
  pLoc = signrank(logRates(ok,1), logRates(ok,2), 'method','approximate');
  pSpread = vartestn(logRates(ok,:), 'TestType','BrownForsythe', 'Display','off');
  centres = median(logRates(ok,:));
else
  statsMean = meanTest(logRates, 'ANOVARM'); pLoc = statsMean.p;
  statsVar = varTest(logRates); pSpread = statsVar.p;
  centres = datamean(logRates);
end
if useMedians && max(centres) + 1.6 > options.xLim(2)   % no room right of the markers: right-aligned to their left
  text(min(centres)-0.15,options.markerVPos, ['p=' num2str(round(pLoc,digit,'significant'))], 'FontSize',20, 'HorizontalAlignment','right');
else
  text(max(centres)+0.2,options.markerVPos, ['p=' num2str(round(pLoc,digit,'significant'))], 'FontSize',20);
end
text(max(centres)-0.6,0.125*max([max(bottomFiringRateDistro) max(bottomFiringRateDistro)]),...
  ['p=' num2str(round(pSpread,digit,'significant'))], 'FontSize',20);

% Save the figure
label = [3.9 3.1];
margin = [0.6 0.55];
width = options.figSize-label(1)-margin(1);
height = options.figSize-label(2)-margin(2);
paperSize = resizeFig(options.fH, gca, width, height, label, margin, 0);
hgsave(options.fH, [distroFolder filesep options.figName '.fig']);
exportFig(options.fH, [distroFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
close(options.fH);

% Draw split +/- data distributions (units without a pupil correlation are left out)
options = struct();
options.xLim = [-3 2];
options.xLabel = 'Log_{10}(firing rate)';
options.yLabel = 'Unit count';
options.figFolder = distroFolder;
options.figSize = 18;
options.pdf = false;
options.colours = [matlabColours(1); matlabColours(2); matlabColours(1); matlabColours(2)];
options.lineStyles = {'--','--','-','-'};
options.markerStyles = {'v','v','^','^'};
options.legendLabels = {'Constricted -','Constricted +','Dilated -','Dilated +'};
options.stats = [];
options.saveFig = false;
if useMedians, options.centre = 'median'; end
logRates2 = {bottomFiringRate(rSpearman < 0, 1), bottomFiringRate(rSpearman >= 0, 1), ...
  topFiringRate(rSpearman < 0, end), topFiringRate(rSpearman >= 0, end)};
fH = histPlotFR(binEdges, logRates2, options);
legend;
options.figName = ['Split_firing_rates_' areaName];
paperSize = resizeFig(fH, gca, width, height, label, margin, 0);
hgsave(fH, [distroFolder filesep options.figName '.fig']);
exportFig(fH, [distroFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
close(fH);
end
