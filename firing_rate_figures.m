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

% Load preprocessed data
if ~exist('infraslowData', 'var')
  preprocessedDataFile = fullfile(processedDataFolder, 'bwmPreprocessedData2.mat');
  load(preprocessedDataFile);
end

% Load data analysis results
if ~exist('infraslowAnalyses', 'var')
  analysisResultsFile = fullfile(processedDataFolder, 'bwmAnalysisResults.mat');
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
    bottomFiringRate_exp = zeros(nUnits,nCentiles);
    topFiringRate_exp = zeros(nUnits,nCentiles);
    if ~isempty(pupilAreaSize)
      percentileValues = prctile(pupilAreaSize,percentiles);
      for iCent = 1:nCentiles
        bottomPupilAreaSizeMask = pupilAreaSize <= percentileValues(iCent);
        topPupilAreaSizeMask = pupilAreaSize > percentileValues(end - iCent + 1);
        bottomFiringRate_exp(activeUnits,iCent) = sum(spikeCounts(activeUnits,bottomPupilAreaSizeMask),2)/expDuration;
        topFiringRate_exp(activeUnits,iCent) = sum(spikeCounts(activeUnits,topPupilAreaSizeMask),2)/expDuration;
      end
      bottomFiringRate_exp(~bottomFiringRate_exp) = nan;
      topFiringRate_exp(~topFiringRate_exp) = nan;
      firingRateChange_exp = topFiringRate_exp - bottomFiringRate_exp;
      logBottomFiringRate_exp = log10(bottomFiringRate_exp);
      logTopFiringRate_exp = log10(topFiringRate_exp);
      logFiringRateChange_exp = logTopFiringRate_exp - logBottomFiringRate_exp;
    else
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
figure;
nAreasToVisualise = numel(areaOrderAll);
for iArea = 1:nAreasToVisualise
  areaFiringRates = firingRates{areaOrderAll(iArea)};
  positiveFraction = sum(areaFiringRates(:,1)>=0)/size(areaFiringRates,1);
  firingRateChange = mean(abs(areaFiringRates(:,30)), 'omitnan');
  plot(positiveFraction, firingRateChange, '.', 'MarkerSize',1);
  if iArea == 1
    hold on
  end
end
hold off
%close all

% Create firing rate distribution figures: All types of units combined
nAreas = numel(areasOI);
firingRatesOI = cell(nAreas,1);
binSize = 0.35;
binEdges = -3.5:binSize:2.5;
bins = binEdges(1:end-1)+binSize/2;
for iArea = 1:nAreas
  areaInds = getAreaInds(areasOI{iArea}, infraslowAnalyses.areaSummaries.areaTable);
  firingRatesOI{iArea} = concatenateCells(firingRates(areaInds));
  rSpearman = firingRatesOI{iArea}(:,1);
  bottomFiringRate = firingRatesOI{iArea}(:,20);
  topFiringRate = firingRatesOI{iArea}(:,25);
  bottomFiringRateDistro = histcounts(bottomFiringRate, binEdges);
  topFiringRateDistro = histcounts(topFiringRate, binEdges);

  % fH = figure;
  % plot(bins, bottomFiringRateDistro);
  % hold on
  % plot(bins, topFiringRateDistro);
  % hold off
  % title(areasOI{iArea})
  % legend('0:12.5%','87.5:100%')
  % close(fH);

  % Plot data distributions
  logRates = [bottomFiringRate topFiringRate];
  options.legendLabels = {'constricted','dilated'};
  options.lineStyles = {'-','-'};
  options.lineWidths = [2 2];
  options.markerStyles = {'.','o'};
  options.displayMeans = true;
  options.displayVars = true;
  options.markerVPos = 1.15*max([max(bottomFiringRateDistro) max(topFiringRateDistro)]);
  options.xLim = [-3 1.5];
  options.xLabel = 'Log_{10}(firing rate)';
  options.yLabel = 'Unit count';
  options.saveFig = false;
  options.figSize = 18;
  options.figName = ['Firing_rates_' areasOI{iArea}];
  options.fH = [];
  options.fH = histPlotFR(binEdges, logRates, options);

  % Display stats
  digit = 3;
  statsMean = meanTest(logRates, 'ANOVARM');
  statsVar = varTest(logRates);
  text(max(datamean(logRates))+0.2,options.markerVPos, ['p=' num2str(round(statsMean.p,digit,'significant'))], 'FontSize',20);
  text(max(datamean(logRates))-0.6,0.125*max([max(bottomFiringRateDistro) max(bottomFiringRateDistro)]),...
    ['p=' num2str(round(statsVar.p,digit,'significant'))], 'FontSize',20);

  % Save the figure
  figFolder = './figures/firing_rate_distros';
  if ~isfolder(figFolder)
    mkdir(figFolder);
  end
  label = [3.9 3.1];
  margin = [0.6 0.55];
  width = options.figSize-label(1)-margin(1);
  height = options.figSize-label(2)-margin(2);
  paperSize = resizeFig(options.fH, gca, width, height, label, margin, 0);
  hgsave(options.fH, [figFolder filesep options.figName '.fig']);
  exportFig(options.fH, [figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
  close(options.fH);

  % Draw split +/- data distributions
  options = struct();
  options.xLim = [-3 1.5];
  options.xLabel = 'Log_{10}(firing rate)';
  options.yLabel = 'Unit count';
  options.figFolder = figFolder;
  options.figSize = 18;
  options.pdf = false;
  options.colours = [matlabColours(1); matlabColours(2); matlabColours(1); matlabColours(2)];
  options.lineStyles = {'--','--','-','-'};
  options.markerStyles = {'v','v','^','^'};
  options.legendLabels = {'Constricted -','Constricted +','Dilated -','Dilated +'};
  options.stats = [];
  options.saveFig = false;
  logRates2 = {bottomFiringRate(rSpearman < 0, 1), bottomFiringRate(rSpearman >= 0, 1), ...
    topFiringRate(rSpearman < 0, end), topFiringRate(rSpearman >= 0, end)};
  fH = histPlotFR(binEdges, logRates2, options);
  l = legend;

  % Save the figure
  options.figName = ['Split_firing_rates_' areasOI{iArea}];
  label = [3.9 3.1];
  margin = [0.6 0.55];
  width = options.figSize-label(1)-margin(1);
  height = options.figSize-label(2)-margin(2);
  paperSize = resizeFig(fH, gca, width, height, label, margin, 0);
  hgsave(fH, [figFolder filesep options.figName '.fig']);
  exportFig(fH, [figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
  close(fH);
end


% Additional data
includeExtraData = true;
uolData_areas = {'CA UoL', 'neoCx UoL', 'DG UoL', 'CA-DG UoL', 'RSP UoL', ...
  'S1 UoL', 'Th UoL'};
uolDdata_incMask = [true false true false true true true];
uolData_firingRateChange = [0.145666 - 0.298024;
                            0.341518 - 0.259733;
                            0.289202 - 0.458551;
                            0.249652 - 0.424188;
                            0.434660 - 0.418259;
                            0.300893 - 0.190590;
                            0.770454 - 0.415658];
uolData_positiveFraction = [0.475467; 0.674299; 0.400523; 0.461996; ...
                            0.607605; 0.69995; 0.874251];

allenData_areas = {'CA Allen', 'DG Allen', 'CA-DG Allen', 'Th Allen', ...
  'V1 Allen', 'V2 Allen', 'Cx Allen'};
allenDdata_incMask = [true true false true true true false];
allenData_firingRateChange = [-(0.675289 - 0.635800);
                              0.539034 - 0.615952;
                              -0.06; %0.647850 - 0.631803;
                              0.933029 - 0.693471;
                              0.452149 - 0.481036;
                              0.476366 - 0.500610;
                              0.471062 - 0.496322];
allenData_positiveFraction = [0.400404; 0.235446; 0.384781; 0.781899; ...
                              0.446529; 0.473333; 0.467418];

% Firing rate summary figures: Firing rate change
incMask = [false false false false false true true false true true true ...
  true true true false false false false false false false false false ...
  false false false false false false false false true true true true ...
  false true true true true true true true true true true];
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

includeAreaInds = [1:14 17:21 24 29:numel(incMask)];
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
      rateChange = mean(topFiringRate, 'omitnan') - mean(bottomFiringRate, 'omitnan');
      rateChange_alt = mean(topFiringRate - bottomFiringRate, 'omitnan');
    end

    if contains(areasOI(areaOrderOI(iArea)), {'Th', 'VP', 'PO', 'LG', 'LP'}) && ~strcmpi(areasOI(areaOrderOI(iArea)), 'STh')
      lineColour = 'b';
    elseif contains(areasOI(areaOrderOI(iArea)), {'Cx', 'SS', 'RSP', 'VIS'})
      lineColour = 'g';
    elseif contains(areasOI(areaOrderOI(iArea)), {'Hp', 'CA', 'DG'})
      lineColour = 'm';
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
xticklabels(areasOI_full(areaOrderOI(ismember(areaOrderOI, includeAreaInds))));
ylabel('Fraction')

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
set(get(ax, 'YAxis'), 'FontWeight', 'bold');
xlim([0.5 lineCount+0.5])
%ylim([0.35 0.85])


% Firing rate summary figures: Firing rate change correlation
% Aggregate firing rates
rateChange = nan(nAreasInit,1);
for iArea = 1:nAreasInit
  areaInds = getAreaInds(areasOI{iArea}, infraslowAnalyses.areaSummaries.areaTable);
  firingRatesOI{iArea} = concatenateCells(firingRates(areaInds)); %#ok<*SAGROW>
  bottomFiringRate = firingRatesOI{iArea}(:,20);
  topFiringRate = firingRatesOI{iArea}(:,25);
  rateChange(iArea) = mean(topFiringRate, 'omitnan') - mean(bottomFiringRate, 'omitnan');
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
set(get(ax, 'XAxis'), 'FontWeight', 'bold');
set(get(ax, 'YAxis'), 'FontWeight', 'bold');