% % Figures quantifying firing rates and firing rate changes in correlated
% % and anticorrelated cells in different brain areas
% 
% % Load parameters
% params
% alpha = 0.05; % Significance level
% excludeMovement = false;
% includeAllAreas = true;
% percentiles = [12.5 25 100/3 37.5 50 62.5 200/3 75 87.5];
% firingRateThr = 50/3600;
% 
% % Load preprocessed data
% if ~exist('infraslowData', 'var')
%   preprocessedDataFile = fullfile(processedDataFolder, 'bwmPreprocessedData2.mat');
%   load(preprocessedDataFile);
% end
% 
% % Load data analysis results
% if ~exist('infraslowAnalyses', 'var')
%   analysisResultsFile = fullfile(processedDataFolder, 'bwmAnalysisResults.mat');
%   load(analysisResultsFile);
% end
% 
% % Compile firing rates
% nAreas = size(areaLabels,1);
% firingRates = cell(nAreas,1);
% for iArea = 1:nAreas
%   areaInds = infraslowAnalyses.areaSummaries.groupedUnitInds{iArea};
%   firingRates_area = nan(size(areaInds,1),34);
%   expInds = unique(areaInds(:,1));
%   nExps = numel(expInds);
%   for iExp = 1:nExps
%     unitMask = areaInds(:,1) == expInds(iExp);
%     nUnits = sum(unitMask);
% 
%     % Get experiment data
%     expData = infraslowData.experimentData{expInds(iExp)};
% 
%     % Get pupil correlation data
%     if isempty(infraslowAnalyses.spikingPupilCorr.singleAreas.rSpearman{iArea})
%       rSpearman = nan(nUnits,1);
%       pvalSpearman = nan(nUnits,1);
%     else
%       rSpearman = infraslowAnalyses.spikingPupilCorr.singleAreas.rSpearman{iArea}(unitMask);
%       pvalSpearman = infraslowAnalyses.spikingPupilCorr.singleAreas.pvalSpearman{iArea}(unitMask);
%       if isempty(rSpearman)
%         rSpearman = nan(nUnits,1);
%         pvalSpearman = nan(nUnits,1);
%       end
%     end
% 
%     % Get the pupil area size
%     if ~isempty(expData.leftPupilAreaSize) && ~isempty(expData.rightPupilAreaSize)
%       pupilAreaSize = mean([expData.leftPupilAreaSize; expData.rightPupilAreaSize]);
%     elseif ~isempty(expData.leftPupilAreaSize)
%       pupilAreaSize = expData.leftPupilAreaSize;
%     elseif ~isempty(expData.rightPupilAreaSize)
%       pupilAreaSize = expData.rightPupilAreaSize;
%     else
%       pupilAreaSize = [];
%     end
% 
%     % Calculate firing rates
%     spikeCounts = expData.spikeCounts(areaInds(unitMask,2),:);
%     nSamples = size(spikeCounts,2);
%     samplingRate = expData.samplingRate;
%     expDuration = nSamples/samplingRate;
%     firingRates_exp = sum(spikeCounts,2)/expDuration;
%     firingRates_exp(~firingRates_exp) = nan;
%     logFiringRates_exp = log10(firingRates_exp);
% 
%     activeUnits = firingRates_exp > firingRateThr;
%     bottomFiringRate_exp = zeros(nUnits,nCentiles);
%     topFiringRate_exp = zeros(nUnits,nCentiles);
%     if ~isempty(pupilAreaSize)
%       percentileValues = prctile(pupilAreaSize,percentiles);
%       nCentiles = ceil(numel(percentiles)/2);
%       for iCent = 1:nCentiles
%         bottomPupilAreaSizeMask = pupilAreaSize <= percentileValues(iCent);
%         topPupilAreaSizeMask = pupilAreaSize > percentileValues(end - iCent + 1);
%         bottomFiringRate_exp(activeUnits,iCent) = sum(spikeCounts(activeUnits,bottomPupilAreaSizeMask),2)/expDuration;
%         topFiringRate_exp(activeUnits,iCent) = sum(spikeCounts(activeUnits,topPupilAreaSizeMask),2)/expDuration;
%       end
%       bottomFiringRate_exp(~bottomFiringRate_exp) = nan;
%       topFiringRate_exp(~topFiringRate_exp) = nan;
%       firingRateChange_exp = topFiringRate_exp - bottomFiringRate_exp;
%       logBottomFiringRate_exp = log10(bottomFiringRate_exp);
%       logTopFiringRate_exp = log10(topFiringRate_exp);
%       logFiringRateChange_exp = logTopFiringRate_exp - logBottomFiringRate_exp;
%     else
%       firingRateChange_exp = nan(nUnits,nCentiles);
%       logBottomFiringRate_exp = nan(nUnits,nCentiles);
%       logTopFiringRate_exp = nan(nUnits,nCentiles);
%       logFiringRateChange_exp = nan(nUnits,nCentiles);
%     end
%     firingRates_area(unitMask,:) = [rSpearman pvalSpearman ...
%       firingRates_exp logFiringRates_exp ...
%       bottomFiringRate_exp topFiringRate_exp firingRateChange_exp ...
%       logBottomFiringRate_exp logTopFiringRate_exp logFiringRateChange_exp];
%   end
%   firingRates{iArea} = firingRates_area;
% end
% 
% % Sort areas large to small positive fraction
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