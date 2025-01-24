% Use this script to correlate unit spiking with pupil area size.
% The script can only be executed on Matlab version R2021a or higher due to
% the use of functions with keyword-value argument pairs.

% Load parameters
params
pupilPassbandFrequency = 1.5;
pupilStopbandFrequency = 2;
averagedPupilDownsampling = true;
alpha = 0.05; % Significance level
excludeMovement = false;

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

% Correlate individual unit activity with the pupil area: organise by recording
nRecs = numel(infraslowData.experimentData);
rPearson = cell(nRecs,1);
pvalPearson = cell(nRecs,1);
rSpearman = cell(nRecs,1);
pvalSpearman = cell(nRecs,1);
for iRec = 1:nRecs
  disp(['Progress: ' num2str(100*iRec/nRecs) '%']);
  expData = infraslowData.experimentData{iRec};
  if isfield(expData, 'spikeCounts')
    nUnits = expData.nUnits;
    rPearson{iRec} = NaN(nUnits,1);
    pvalPearson{iRec} = NaN(nUnits,1);
    rSpearman{iRec} = NaN(nUnits,1);
    pvalSpearman{iRec} = NaN(nUnits,1);

    if ~isempty(expData.leftPupilAreaSize) || ~isempty(expData.rightPupilAreaSize)
      unitSpikeCounts = full(expData.spikeCounts);
      spikeTimeBins = expData.spikeTimeBins;
      if excludeMovement
        [spikeTimeBins, noMovementInds] = selectArrayValues( ...
          spikeTimeBins, expData.noMovementIntervals);
        unitSpikeCounts = unitSpikeCounts(noMovementInds);
      end
      [unitSpikeCounts, downsampledTimes] = resampleSpikeCounts( ...
        unitSpikeCounts, stepsize=1/expData.samplingRate, ...
        newStepsize=1/downsampledRate); % Downsample spiking data

      % Get the pupil area size
      if ~isempty(expData.leftPupilAreaSize) && ~isempty(expData.rightPupilAreaSize)
        pupilAreaSize = mean([expData.leftPupilAreaSize; expData.rightPupilAreaSize]);
      elseif ~isempty(expData.leftPupilAreaSize)
        pupilAreaSize = expData.leftPupilAreaSize;
      elseif ~isempty(expData.rightPupilAreaSize)
        pupilAreaSize = expData.rightPupilAreaSize;
      else
        continue
      end

      % Filter the pupil area size
      pupilAreaSizeFilt = pupilAreaSize;
      valueExists = ~isnan(pupilAreaSize);
      sr = 1/mean(diff(expData.spikeTimeBins));
      d = designfilt('lowpassiir', ...
        'PassbandFrequency',pupilPassbandFrequency, ...
        'StopbandFrequency',pupilStopbandFrequency, ...
        'PassbandRipple',0.5, 'StopbandAttenuation',1, ...
        'DesignMethod','butter', 'SampleRate',sr);
      pupilAreaSizeFilt(valueExists) = filtfilt(d,pupilAreaSize(valueExists));

      % Exclude movement periods
      if excludeMovement
        pupilAreaSizeFilt = pupilAreaSizeFilt(noMovementInds);
      end

      if averagedPupilDownsampling
        % Average the pupil area size (most accurate downsampling)
        averagedPupilAreaSize = movmean(pupilAreaSizeFilt, ...
          round(expData.samplingRate/downsampledRate), 'omitnan');
        downsampledPupilAreaSize = interp1(spikeTimeBins, averagedPupilAreaSize, ...
          downsampledTimes, 'linear', 'extrap');
      else
        % Downsample the pupil area size
        downsampledPupilAreaSize = interp1(spikeTimeBins, pupilAreaSizeFilt, ...
          downsampledTimes, 'linear', 'extrap'); %#ok<*UNRCH>
      end

      %fH = figure; plot(spikeTimeBins, pupilAreaSize, 'LineWidth',0.5); hold on
      %plot(spikeTimeBins, pupilAreaSizeFilt, 'LineWidth',0.5);
      %plot(downsampledTimes, downsampledPupilAreaSize, 'LineWidth',1.5);
      %plot(downsampledTimes, averagedPupilAreaSize, 'LineWidth',1.5); hold off
      %legend('smoothed','filtered','downsampled','averaged');
      %close(fH);

      if sum(isnan(downsampledPupilAreaSize)) == numel(downsampledPupilAreaSize) ...
          || ~any(downsampledPupilAreaSize)
        continue
      end

      % Correlate the two signals
      [rPearson{iRec}, pvalPearson{iRec}] = ...
        corrMulti(downsampledPupilAreaSize, unitSpikeCounts, 'Pearson');
      rPearson{iRec} = rPearson{iRec}'; pvalPearson{iRec} = pvalPearson{iRec}';
      [rSpearman{iRec}, pvalSpearman{iRec}] = ...
        corrMulti(downsampledPupilAreaSize, unitSpikeCounts, 'Spearman');
      rSpearman{iRec} = rSpearman{iRec}'; pvalSpearman{iRec} = pvalSpearman{iRec}';
    end
  end
end
spikingPupilCorr.recordings.rPearson = rPearson;
spikingPupilCorr.recordings.pvalPearson = pvalPearson;
spikingPupilCorr.recordings.rSpearman = rSpearman;
spikingPupilCorr.recordings.pvalSpearman = pvalSpearman;


% Correlate individual unit activity with the pupil area: organise by area
nAreas = numel(infraslowAnalyses.areaSummaries.groupedUnitInds);
rPearson = cell(nAreas,1);
pvalPearson = cell(nAreas,1);
rSpearman = cell(nAreas,1);
pvalSpearman = cell(nAreas,1);
for iArea = 1:nAreas
  disp(['Progress: ' num2str(100*iArea/nAreas) '%']);
  areaGroup = infraslowAnalyses.areaSummaries.areaTable.Brain_area_group{iArea};
  if ~strcmpi(areaGroup, '???')
    nUnits = size(infraslowAnalyses.areaSummaries.groupedUnitInds{iArea}, 1);
    rPearson{iArea} = NaN(nUnits,1);
    pvalPearson{iArea} = NaN(nUnits,1);
    rSpearman{iArea} = NaN(nUnits,1);
    pvalSpearman{iArea} = NaN(nUnits,1);
    prevExpInd = 0;
    for iUnit = 1:nUnits

      % Get the single unit spike count
      expInd = infraslowAnalyses.areaSummaries.groupedUnitInds{iArea}(iUnit,1);
      expData = infraslowData.experimentData{expInd};
      if ~isempty(expData.leftPupilAreaSize) || ~isempty(expData.rightPupilAreaSize)
        unitInd = infraslowAnalyses.areaSummaries.groupedUnitInds{iArea}(iUnit,2);
        unitSpikeCounts = full(expData.spikeCounts(unitInd,:));
        spikeTimeBins = expData.spikeTimeBins;
        if excludeMovement
          [spikeTimeBins, noMovementInds] = selectArrayValues( ...
            spikeTimeBins, expData.noMovementIntervals);
          unitSpikeCounts = unitSpikeCounts(noMovementInds);
        end
        [unitSpikeCounts, downsampledTimes] = resampleSpikeCounts( ...
          unitSpikeCounts, stepsize=1/expData.samplingRate, ...
          newStepsize=1/downsampledRate); % Downsample spiking data
        if size(unitSpikeCounts,2) == 1
          unitSpikeCounts = unitSpikeCounts';
        end

        % Get the pupil area size
        if expInd ~= prevExpInd
          if ~isempty(expData.leftPupilAreaSize) && ~isempty(expData.rightPupilAreaSize)
            pupilAreaSize = mean([expData.leftPupilAreaSize; expData.rightPupilAreaSize]);
          elseif ~isempty(expData.leftPupilAreaSize)
            pupilAreaSize = expData.leftPupilAreaSize;
          elseif ~isempty(expData.rightPupilAreaSize)
            pupilAreaSize = expData.rightPupilAreaSize;
          else
            break
          end

          % Filter the pupil area size
          pupilAreaSizeFilt = pupilAreaSize;
          valueExists = ~isnan(pupilAreaSize);
          sr = 1/mean(diff(expData.spikeTimeBins));
          d = designfilt('lowpassiir', ...
            'PassbandFrequency',pupilPassbandFrequency, ...
            'StopbandFrequency',pupilStopbandFrequency, ...
            'PassbandRipple',0.5, 'StopbandAttenuation',1, ...
            'DesignMethod','butter', 'SampleRate',sr);
          pupilAreaSizeFilt(valueExists) = filtfilt(d,pupilAreaSize(valueExists));

          % Exclude movement periods
          if excludeMovement
            pupilAreaSizeFilt = pupilAreaSizeFilt(noMovementInds);
          end

          if averagedPupilDownsampling
            % Average the pupil area size (most accurate downsampling)
            averagedPupilAreaSize = movmean(pupilAreaSizeFilt, ...
              round(expData.samplingRate/downsampledRate), 'omitnan');
            downsampledPupilAreaSize = interp1(spikeTimeBins, averagedPupilAreaSize, ...
              downsampledTimes, 'linear', 'extrap');
          else
            % Downsample the pupil area size
            downsampledPupilAreaSize = interp1(spikeTimeBins, pupilAreaSizeFilt, ...
              downsampledTimes, 'linear', 'extrap'); %#ok<*UNRCH>
          end

          %fH = figure; plot(spikeTimeBins, pupilAreaSize, 'LineWidth',0.5); hold on
          %plot(spikeTimeBins, pupilAreaSizeFilt, 'LineWidth',0.5);
          %plot(downsampledTimes, downsampledPupilAreaSize, 'LineWidth',1.5);
          %plot(downsampledTimes, averagedPupilAreaSize, 'LineWidth',1.5); hold off
          %legend('smoothed','filtered','downsampled','averaged');
          %close(fH);
        end
        if sum(isnan(downsampledPupilAreaSize)) == numel(downsampledPupilAreaSize) ...
            || ~any(downsampledPupilAreaSize)
          continue
        end

        % Correlate the two signals
        [rPearson{iArea}(iUnit), pvalPearson{iArea}(iUnit)] = ...
          corrMulti(downsampledPupilAreaSize, unitSpikeCounts, 'Pearson');
        [rSpearman{iArea}(iUnit), pvalSpearman{iArea}(iUnit)] = ...
          corrMulti(downsampledPupilAreaSize, unitSpikeCounts, 'Spearman');
      end

      prevExpInd = expInd;
    end
  end
end

% Work out correlated unit proportions in single areas
spikingPupilCorr.singleAreas = calcPupilCorrFractions( ...
  rPearson, pvalPearson, rSpearman, pvalSpearman, ...
  alpha, infraslowAnalyses.areaSummaries.groupedUnitInds);

% Group areas
areaGroups = unique(infraslowAnalyses.areaSummaries.areaTable.Brain_area_group);
nAreas = numel(areaGroups);
rPearsonGroups = cell(nAreas,1);
pvalPearsonGroups = cell(nAreas,1);
rSpearmanGroups = cell(nAreas,1);
pvalSpearmanGroups = cell(nAreas,1);
groupedUnitInds = cell(nAreas,1);
for iArea = 1:nAreas
  areaInds = ismember( ...
    infraslowAnalyses.areaSummaries.areaTable.Brain_area_group, areaGroups{iArea});
  rPearsonGroups{iArea} = concatenateCells(rPearson(areaInds));
  pvalPearsonGroups{iArea} = concatenateCells(pvalPearson(areaInds));
  rSpearmanGroups{iArea} = concatenateCells(rSpearman(areaInds));
  pvalSpearmanGroups{iArea} = concatenateCells(pvalSpearman(areaInds));
  groupedUnitInds{iArea} = concatenateCells( ...
    infraslowAnalyses.areaSummaries.groupedUnitInds(areaInds));
end
areas2keep = ~ismember(areaGroups, '???');
rPearsonGroups = rPearsonGroups(areas2keep);
pvalPearsonGroups = pvalPearsonGroups(areas2keep);
rSpearmanGroups = rSpearmanGroups(areas2keep);
pvalSpearmanGroups = pvalSpearmanGroups(areas2keep);
groupedUnitInds = groupedUnitInds(areas2keep);
areaGroups = areaGroups(areas2keep);

% Work out correlated unit proportions in area groups
spikingPupilCorr.areaGroups = calcPupilCorrFractions( ...
  rPearsonGroups, pvalPearsonGroups, rSpearmanGroups, pvalSpearmanGroups, ...
  alpha, groupedUnitInds);
spikingPupilCorr.areaGroups.areaAcronyms = areaGroups;

% Select areas of interest
nAreas = numel(areasOI);
rPearsonOI = cell(nAreas,1);
pvalPearsonOI = cell(nAreas,1);
rSpearmanOI = cell(nAreas,1);
pvalSpearmanOI = cell(nAreas,1);
groupedUnitInds = cell(nAreas,1);
for iArea = 1:nAreas
  areaInds = getAreaInds(areasOI{iArea}, infraslowAnalyses.areaSummaries.areaTable);
  rPearsonOI{iArea} = concatenateCells(rPearson(areaInds));
  pvalPearsonOI{iArea} = concatenateCells(pvalPearson(areaInds));
  rSpearmanOI{iArea} = concatenateCells(rSpearman(areaInds));
  pvalSpearmanOI{iArea} = concatenateCells(pvalSpearman(areaInds));
  groupedUnitInds{iArea} = concatenateCells( ...
    infraslowAnalyses.areaSummaries.groupedUnitInds(areaInds));
end

% Work out correlated unit proportions in area groups
spikingPupilCorr.areasOI = calcPupilCorrFractions( ...
  rPearsonOI, pvalPearsonOI, rSpearmanOI, pvalSpearmanOI, ...
  alpha, groupedUnitInds);
spikingPupilCorr.areasOI.areaAcronyms = areasOI;

% Sort areas large to small positive fraction
[sortedAreas, areaOrder] = sort( ...
  spikingPupilCorr.singleAreas.positiveSpearmanFractionsMeans(:,1), 'descend');
areaOrder = areaOrder(~isnan(sortedAreas));
spikingPupilCorr.singleAreas.fractionTable = table( ...
  infraslowAnalyses.areaSummaries.areaTable.Brain_area_name(areaOrder), ...
  sortedAreas(~isnan(sortedAreas)), ...
  'VariableNames', {'Brain_area', 'Positive_cell_fraction'});
disp(spikingPupilCorr.singleAreas.fractionTable);
%disp(infraslowAnalyses.spikingPupilCorr.singleAreas.fractionTable);

[sortedAreas, areaOrder] = sort( ...
  spikingPupilCorr.areasOI.positiveSpearmanFractionsMeans(:,1), 'descend');
areaOrder = areaOrder(~isnan(sortedAreas));
spikingPupilCorr.areasOI.fractionTable = table(areasOI(areaOrder), ...
  spikingPupilCorr.areasOI.positiveSpearmanFractionsMeans(areaOrder,1), ...
  'VariableNames', {'Brain_area', 'Positive_cell_fraction'});
disp(spikingPupilCorr.areasOI.fractionTable);
%disp(infraslowAnalyses.spikingPupilCorr.areasOI.fractionTable);

% Save the data
spikingPupilCorr.timeOfCompletion = datetime;
if excludeMovement
  infraslowAnalyses.spikingPupilCorr_noMovement = spikingPupilCorr;
else
  infraslowAnalyses.spikingPupilCorr = spikingPupilCorr;
end
save(analysisResultsFile, 'infraslowAnalyses', '-v7.3');



%% Local functions
function dataSummary = calcPupilCorrFractions(rPearson, pvalPearson, ...
  rSpearman, pvalSpearman, alpha, groupedUnitInds)

% Initialise data containers
nAreas = numel(rPearson);
positivePearsonFractions = NaN(nAreas,2);
negativePearsonFractions = NaN(nAreas,2);
neutralPearsonFractions = NaN(nAreas,2);
positiveSpearmanFractions = NaN(nAreas,2);
negativeSpearmanFractions = NaN(nAreas,2);
neutralSpearmanFractions = NaN(nAreas,2);
positivePearsonFractionsPerRec = cell(nAreas,1);
negativePearsonFractionsPerRec = cell(nAreas,1);
neutralPearsonFractionsPerRec = cell(nAreas,1);
positiveSpearmanFractionsPerRec = cell(nAreas,1);
negativeSpearmanFractionsPerRec = cell(nAreas,1);
neutralSpearmanFractionsPerRec = cell(nAreas,1);
positivePearsonFractionsMeans = NaN(nAreas,2);
negativePearsonFractionsMeans = NaN(nAreas,2);
neutralPearsonFractionsMeans = NaN(nAreas,2);
positiveSpearmanFractionsMeans = NaN(nAreas,2);
negativeSpearmanFractionsMeans = NaN(nAreas,2);
neutralSpearmanFractionsMeans = NaN(nAreas,2);
positivePearsonFractionsCI95 = NaN(nAreas,2);
negativePearsonFractionsCI95 = NaN(nAreas,2);
neutralPearsonFractionsCI95 = NaN(nAreas,2);
positiveSpearmanFractionsCI95 = NaN(nAreas,2);
negativeSpearmanFractionsCI95 = NaN(nAreas,2);
neutralSpearmanFractionsCI95 = NaN(nAreas,2);
for iArea = 1:nAreas
  disp(['Progress: ' num2str(100*iArea/nAreas) '%']);

  % Work out overall proportions
  if ~isempty(rPearson{iArea}) && any(~isnan(rPearson{iArea}))
    positivePearsonFractions(iArea,1) = ...
      sum(rPearson{iArea} > 0)/sum(~isnan(rPearson{iArea}));
    negativePearsonFractions(iArea,1) = ...
      sum(rPearson{iArea} < 0)/sum(~isnan(rPearson{iArea}));
    neutralPearsonFractions(iArea,1) = ...
      sum(~isnan(rPearson{iArea}))/numel(rPearson{iArea}) - ...
      positivePearsonFractions(iArea,1) - negativePearsonFractions(iArea,1);
    positivePearsonFractions(iArea,2) = ...
      sum(rPearson{iArea} > 0 & pvalPearson{iArea} < alpha)/sum(~isnan(rPearson{iArea}));
    negativePearsonFractions(iArea,2) = ...
      sum(rPearson{iArea} < 0 & pvalPearson{iArea} < alpha)/sum(~isnan(rPearson{iArea}));
    neutralPearsonFractions(iArea,2) = ...
      sum(~isnan(rPearson{iArea}))/numel(rPearson{iArea}) - ...
      positivePearsonFractions(iArea,2) - negativePearsonFractions(iArea,2);

    positiveSpearmanFractions(iArea,1) = ...
      sum(rSpearman{iArea} > 0)/sum(~isnan(rSpearman{iArea}));
    negativeSpearmanFractions(iArea,1) = ...
      sum(rSpearman{iArea} < 0)/sum(~isnan(rSpearman{iArea}));
    neutralSpearmanFractions(iArea,1) = ...
      sum(~isnan(rSpearman{iArea}))/numel(rSpearman{iArea}) - ...
      positiveSpearmanFractions(iArea,1) - negativeSpearmanFractions(iArea,1);
    positiveSpearmanFractions(iArea,2) = ...
      sum(rSpearman{iArea} > 0 & pvalSpearman{iArea} < alpha)/sum(~isnan(rSpearman{iArea}));
    negativeSpearmanFractions(iArea,2) = ...
      sum(rSpearman{iArea} < 0 & pvalSpearman{iArea} < alpha)/sum(~isnan(rSpearman{iArea}));
    neutralSpearmanFractions(iArea,2) = ...
      sum(~isnan(rSpearman{iArea}))/numel(rSpearman{iArea}) - ...
      positiveSpearmanFractions(iArea,2) - negativeSpearmanFractions(iArea,2);
  end

  % Work out proportions per recording
  if ~isempty(rPearson{iArea}) && any(~isnan(rPearson{iArea}))
    recs = unique(groupedUnitInds{iArea}(~isnan(rPearson{iArea}),1));
    nRecs = numel(recs);
    positivePearsonFractionsPerRec{iArea} = NaN(nRecs,2);
    negativePearsonFractionsPerRec{iArea} = NaN(nRecs,2);
    neutralPearsonFractionsPerRec{iArea} = NaN(nRecs,2);
    for iRec = 1:nRecs
      recInds = groupedUnitInds{iArea}(:,1) == recs(iRec) & ...
        ~isnan(rPearson{iArea});
      positivePearsonFractionsPerRec{iArea}(iRec,1) = ...
        sum(rPearson{iArea} > 0 & recInds)/sum(recInds);
      negativePearsonFractionsPerRec{iArea}(iRec,1) = ...
        sum(rPearson{iArea} < 0 & recInds)/sum(recInds);
      neutralPearsonFractionsPerRec{iArea}(iRec,1) = ...
        1 - positivePearsonFractionsPerRec{iArea}(iRec,1) - ...
        negativePearsonFractionsPerRec{iArea}(iRec,1);
      positivePearsonFractionsPerRec{iArea}(iRec,2) = ...
        sum(rPearson{iArea} > 0 & pvalPearson{iArea} < alpha & recInds)/sum(recInds);
      negativePearsonFractionsPerRec{iArea}(iRec,2) = ...
        sum(rPearson{iArea} < 0 & pvalPearson{iArea} < alpha & recInds)/sum(recInds);
      neutralPearsonFractionsPerRec{iArea}(iRec,2) = ...
        1 - positivePearsonFractionsPerRec{iArea}(iRec,2) - ...
        negativePearsonFractionsPerRec{iArea}(iRec,2);
    end
    if nRecs == 1
      positivePearsonFractionsMeans(iArea,:) = positivePearsonFractionsPerRec{iArea};
      positivePearsonFractionsCI95(iArea,:) = [0 0];
      negativePearsonFractionsMeans(iArea,:) = negativePearsonFractionsPerRec{iArea};
      negativePearsonFractionsCI95(iArea,:) = [0 0];
      neutralPearsonFractionsMeans(iArea,:) = neutralPearsonFractionsPerRec{iArea};
      neutralPearsonFractionsCI95(iArea,:) = [0 0];
    else
      [positivePearsonFractionsMeans(iArea,:), positivePearsonFractionsCI95Temp] = ...
        datamean(positivePearsonFractionsPerRec{iArea});
      positivePearsonFractionsCI95(iArea,:) = positivePearsonFractionsCI95Temp(2,:);
      [negativePearsonFractionsMeans(iArea,:), negativePearsonFractionsCI95Temp] = ...
        datamean(negativePearsonFractionsPerRec{iArea});
      negativePearsonFractionsCI95(iArea,:) = negativePearsonFractionsCI95Temp(2,:);
      [neutralPearsonFractionsMeans(iArea,:), neutralPearsonFractionsCI95Temp] = ...
        datamean(neutralPearsonFractionsPerRec{iArea});
      neutralPearsonFractionsCI95(iArea,:) = neutralPearsonFractionsCI95Temp(2,:);
    end

    positiveSpearmanFractionsPerRec{iArea} = NaN(nRecs,2);
    negativeSpearmanFractionsPerRec{iArea} = NaN(nRecs,2);
    neutralSpearmanFractionsPerRec{iArea} = NaN(nRecs,2);
    for iRec = 1:nRecs
      recInds = groupedUnitInds{iArea}(:,1) == recs(iRec);
      positiveSpearmanFractionsPerRec{iArea}(iRec,1) = ...
        sum(rSpearman{iArea} > 0 & recInds)/sum(recInds);
      negativeSpearmanFractionsPerRec{iArea}(iRec,1) = ...
        sum(rSpearman{iArea} < 0 & recInds)/sum(recInds);
      neutralSpearmanFractionsPerRec{iArea}(iRec,1) = ...
        1 - positiveSpearmanFractionsPerRec{iArea}(iRec,1) - ...
        negativeSpearmanFractionsPerRec{iArea}(iRec,1);
      positiveSpearmanFractionsPerRec{iArea}(iRec,2) = ...
        sum(rSpearman{iArea} > 0 & pvalSpearman{iArea} < alpha & recInds)/sum(recInds);
      negativeSpearmanFractionsPerRec{iArea}(iRec,2) = ...
        sum(rSpearman{iArea} < 0 & pvalSpearman{iArea} < alpha & recInds)/sum(recInds);
      neutralSpearmanFractionsPerRec{iArea}(iRec,2) = ...
        1 - positiveSpearmanFractionsPerRec{iArea}(iRec,2) - ...
        negativeSpearmanFractionsPerRec{iArea}(iRec,2);
    end
    if nRecs == 1
      positiveSpearmanFractionsMeans(iArea,:) = positiveSpearmanFractionsPerRec{iArea};
      positiveSpearmanFractionsCI95(iArea,:) = [0 0];
      negativeSpearmanFractionsMeans(iArea,:) = negativeSpearmanFractionsPerRec{iArea};
      negativeSpearmanFractionsCI95(iArea,:) = [0 0];
      neutralSpearmanFractionsMeans(iArea,:) = neutralSpearmanFractionsPerRec{iArea};
      neutralSpearmanFractionsCI95(iArea,:) = [0 0];
    else
      [positiveSpearmanFractionsMeans(iArea,:), positiveSpearmanFractionsCI95Temp] = ...
        datamean(positiveSpearmanFractionsPerRec{iArea});
      positiveSpearmanFractionsCI95(iArea,:) = positiveSpearmanFractionsCI95Temp(2,:);
      [negativeSpearmanFractionsMeans(iArea,:), negativeSpearmanFractionsCI95Temp] = ...
        datamean(negativeSpearmanFractionsPerRec{iArea});
      negativeSpearmanFractionsCI95(iArea,:) = negativeSpearmanFractionsCI95Temp(2,:);
      [neutralSpearmanFractionsMeans(iArea,:), neutralSpearmanFractionsCI95Temp] = ...
        datamean(neutralSpearmanFractionsPerRec{iArea});
      neutralSpearmanFractionsCI95(iArea,:) = neutralSpearmanFractionsCI95Temp(2,:);
    end
  end
end

% Store single area data
dataSummary.rPearson = rPearson;
dataSummary.pvalPearson = pvalPearson;
dataSummary.rSpearman = rSpearman;
dataSummary.pvalSpearman = pvalSpearman;

dataSummary.positivePearsonFractions = positivePearsonFractions;
dataSummary.negativePearsonFractions = negativePearsonFractions;
dataSummary.neutralPearsonFractions = neutralPearsonFractions;
dataSummary.positiveSpearmanFractions = positiveSpearmanFractions;
dataSummary.negativeSpearmanFractions = negativeSpearmanFractions;
dataSummary.neutralSpearmanFractions = neutralSpearmanFractions;

dataSummary.positivePearsonFractionsPerRec = positivePearsonFractionsPerRec;
dataSummary.negativePearsonFractionsPerRec = negativePearsonFractionsPerRec;
dataSummary.neutralPearsonFractionsPerRec = neutralPearsonFractionsPerRec;
dataSummary.positiveSpearmanFractionsPerRec = positiveSpearmanFractionsPerRec;
dataSummary.negativeSpearmanFractionsPerRec = negativeSpearmanFractionsPerRec;
dataSummary.neutralSpearmanFractionsPerRec = neutralSpearmanFractionsPerRec;

dataSummary.positivePearsonFractionsMeans = positivePearsonFractionsMeans;
dataSummary.negativePearsonFractionsMeans = negativePearsonFractionsMeans;
dataSummary.neutralPearsonFractionsMeans = neutralPearsonFractionsMeans;
dataSummary.positiveSpearmanFractionsMeans = positiveSpearmanFractionsMeans;
dataSummary.negativeSpearmanFractionsMeans = negativeSpearmanFractionsMeans;
dataSummary.neutralSpearmanFractionsMeans = neutralSpearmanFractionsMeans;

dataSummary.positivePearsonFractionsCI95 = positivePearsonFractionsCI95;
dataSummary.negativePearsonFractionsCI95 = negativePearsonFractionsCI95;
dataSummary.neutralPearsonFractionsCI95 = neutralPearsonFractionsCI95;
dataSummary.positiveSpearmanFractionsCI95 = positiveSpearmanFractionsCI95;
dataSummary.negativeSpearmanFractionsCI95 = negativeSpearmanFractionsCI95;
dataSummary.neutralSpearmanFractionsCI95 = neutralSpearmanFractionsCI95;
end