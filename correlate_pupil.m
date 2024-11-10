% Correlate neural activity with pupil area

% Load parameters
params
pupilPassbandFrequency = 0.2; % 1.5
pupilStopbandFrequency = 0.25; % 2
averagedPupilDownsampling = true;
alpha = 0.05; % Significance level

% Load preprocessed data
preprocessedDataFile = fullfile(processedDataFolder, 'bwmPreprocessedData2.mat');
load(preprocessedDataFile);

% Load data analysis results
analysisResultsFile = fullfile(processedDataFolder, 'bwmAnalysisResults.mat');
load(analysisResultsFile);

% Correlate individual unit activity with the pupil area
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
        [unitSpikeCounts, downsampledTimes] = resampleSpikeCounts( ...
          unitSpikeCounts, stepsize=1/expData.samplingRate, ...
          newStepsize=1/downsampledRate); % Downsample spiking data
        downsampledTimes = downsampledTimes - 0.5/downsampledRate;

        % Get the pupil area size
        if expInd ~= prevExpInd
          if ~isempty(expData.leftPupilAreaSize) && ~isempty(expData.rightPupilAreaSize)
            pupilAreaSize = mean([expData.leftPupilAreaSize; expData.rightPupilAreaSize]);
          elseif ~isempty(expData.leftPupilAreaSize)
            pupilAreaSize = expData.leftPupilAreaSize;
          elseif ~isempty(expData.rightPupilAreaSize)
            pupilAreaSize = expData.rightPupilAreaSize;
          end

          % Filter pupil area size
          pupilAreaSizeFilt = pupilAreaSize;
          valueExists = ~isnan(pupilAreaSize);
          sr = 1/mean(diff(spikeTimeBins));
          d = designfilt('lowpassiir', ...
            'PassbandFrequency',pupilPassbandFrequency, ...
            'StopbandFrequency',pupilStopbandFrequency, ...
            'PassbandRipple',0.5, 'StopbandAttenuation',1, ...
            'DesignMethod','butter', 'SampleRate',sr);
          pupilAreaSizeFilt(valueExists) = filtfilt(d,pupilAreaSize(valueExists));

          if averagedPupilDownsampling
            % Average pupil area size (most accurate downsampling)
            averagedPupilAreaSize = movmean(pupilAreaSizeFilt, ...
              round(expData.samplingRate/downsampledRate), 'omitnan');
            downsampledPupilAreaSize = interp1(spikeTimeBins, averagedPupilAreaSize, ...
              downsampledTimes, 'linear', 'extrap');
          else
            % Downsample pupil area size
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
        if sum(isnan(downsampledPupilAreaSize)) == numel(downsampledPupilAreaSize)
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

% Work out correlated unit proportions
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
  if ~isempty(rPearson{iArea})
    positivePearsonFractions(iArea,1) = ...
      sum(rPearson{iArea} > 0)/numel(rPearson{iArea});
    negativePearsonFractions(iArea,1) = ...
      sum(rPearson{iArea} < 0)/numel(rPearson{iArea});
    neutralPearsonFractions(iArea,1) = ...
      1 - positivePearsonFractions(iArea,1) - negativePearsonFractions(iArea,1);
    positivePearsonFractions(iArea,2) = ...
      sum(rPearson{iArea} > 0 & pvalPearson{iArea} < alpha)/numel(rPearson{iArea});
    negativePearsonFractions(iArea,2) = ...
      sum(rPearson{iArea} < 0 & pvalPearson{iArea} < alpha)/numel(rPearson{iArea});
    neutralPearsonFractions(iArea,2) = ...
      1 - positivePearsonFractions(iArea,2) - negativePearsonFractions(iArea,2);

    positiveSpearmanFractions(iArea,1) = ...
      sum(rSpearman{iArea} > 0)/numel(rSpearman{iArea});
    negativeSpearmanFractions(iArea,1) = ...
      sum(rSpearman{iArea} < 0)/numel(rSpearman{iArea});
    neutralSpearmanFractions(iArea,1) = ...
      1 - positiveSpearmanFractions(iArea,1) - negativeSpearmanFractions(iArea,1);
    positiveSpearmanFractions(iArea,2) = ...
      sum(rSpearman{iArea} > 0 & pvalSpearman{iArea} < alpha)/numel(rSpearman{iArea});
    negativeSpearmanFractions(iArea,2) = ...
      sum(rSpearman{iArea} < 0 & pvalSpearman{iArea} < alpha)/numel(rSpearman{iArea});
    neutralSpearmanFractions(iArea,2) = ...
      1 - positiveSpearmanFractions(iArea,2) - negativeSpearmanFractions(iArea,2);
  end

  % Work out proportions per recording
  if ~isempty(rPearson{iArea})
    recs = unique(infraslowAnalyses.areaSummaries.groupedUnitInds{iArea}(:,1));
    nRecs = numel(recs);
    positivePearsonFractionsPerRec{iArea} = NaN(nRecs,2);
    negativePearsonFractionsPerRec{iArea} = NaN(nRecs,2);
    neutralPearsonFractionsPerRec{iArea} = NaN(nRecs,2);
    for iRec = 1:nRecs
      recInds = infraslowAnalyses.areaSummaries.groupedUnitInds{iArea}(:,1) == recs(iRec);
      positivePearsonFractionsPerRec{iArea}(iRec,1) = ...
        sum(rPearson{iArea} > 0 & recInds)/numel(rPearson{iArea});
      negativePearsonFractionsPerRec{iArea}(iRec,1) = ...
        sum(rPearson{iArea} < 0 & recInds)/numel(rPearson{iArea});
      neutralPearsonFractionsPerRec{iArea}(iRec,1) = ...
        1 - positivePearsonFractionsPerRec{iArea}(iRec,1) - negativePearsonFractionsPerRec{iArea}(iRec,1);
      positivePearsonFractionsPerRec{iArea}(iRec,2) = ...
        sum(rPearson{iArea} > 0 & pvalPearson{iArea} < alpha & recInds)/numel(rPearson{iArea});
      negativePearsonFractionsPerRec{iArea}(iRec,2) = ...
        sum(rPearson{iArea} < 0 & pvalPearson{iArea} < alpha & recInds)/numel(rPearson{iArea});
      neutralPearsonFractionsPerRec{iArea}(iRec,2) = ...
        1 - positivePearsonFractionsPerRec{iArea}(iRec,2) - negativePearsonFractionsPerRec{iArea}(iRec,2);
    end
    if nRecs == 1
      positivePearsonFractionsMeans(iArea,:) = positivePearsonFractionsPerRec{iArea};
      positivePearsonFractionsCI95Temp = [0 0];
      negativePearsonFractionsMeans(iArea,:) = negativePearsonFractionsPerRec{iArea};
      negativePearsonFractionsCI95Temp = [0 0];
      neutralPearsonFractionsMeans(iArea,:) = neutralPearsonFractionsPerRec{iArea};
      neutralPearsonFractionsCI95Temp = [0 0];
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
      recInds = infraslowAnalyses.areaSummaries.groupedUnitInds{iArea}(:,1) == recs(iRec);
      positiveSpearmanFractionsPerRec{iArea}(iRec,1) = ...
        sum(rSpearman{iArea} > 0 & recInds)/numel(rSpearman{iArea});
      negativeSpearmanFractionsPerRec{iArea}(iRec,1) = ...
        sum(rSpearman{iArea} < 0 & recInds)/numel(rSpearman{iArea});
      neutralSpearmanFractionsPerRec{iArea}(iRec,1) = ...
        1 - positiveSpearmanFractionsPerRec{iArea}(iRec,1) - negativeSpearmanFractionsPerRec{iArea}(iRec,1);
      positiveSpearmanFractionsPerRec{iArea}(iRec,2) = ...
        sum(rSpearman{iArea} > 0 & pvalSpearman{iArea} < alpha & recInds)/numel(rSpearman{iArea});
      negativeSpearmanFractionsPerRec{iArea}(iRec,2) = ...
        sum(rSpearman{iArea} < 0 & pvalSpearman{iArea} < alpha & recInds)/numel(rSpearman{iArea});
      neutralSpearmanFractionsPerRec{iArea}(iRec,2) = ...
        1 - positiveSpearmanFractionsPerRec{iArea}(iRec,2) - negativeSpearmanFractionsPerRec{iArea}(iRec,2);
    end
    if nRecs == 1
      positiveSpearmanFractionsMeans(iArea,:) = positiveSpearmanFractionsPerRec{iArea};
      positiveSpearmanFractionsCI95Temp = [0 0];
      negativeSpearmanFractionsMeans(iArea,:) = negativeSpearmanFractionsPerRec{iArea};
      negativeSpearmanFractionsCI95Temp = [0 0];
      neutralSpearmanFractionsMeans(iArea,:) = neutralSpearmanFractionsPerRec{iArea};
      neutralSpearmanFractionsCI95Temp = [0 0];
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

% Save the data
infraslowAnalyses.spikingPupilCorr.rPearson = rPearson;
infraslowAnalyses.spikingPupilCorr.pvalPearson = pvalPearson;
infraslowAnalyses.spikingPupilCorr.rSpearman = rSpearman;
infraslowAnalyses.spikingPupilCorr.pvalSpearman = pvalSpearman;

infraslowAnalyses.spikingPupilCorr.positivePearsonFractions = positivePearsonFractions;
infraslowAnalyses.spikingPupilCorr.negativePearsonFractions = negativePearsonFractions;
infraslowAnalyses.spikingPupilCorr.neutralPearsonFractions = neutralPearsonFractions;
infraslowAnalyses.spikingPupilCorr.positiveSpearmanFractions = positiveSpearmanFractions;
infraslowAnalyses.spikingPupilCorr.negativeSpearmanFractions = negativeSpearmanFractions;
infraslowAnalyses.spikingPupilCorr.neutralSpearmanFractions = neutralSpearmanFractions;

infraslowAnalyses.spikingPupilCorr.positivePearsonFractionsPerRec = positivePearsonFractionsPerRec;
infraslowAnalyses.spikingPupilCorr.negativePearsonFractionsPerRec = negativePearsonFractionsPerRec;
infraslowAnalyses.spikingPupilCorr.neutralPearsonFractionsPerRec = neutralPearsonFractionsPerRec;
infraslowAnalyses.spikingPupilCorr.positiveSpearmanFractionsPerRec = positiveSpearmanFractionsPerRec;
infraslowAnalyses.spikingPupilCorr.negativeSpearmanFractionsPerRec = negativeSpearmanFractionsPerRec;
infraslowAnalyses.spikingPupilCorr.neutralSpearmanFractionsPerRec = neutralSpearmanFractionsPerRec;

infraslowAnalyses.spikingPupilCorr.positivePearsonFractionsMeans = positivePearsonFractionsMeans;
infraslowAnalyses.spikingPupilCorr.negativePearsonFractionsMeans = negativePearsonFractionsMeans;
infraslowAnalyses.spikingPupilCorr.neutralPearsonFractionsMeans = neutralPearsonFractionsMeans;
infraslowAnalyses.spikingPupilCorr.positiveSpearmanFractionsMeans = positiveSpearmanFractionsMeans;
infraslowAnalyses.spikingPupilCorr.negativeSpearmanFractionsMeans = negativeSpearmanFractionsMeans;
infraslowAnalyses.spikingPupilCorr.neutralSpearmanFractionsMeans = neutralSpearmanFractionsMeans;

infraslowAnalyses.spikingPupilCorr.positivePearsonFractionsCI95 = positivePearsonFractionsCI95;
infraslowAnalyses.spikingPupilCorr.negativePearsonFractionsCI95 = negativePearsonFractionsCI95;
infraslowAnalyses.spikingPupilCorr.neutralPearsonFractionsCI95 = neutralPearsonFractionsCI95;
infraslowAnalyses.spikingPupilCorr.positiveSpearmanFractionsCI95 = positiveSpearmanFractionsCI95;
infraslowAnalyses.spikingPupilCorr.negativeSpearmanFractionsCI95 = negativeSpearmanFractionsCI95;
infraslowAnalyses.spikingPupilCorr.neutralSpearmanFractionsCI95 = neutralSpearmanFractionsCI95;

save(analysisResultsFile, 'infraslowAnalyses', '-v7.3');