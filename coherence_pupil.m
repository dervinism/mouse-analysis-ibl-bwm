% Coherence analysis script comparing unit spiking activity with respect to the pupil area  size

% Load parameters
addDependencies
params
parallelCores = 32;
pupilPassbandFrequency = 1.5;
pupilStopbandFrequency = 2;
excludeMovement = false;
resumeInd = 1;

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

% Set up parallelisation
warning('off', 'all');
if parallelCores > 1
  parallelise = true;
  p = gcp('nocreate');
  if isempty(p)
    parpool(parallelCores);
  end
  parfevalOnAll(@warning,0,'off','all');
else
  parallelise = false;
end

% Carry out coherence analysis of individual units wrt the pupil area
warning('off', 'all');
nAreas = numel(infraslowAnalyses.areaSummaries.groupedUnitInds);
spikingPupilCoh = struct();
for iArea = resumeInd:nAreas
  disp(['Progress: ' num2str(100*iArea/nAreas) '%']);
  areaGroup = infraslowAnalyses.areaSummaries.areaTable.Brain_area_group{iArea};
  areaInds = infraslowAnalyses.areaSummaries.groupedUnitInds{iArea};
  if ~strcmpi(areaGroup, '???')
    areaAcronym = strrep(strrep( ...
      infraslowAnalyses.areaSummaries.areaTable.Brain_area_acronym{iArea}, ' ', '_'), '-', '_');
    recs = unique(areaInds(:,1));
    for iRec = 1:numel(recs)
      recID = recs(iRec);
      expData = infraslowData.experimentData{recID};

      % Get all unit spike counts
      spikeCounts = expData.spikeCounts;
      spikeTimeBins = expData.spikeTimeBins;
      effectiveSR = 1/mean(diff(spikeTimeBins));
      if excludeMovement
        [spikeTimeBins, noMovementInds] = selectArrayValues( ...
          spikeTimeBins, expData.noMovementIntervals);
        spikeCounts = spikeCounts(:,noMovementInds);
      end

      % Get the pupil are size
      if ~isempty(expData.leftPupilAreaSize) && ~isempty(expData.rightPupilAreaSize)
        pupilAreaSize = mean([expData.leftPupilAreaSize; expData.rightPupilAreaSize]);
      elseif ~isempty(expData.leftPupilAreaSize)
        pupilAreaSize = expData.leftPupilAreaSize;
      elseif ~isempty(expData.rightPupilAreaSize)
        pupilAreaSize = expData.rightPupilAreaSize;
      else
        spikingPupilCoh.(areaAcronym){iRec}.experimentIndex = recID;
        spikingPupilCoh.(areaAcronym){iRec}.timeOfCompletion = datetime;
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

      % Exclude periods with no values
      nanInds = isnan(pupilAreaSizeFilt);
      if sum(nanInds) == numel(pupilAreaSizeFilt)
        spikingPupilCoh.(areaAcronym){iRec}.experimentIndex = recID;
        spikingPupilCoh.(areaAcronym){iRec}.timeOfCompletion = datetime;
        continue
      else
        pupilAreaSizeFilt = pupilAreaSizeFilt(~nanInds); 
      end

      % Put spike counts into a cell array for coherence analysis
      nUnits = size(spikeCounts,1);
      spikeCountsCell = cell(nUnits,1);
      for iUnit = 1:nUnits
        spikeCountsCell{iUnit} = full(spikeCounts(iUnit,~nanInds));
      end

      % Coherence analysis
      [spikingPupilCoh.(areaAcronym){iRec}.fullCoherence, ...
        spikingPupilCoh.(areaAcronym){iRec}.half1Coherence, ...
        spikingPupilCoh.(areaAcronym){iRec}.half2Coherence, ...
        spikingPupilCoh.(areaAcronym){iRec}.fullInterpCoherence, ...
        spikingPupilCoh.(areaAcronym){iRec}.half1InterpCoherence, ...
        spikingPupilCoh.(areaAcronym){iRec}.half2InterpCoherence] = ...
        coherence(spikeCountsCell, pupilAreaSizeFilt, ...
        stepsize=1/effectiveSR, startTime=1/effectiveSR, freqGrid=FOI, ...
        typespk1='pbc', typespk2='c', winfactor=winfactor, ...
        freqfactor=freqfactor, tapers=tapers, halfCoherence=true, ...
        parallelise=parallelise);
      spikingPupilCoh.(areaAcronym){iRec}.experimentIndex = recID;
      spikingPupilCoh.(areaAcronym){iRec}.timeOfCompletion = datetime;
    end

    % Save data analysis results
    if excludeMovement
      infraslowAnalyses.spikingPupilCoh_noMovement.(areaAcronym) = spikingPupilCoh.(areaAcronym);
    else
      infraslowAnalyses.spikingPupilCoh.(areaAcronym) = spikingPupilCoh.(areaAcronym); %#ok<*UNRCH>
    end
    save(analysisResultsFile, 'infraslowAnalyses', '-v7.3');
  end
end

% Switch on warnings
if parallelCores > 1
  parfevalOnAll(@warning,0,'on','all');
end
warning('on', 'all');