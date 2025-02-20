% Coherence analysis script comparing unit spiking activity with respect to
% the pupil area size. The script can only be executed on Matlab version
% R2020a or higher due to the use of arguments section in functions
% (adapted to run on a computing cluster).

% Load parameters
addDependencies
params
parallelCores = 1; %32;
pupilPassbandFrequency = 1.5;
pupilStopbandFrequency = 2;
population = 'all'; %'all', 'positive', 'negative'
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
    areaAcronymInit = infraslowAnalyses.areaSummaries.areaTable.Brain_area_acronym{iArea};
    areaAcronym = strrep(strrep(areaAcronymInit, ' ', '_'), '-', '_');
    recs = unique(areaInds(:,1));
    for iRec = 1:numel(recs)
      recID = recs(iRec);
      expData = infraslowData.experimentData{recID};

      % Get all unit spike counts
      unitInds = ismember(expData.unitBrainAreas, areaAcronymInit);
      if strcmpi(population, 'positive') || strcmpi(population, 'negative')
        expPupilCorrData = infraslowAnalyses.spikingPupilCorr.recordings.rPearson{recID};
        if ~isempty(expPupilCorrData) && any(~isnan(expPupilCorrData))
          if strcmpi(population, 'positive')
            unitInds = unitInds & expPupilCorrData > 0;
          elseif strcmpi(population, 'negative')
            unitInds = unitInds & expPupilCorrData < 0;
          end
        else
          continue
        end
      end
      if any(unitInds)
        spikeCounts = expData.spikeCounts(unitInds,:);
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
        %if exist('dependenciesAdded', 'var') && dependenciesAdded
        [spikingPupilCoh.(areaAcronym){iRec}.fullCoherence, ...
          spikingPupilCoh.(areaAcronym){iRec}.half1Coherence, ...
          spikingPupilCoh.(areaAcronym){iRec}.half2Coherence, ...
          spikingPupilCoh.(areaAcronym){iRec}.fullInterpCoherence, ...
          spikingPupilCoh.(areaAcronym){iRec}.half1InterpCoherence, ...
          spikingPupilCoh.(areaAcronym){iRec}.half2InterpCoherence] = ...
          coherence(spikeCountsCell, pupilAreaSizeFilt, [], ...
          1/effectiveSR, 1/effectiveSR, [0 0], FOI, ...
          'pbc', 'c', winfactor, ...
          freqfactor, tapers, false, true, false, 0, true, true, true, ...
          parallelise);
        %else
        %  [spikingPupilCoh.(areaAcronym){iRec}.fullCoherence, ...
        %    spikingPupilCoh.(areaAcronym){iRec}.half1Coherence, ...
        %    spikingPupilCoh.(areaAcronym){iRec}.half2Coherence, ...
        %    spikingPupilCoh.(areaAcronym){iRec}.fullInterpCoherence, ...
        %    spikingPupilCoh.(areaAcronym){iRec}.half1InterpCoherence, ...
        %    spikingPupilCoh.(areaAcronym){iRec}.half2InterpCoherence] = ...
        %    coherence(spikeCountsCell, pupilAreaSizeFilt, ...
        %    stepsize=1/effectiveSR, startTime=1/effectiveSR, freqGrid=FOI, ...
        %    typespk1='pbc', typespk2='c', winfactor=winfactor, ...
        %    freqfactor=freqfactor, tapers=tapers, halfCoherence=true, ...
        %    parallelise=parallelise);
        %end
        spikingPupilCoh.(areaAcronym){iRec}.experimentIndex = recID;
        spikingPupilCoh.(areaAcronym){iRec}.timeOfCompletion = datetime;
      end
    end

    % Save data analysis results
    if ~isfield(spikingPupilCoh, areaAcronym)
      spikingPupilCoh.(areaAcronym) = {};
    end
    if strcmpi(population, 'all')
      if excludeMovement
        infraslowAnalyses.spikingPupilCoh_noMovement.(areaAcronym) = spikingPupilCoh.(areaAcronym);
      else
        infraslowAnalyses.spikingPupilCoh.(areaAcronym) = spikingPupilCoh.(areaAcronym); %#ok<*UNRCH>
      end
    elseif strcmpi(population, 'positive')
      if excludeMovement
        infraslowAnalyses.spikingPupilCohPositive_noMovement.(areaAcronym) = spikingPupilCoh.(areaAcronym);
      else
        infraslowAnalyses.spikingPupilCohPositive.(areaAcronym) = spikingPupilCoh.(areaAcronym); %#ok<*UNRCH>
      end
    elseif strcmpi(population, 'negative')
      if excludeMovement
        infraslowAnalyses.spikingPupilCohNegative_noMovement.(areaAcronym) = spikingPupilCoh.(areaAcronym);
      else
        infraslowAnalyses.spikingPupilCohNegative.(areaAcronym) = spikingPupilCoh.(areaAcronym); %#ok<*UNRCH>
      end
    end
    save(analysisResultsFile, 'infraslowAnalyses', '-v7.3');
  end
end

% Switch on warnings
if parallelCores > 1
  parfevalOnAll(@warning,0,'on','all');
end
warning('on', 'all');