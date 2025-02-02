% Coherence analysis script comparing unit spiking activity with respect to
% the population rates of different brain areas. The script can only be
% executed on Matlab version R2020a or higher due to the use of arguments
% section in functions (adapted to run on a computing cluster).

% Load parameters
addDependencies
params
parallelCores = 32;
excludeMovement = false;
population = 'all'; %'all', 'positive', 'negative'
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

% Go through all areas of interest, then through every experiment selecting
% units, and then through areas again selecting populations
nRecs = numel(infraslowData.experimentData);

% Proceed area by area
areasOI = areasMinimal;
nAreas = numel(areasOI);
spikingSpikingCoh = struct();
for iUnitArea = resumeInd:nAreas
  disp(['Progress: ' num2str(100*iUnitArea/nAreas) '% (' num2str(iUnitArea) ')']);
  unitArea = areasOI{iUnitArea};
  areaInds = getAreaInds(unitArea, infraslowAnalyses.areaSummaries.areaTable);
  subAreasOI = infraslowAnalyses.areaSummaries.areaTable.Brain_area_acronym(areaInds);

  % Proceed experiment by experiment
  for iRec = 1:nRecs
    expData = infraslowData.experimentData{iRec};
    if isfield(expData, 'unitBrainAreas')
      unitsOI = ismember(expData.unitBrainAreas, subAreasOI);
      if strcmpi(population, 'positive') || strcmpi(population, 'negative')
        expPupilCorrData = infraslowAnalyses.spikingPupilCorr.recordings.rPearson{iRec};
        if ~isempty(expPupilCorrData) && any(~isnan(expPupilCorrData))
          if strcmpi(population, 'positive')
            unitsOI = unitsOI & expPupilCorrData > 0;
          elseif strcmpi(population, 'negative')
            unitsOI = unitsOI & expPupilCorrData < 0;
          end
        else
          continue
        end
      end
      if any(unitsOI)

        % Get all unit spike counts
        unitSpikeCounts = expData.spikeCounts(unitsOI,:);
        spikeTimeBins = expData.spikeTimeBins;
        effectiveSR = 1/mean(diff(spikeTimeBins));
        if excludeMovement
          [spikeTimeBins, noMovementInds] = selectArrayValues( ...
            spikeTimeBins, expData.noMovementIntervals); %#ok<*UNRCH>
          unitSpikeCounts = unitSpikeCounts(:,noMovementInds);
        end

        % Put spike counts into a cell array for coherence analysis
        nUnits = size(unitSpikeCounts,1);
        unitSpikeCountsCell = cell(nUnits,1);
        for iUnit = 1:nUnits
          unitSpikeCountsCell{iUnit} = full(unitSpikeCounts(iUnit,:));
        end

        % Proceed area by area again
        for iPopulationArea = 1:nAreas
          populationArea = areasOI{iPopulationArea};
          areaInds = getAreaInds(populationArea, infraslowAnalyses.areaSummaries.areaTable);
          subAreasOI_population = infraslowAnalyses.areaSummaries.areaTable.Brain_area_acronym(areaInds);
          unitArea = strrep(strrep(unitArea, ' ', '_'), '-', '_');
          populationArea = strrep(strrep(populationArea, ' ', '_'), '-', '_');
          if strcmpi(unitArea, populationArea)
            continue
          end

          % Get the population rate
          unitInds = ismember(expData.unitBrainAreas, subAreasOI_population);
          if strcmpi(population, 'positive')
            unitInds = unitInds & expPupilCorrData > 0;
          elseif strcmpi(population, 'negative')
            unitInds = unitInds & expPupilCorrData < 0;
          end
          if ~any(unitInds)
            continue
          end
          populationRate = full(sum(expData.spikeCounts(unitInds,:)));
          if excludeMovement
            populationRate = populationRate(noMovementInds);
          end

          % Coherence analysis
          if ~isfield(spikingSpikingCoh, unitArea) || ...
              ~isfield(spikingSpikingCoh.(unitArea), populationArea)
            spikingSpikingCoh.(unitArea).(populationArea) = {};
          end
          recInd = numel(spikingSpikingCoh.(unitArea).(populationArea))+1;
          if sum(populationRate)
            %if exist('dependenciesAdded', 'var') && dependenciesAdded
              [~, ~, ~, ...
                spikingSpikingCoh.(unitArea).(populationArea){recInd}.fullInterpCoherence, ...
                spikingSpikingCoh.(unitArea).(populationArea){recInd}.half1InterpCoherence, ...
                spikingSpikingCoh.(unitArea).(populationArea){recInd}.half2InterpCoherence] = ...
                coherence(unitSpikeCountsCell, populationRate, [], ...
                1/effectiveSR, 1/effectiveSR, [0 0], FOI, ...
                'pbc', 'pbc', winfactor, ...
                freqfactor, tapers, false, true, false, 0, true, true, true, ...
                parallelise);
            %else
            %  [~, ~, ~, ...
            %    spikingSpikingCoh.(unitArea).(populationArea){recInd}.fullInterpCoherence, ...
            %    spikingSpikingCoh.(unitArea).(populationArea){recInd}.half1InterpCoherence, ...
            %    spikingSpikingCoh.(unitArea).(populationArea){recInd}.half2InterpCoherence] = ...
            %    coherence(unitSpikeCountsCell, populationRate, ...
            %    stepsize=1/effectiveSR, startTime=1/effectiveSR, freqGrid=FOI, ...
            %    typespk1='pbc', typespk2='pbc', winfactor=winfactor, ...
            %    freqfactor=freqfactor, tapers=tapers, halfCoherence=true, ...
            %    parallelise=parallelise);
            %end

            % Record metadata
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.units = expData.units(unitsOI);
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.nUnits = sum(unitsOI);
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.unitChannels = expData.unitChannels(unitsOI);
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.horzCoords = expData.horzCoords(unitsOI);
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.vertCoords = expData.vertCoords(unitsOI);
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.unitBrainAreas = expData.unitBrainAreas(unitsOI);
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.goodUnits = expData.goodUnits(unitsOI);
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.unitsPR = expData.units(unitInds);
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.populationBrainAreas = expData.unitBrainAreas(unitInds);
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.recordingNumber = iRec;
            spikingSpikingCoh.(unitArea).(populationArea){recInd}.timeOfCompletion = datetime;

            % Visualise data
            %fH = figure;
            %hold on
            %for iUnit = 1:spikingSpikingCoh.(unitArea).(populationArea){recInd}.nUnits
            %  plot(spikingSpikingCoh.(unitArea).(populationArea){recInd}.fullInterpCoherence.frequency(iUnit,:), ...
            %    spikingSpikingCoh.(unitArea).(populationArea){recInd}.fullInterpCoherence.phase(iUnit,:))
            %end
            %plot(spikingSpikingCoh.(unitArea).(populationArea){recInd}.fullInterpCoherence.frequency(iUnit,:), ...
            %  datamean(spikingSpikingCoh.(unitArea).(populationArea){recInd}.fullInterpCoherence.phase, ...
            %  'circularNP'), 'LineWidth',2)
            %hold off
            %ax = gca;
            %ax.XScale = 'log';
            %ylim([-pi pi]);
            %fTitle = strrep([unitArea ' wrt ' populationArea ': rec' num2str(recInd)], '_', '-');
            %title(fTitle)
            %try; close(fH); catch; end; %#ok<NOSEMI>
          end
        end
      end
    end
  end

  % Save data analysis results
  if strcmpi(population, 'all')
    if excludeMovement
      infraslowAnalyses.spikingSpikingCoh_noMovement.(unitArea) = spikingSpikingCoh.(unitArea);
    else
      infraslowAnalyses.spikingSpikingCoh.(unitArea) = spikingSpikingCoh.(unitArea);
    end
  elseif strcmpi(population, 'positive')
    if excludeMovement
      infraslowAnalyses.spikingSpikingCohPositive_noMovement.(unitArea) = spikingSpikingCoh.(unitArea);
    else
      infraslowAnalyses.spikingSpikingCohPositive.(unitArea) = spikingSpikingCoh.(unitArea);
    end
  elseif strcmpi(population, 'negative')
    if excludeMovement
      infraslowAnalyses.spikingSpikingCohNegative_noMovement.(unitArea) = spikingSpikingCoh.(unitArea);
    else
      infraslowAnalyses.spikingSpikingCohNegative.(unitArea) = spikingSpikingCoh.(unitArea);
    end
  end
  save(analysisResultsFile, 'infraslowAnalyses', '-v7.3');
end

% Switch on warnings
if parallelCores > 1
  parfevalOnAll(@warning,0,'on','all');
end
warning('on', 'all');