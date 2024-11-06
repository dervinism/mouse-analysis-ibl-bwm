% Parameters
samplingRate = 60;
binSpikes = true;

% I/O
preprocessedDataFile = 'C:\Users\44079\Work\Leicester\infraslow-dynamics\04_data_analysis\005_ibl_bwm\bwmPreprocessedData.mat';
outputDataFile = 'C:\Users\44079\Work\Leicester\infraslow-dynamics\04_data_analysis\005_ibl_bwm\bwmPreprocessedData2.mat';

% Load preprocessed data
load(preprocessedDataFile);

% Convert preprocessed data into a more structured format
nExperiments = numel(dataStruct.experiment_data);
infraslowData.experimentData = cell(nExperiments,1);
infraslowData.probeTable = dataStruct.probeTable;
infraslowData.unitTable = dataStruct.unitTable;
for iExp = 1:nExperiments

  % Convert electrophysiology data
  probe0Data = extractEphysData(dataStruct.experiment_data{iExp}.probe0, ...
    dataStruct.experiment_data{iExp}.spontaneous_activity_times, samplingRate, binSpikes);
  probe1Data = extractEphysData(dataStruct.experiment_data{iExp}.probe0, ...
    dataStruct.experiment_data{iExp}.spontaneous_activity_times, samplingRate, binSpikes);
  if ~isempty(probe0Data) && isempty(probe1Data)
    expData = probe0Data;
  elseif isempty(probe0Data) && ~isempty(probe1Data)
    expData = probe1Data;
  elseif ~isempty(probe0Data) && ~isempty(probe1Data)
    expData = probe0Data;
    expData.units = [expData.units; probe1Data.units];
    expData.nUnits = expData.nUnits + probe1Data.nUnits;
    if binSpikes
      expData.spikeCounts = sparse([expData.spikeCounts; probe1Data.spikeCounts]);
    else
      expData.spikeTimes = [expData.spikeTimes; probe1Data.spikeTimes]; %#ok<*UNRCH>
    end
    expData.unitChannels = [expData.unitChannels; probe1Data.unitChannels];
    expData.unitWaveforms = [expData.unitWaveforms; probe1Data.unitWaveforms];
    expData.unitUIDs = [expData.unitUIDs; probe1Data.unitUIDs];
    expData.horzCoords = [expData.horzCoords; probe1Data.horzCoords];
    expData.vertCoords = [expData.vertCoords; probe1Data.vertCoords];
    expData.unitBrainLocs = [expData.unitBrainLocs; probe1Data.unitBrainLocs];
    expData.unitChannelLabels = [expData.unitChannelLabels; probe1Data.unitChannelLabels];
    expData.unitBrainAreas = [expData.unitBrainAreas; probe1Data.unitBrainAreas];
    expData.probeIDs = [expData.probeIDs; probe1Data.probeIDs];
  else
    expData = struct();
  end
  if ~isempty(expData) && ~isempty(fieldnames(expData))
    expData.goodUnits = ismember(expData.unitUIDs, dataStruct.unitTable.uuids);
    expData.samplingRate = samplingRate;
  end

  % Convert video data
  if ~isempty(expData) && ~isempty(fieldnames(expData))
    leftCameraData = extractVideoData(dataStruct.experiment_data{iExp}.left_camera, ...
      expData.spikeTimeBins + dataStruct.experiment_data{1, 2}.spontaneous_activity_times(1));
    rightCameraData = extractVideoData(dataStruct.experiment_data{iExp}.right_camera, ...
      expData.spikeTimeBins + dataStruct.experiment_data{1, 2}.spontaneous_activity_times(1));
    if ~isempty(leftCameraData)
      expData.leftPupilAreaSize = leftCameraData.pupilAreaSize;
    else
      expData.leftPupilAreaSize = [];
    end
    if ~isempty(rightCameraData)
      expData.rightPupilAreaSize = rightCameraData.pupilAreaSize;
    else
      expData.rightPupilAreaSize = [];
    end
  end

  % Attach the remaining data
  expData.recordingTimeRange = ...
    dataStruct.experiment_data{iExp}.spontaneous_activity_times;
  expData.experimentID = dataStruct.experiment_data{iExp}.eid;
  if numel(dataStruct.experiment_data{iExp}.subject_id) > 1
    expData.subjectID = dataStruct.experiment_data{iExp}.subject_id(2);
  else
    expData.subjectID = [];
  end

  infraslowData.experimentData{iExp} = expData;
end

% Save the converted data
save(outputDataFile, 'infraslowData', '-v7.3');



%% Local functions
function unitEphysData = extractEphysData(probeEphysData, timeRange, samplingRate, binSpikes)
% unitEphysData = extractEphysData(probeEphysData, timeRange, samplingRate, binSpikes)
%
% A local helper function to convert_to_infraslow_data_format script.

if ~isempty(probeEphysData) && ~isempty(fieldnames(probeEphysData))

  % Timing
  timeStep = 1/samplingRate;
  startTime = round((timeRange(1) + timeStep/2)/timeStep)*timeStep;
  endTime = ceil(timeRange(2)/timeStep)*timeStep;

  % Data extraction
  unitEphysData = struct();
  unitEphysData.units = unique(probeEphysData.spike_clusters)';
  unitEphysData.nUnits = numel(unitEphysData.units);
  unitEphysData.spikeTimeBins = (startTime:timeStep:endTime) - startTime;
  if binSpikes
    unitEphysData.spikeCounts = zeros(unitEphysData.nUnits, numel(unitEphysData.spikeTimeBins));
  else
    unitEphysData.spikeTimes = cell(unitEphysData.nUnits,1);
  end
  for iUnit = 1:unitEphysData.nUnits
    spikeInds = probeEphysData.spike_clusters == unitEphysData.units(iUnit);
    if binSpikes
      unitSpikeCounts = resampleSpikes(probeEphysData.spike_times(spikeInds), ...
        startTime=startTime, stepsize=timeStep);
      unitEphysData.spikeCounts(iUnit,1:numel(unitSpikeCounts)) = unitSpikeCounts;
    else
      unitEphysData.spikeTimes{iUnit} = probeEphysData.spike_times(spikeInds) ...
        - startTime - timeStep/2;
    end
  end
  originalUnitIDs = 0:numel(probeEphysData.cluster_uuids)-1;
  remainingUnitInds = ismember(originalUnitIDs, unitEphysData.units);
  unitEphysData.unitChannels = probeEphysData.cluster_channels(remainingUnitInds)' + 1;
  unitEphysData.unitWaveforms = probeEphysData.cluster_waveforms(remainingUnitInds,:,:);
  unitEphysData.unitUIDs = probeEphysData.cluster_uuids(remainingUnitInds);
  unitEphysData.horzCoords = ...
    probeEphysData.channels_localCoordinates(unitEphysData.unitChannels,1);
  unitEphysData.vertCoords = ...
    probeEphysData.channels_localCoordinates(unitEphysData.unitChannels,2);
  unitEphysData.unitBrainLocs = ...
    probeEphysData.channels_brainLocationIds_ccf_2017(unitEphysData.unitChannels)';
  unitEphysData.unitChannelLabels = ...
    probeEphysData.channels_labels(unitEphysData.unitChannels)';
  unitEphysData.unitBrainAreas = ...
    probeEphysData.channels_brain_areas(unitEphysData.unitChannels)';
  unitEphysData.probeIDs = repmat(probeEphysData.pid(2), unitEphysData.nUnits, 1);
else
  unitEphysData = [];
end
end


function convVideoData = extractVideoData(videoData, ephysTimes)
% convVideoData = extractVideoData(videoData, ephysTimes)
%
% A local helper function to convert_to_infraslow_data_format script.

if ~isempty(videoData) && ~isempty(fieldnames(videoData)) && ...
    ~isempty(videoData.pupilDiameter_smooth)

  % Data extraction
  convVideoData = struct();
  convVideoData.pupilAreaSize = interp1(videoData.times, ...
    double(videoData.pupilDiameter_smooth(:,2)), ephysTimes);
  convVideoData.pupilAreaSize = (convVideoData.pupilAreaSize./2).^2;
else
  convVideoData = [];
end
end