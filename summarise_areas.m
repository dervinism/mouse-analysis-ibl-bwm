% Use this script to summarise brain areas and produce a brain area table.
% The script can only be executed on Matlab version R2021a or higher due to
% the use of functions with keyword-value argument pairs.

% Load parameters
params

% Load preprocessed data
preprocessedDataFile = fullfile(processedDataFolder, 'bwmPreprocessedData2.mat');
load(preprocessedDataFile);

% Extract area data
nExp = numel(infraslowData.experimentData);
brainAreas = [];
unitCounts = [];
unitInds = {};
goodUnitInds = {};
areaCount = 0;
for iExp = 1:nExp
  expData = infraslowData.experimentData{iExp};
  if isfield(expData, 'unitBrainAreas')
    expBrainAreas = unique(expData.unitBrainAreas);
    nAreas = numel(expBrainAreas);
    expUnitCounts = zeros(nAreas,2);
    for iArea = 1:nAreas
      areaCount = areaCount + 1;
      unitAreaAllInds = ismember(expData.unitBrainAreas, expBrainAreas{iArea});
      unitAreaGoodInds = unitAreaAllInds & expData.goodUnits;
      expUnitCounts(iArea,:) = [sum(unitAreaAllInds) sum(unitAreaGoodInds)];
      unitAreaAllInds = find(unitAreaAllInds);
      unitInds{areaCount} = [repmat(iExp, numel(unitAreaAllInds), 1) unitAreaAllInds]; %#ok<*SAGROW>
      unitAreaGoodInds = find(unitAreaGoodInds);
      goodUnitInds{areaCount} = [repmat(iExp, numel(unitAreaGoodInds), 1) unitAreaGoodInds];
    end
    brainAreas = [brainAreas; expBrainAreas]; %#ok<*AGROW>
    unitCounts = [unitCounts; expUnitCounts];
  end
end

% Sum units
uniqueBrainAreas = unique(brainAreas);
nAreas = numel(uniqueBrainAreas);
unitCountSum = zeros(nAreas,2);
groupedUnitInds = cell(nAreas,1);
groupedGoodUnitInds = cell(nAreas,1);
for iArea = 1:nAreas
  areaInds = ismember(brainAreas, uniqueBrainAreas{iArea});
  unitCountSum(iArea,:) = sum(unitCounts(areaInds,:));
  groupedUnitInds{iArea} = concatenateCells(unitInds(areaInds));
  groupedGoodUnitInds{iArea} = concatenateCells(goodUnitInds(areaInds));
end

% Load full area names
% Refer to: https://connectivity.brain-map.org/3d-viewer?v=1
fullAreaList = getAllenStructureList;
fullAreaNames = cell(nAreas,1);
for iArea = 1:nAreas
  try
    fullAreaNames{iArea} = fullAreaList.name{ismember(fullAreaList.acronym, ...
      uniqueBrainAreas{iArea})};
  catch
    try
      fullAreaNames{iArea} = fullAreaList.name{ismember(fullAreaList.acronym, ...
        strrep(uniqueBrainAreas{iArea}, ' ', ', '))};
    catch
      if strcmpi(uniqueBrainAreas{iArea}, 'void')
        fullAreaNames{iArea} = uniqueBrainAreas{iArea};
      else
        error('Unknown brain area acronym');
      end
    end
  end
end

% Area table
areaTable = table(areaLabels(:,1), uniqueBrainAreas, fullAreaNames, areaLabels(:,4), unitCountSum(:,1), unitCountSum(:,2), ...
  'VariableNames', {'Brain_area_group', 'Brain_area_acronym', 'Brain_area_name', 'Brain_area_type', 'Total_unit_count', 'Good_unit_count'});

% Save area table and unit indices
analysisResultsFile = fullfile(processedDataFolder, 'bwmAnalysisResults.mat');
if exist(analysisResultsFile, 'file')
  load(analysisResultsFile);
end
infraslowAnalyses.areaSummaries.areaTable = areaTable;
infraslowAnalyses.areaSummaries.groupedUnitInds = groupedUnitInds;
infraslowAnalyses.areaSummaries.groupedGoodUnitInds = groupedGoodUnitInds;
save(analysisResultsFile, 'infraslowAnalyses', '-v7.3');