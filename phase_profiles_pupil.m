% Script for producing phase profile figures with pupil area as the
% reference. The script can only be executed on Matlab version R2021a or
% higher due to the use of functions with keyword-value argument pairs.

% Load parameters
params
population = 'all'; %'all', 'positive', 'negative'
areas = {'VPL_VPM_LG_PO_LP', 'SSp', 'RSP', 'CA_DG'};

% Load data analysis results
if ~exist('infraslowAnalyses', 'var')
  analysisResultsFile = fullfile(processedDataFolder, 'bwmAnalysisResults.mat');
  load(analysisResultsFile); 
end

% Go through area comparisons and produce figures
nComp = numel(areas);

% Confidence intervals
for iComp = 1:nComp

  % Get the data
  areaInds = find(getAreaInds(strrep(areas{iComp}, '_', '-'), infraslowAnalyses.areaSummaries.areaTable));
  compData = {};
  for iArea = 1:numel(areaInds)
    if strcmpi(population, 'all')
      structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCoh);
      compData = [compData, infraslowAnalyses.spikingPupilCoh.(structFieldNames{areaInds(iArea)})];
    elseif strcmpi(population, 'positive')
      structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohPositive);
      compData = [compData, infraslowAnalyses.spikingPupilCohPositive.(structFieldNames{areaInds(iArea)})];
    elseif strcmpi(population, 'negative')
      structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohNegative);
      compData = [compData, infraslowAnalyses.spikingPupilCohNegative.(structFieldNames{areaInds(iArea)})];
    end
  end

  % Combine recording phases
  nRecs = numel(compData);
  phase = [];
  for iRec = 1:nRecs
    if isfield(compData{iRec}, 'fullInterpCoherence')
      recData = compData{iRec}.fullInterpCoherence;
      freqInds = ismember(recData.frequency(1,:), FOI);
      phase = [phase; recData.phase(:,freqInds)]; %#ok<*AGROW>
    end
  end
  [phaseProfile, phaseProfileCI] = datamean(phase);
  freq = recData.frequency(1, freqInds);

  % Draw CIs
  if iComp == 1
    fH = figure;
    colour = 'g';
  elseif iComp == 2
    colour = 'r';
  elseif iComp == 3
    colour = 'b';
  elseif iComp == 4
    colour = 'c';
  end
  ciplot(phaseProfile+phaseProfileCI(1,:), ...
    phaseProfile+phaseProfileCI(2,:), freq, colour, 0.1);
  ax = gca;
  ax.XScale = 'log';
  if iComp == 1
    hold on
  end
  semilogx(freq, phaseProfile, colour, 'LineWidth',1.5);
  if iComp == nComp
    hold off
  end
  %disp(phaseProfile);
end
ylim([-pi pi]);