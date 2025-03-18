% Script for producing phase profile figures with pupil area as the
% reference. The script can only be executed on Matlab version R2021a or
% higher due to the use of functions with keyword-value argument pairs.

% Load parameters
params
population = 'negative'; %'all', 'positive', 'negative'
significantOnly = false;
%areas = {'VPL_VPM_LG_PO_LP', 'SSp', 'RSP', 'CA_DG'};
%areas = {'sensory_Th', 'sensory_nCx', 'paCx', 'CA_DG'};
%areas = {'sensory_Th', 'motor_Th', 'association_Th'};
%areas = {'sensory_nCx', 'motor_nCx', 'association_nCx'};
%areas = {'Th', 'nCx', 'RSPd', 'CA_DG'};
areas = {'Th', 'Cx', 'CA_DG'};
%areas = {'sensory_Th', 'association_Th', 'sensory_nCx', 'association_nCx', 'CA_DG'};

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
      if significantOnly
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohSignificant); %#ok<*UNRCH>
        compData = [compData, infraslowAnalyses.spikingPupilCohSignificant.(structFieldNames{areaInds(iArea)})];
      else
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCoh);
        compData = [compData, infraslowAnalyses.spikingPupilCoh.(structFieldNames{areaInds(iArea)})];
      end
    elseif strcmpi(population, 'positive')
      if significantOnly
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohPositiveSignificant);
        compData = [compData, infraslowAnalyses.spikingPupilCohPositiveSignificant.(structFieldNames{areaInds(iArea)})];
      else
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohPositive);
        compData = [compData, infraslowAnalyses.spikingPupilCohPositive.(structFieldNames{areaInds(iArea)})];
      end
    elseif strcmpi(population, 'negative')
      if significantOnly
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohNegativeSignificant);
        compData = [compData, infraslowAnalyses.spikingPupilCohNegativeSignificant.(structFieldNames{areaInds(iArea)})];
      else
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohNegative);
        compData = [compData, infraslowAnalyses.spikingPupilCohNegative.(structFieldNames{areaInds(iArea)})];
      end
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
    p = zeros(nComp, 1);
  elseif iComp == 2
    colour = 'r';
  elseif iComp == 3
    colour = 'b';
  elseif iComp == 4
    colour = 'c';
  elseif iComp == 5
    colour = 'm';
  end
  ciplot(phaseProfile+phaseProfileCI(1,:), ...
    phaseProfile+phaseProfileCI(2,:), freq, colour, 0.1);
  ax = gca;
  ax.XScale = 'log';
  if iComp == 1
    hold on
  end
  p(iComp) = semilogx(freq, phaseProfile, colour, 'LineWidth',1.5);
  if iComp == nComp
    hold off
    legend(p, areas);
  end
end
ylim([-pi pi]);