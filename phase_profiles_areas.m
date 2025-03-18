% Script for producing phase profile figures for cross-area comparisons.
% The script can only be executed on Matlab version R2021a or higher due to
% the use of functions with keyword-value argument pairs.

% Load parameters
params
population = 'positive'; %'all', 'positive', 'negative'
%areas = {'VPL_VPM_LG_PO_LP', 'SSp', 'RSP', 'CA_DG'};
areas = {'VPL_VPM_LG_PO_LP', 'SSp', 'CA_DG'};
%areas = {'Th', 'nCx', 'RSP', 'CA_DG'};
direction = 'forward'; %'forward', 'backward'

% Load data analysis results
if ~exist('infraslowAnalyses', 'var')
  analysisResultsFile = fullfile(processedDataFolder, 'bwmAnalysisResults.mat');
  load(analysisResultsFile); 
end

% Go through area comparisons and produce figures
nComp = numel(areas)-1;

% Confidence intervals
for iComp = 1:nComp

  % Get the data
  if strcmpi(direction, 'forward')
    areaSignal = areas{1};
    areaReference = areas{iComp+1};
  elseif strcmpi(direction, 'backward')
    areaSignal = areas{iComp+1};
    areaReference = areas{1};
  end
  if strcmpi(population, 'all')
    compData = infraslowAnalyses.spikingSpikingCoh.(areaSignal).(areaReference);
  elseif strcmpi(population, 'positive')
    compData = infraslowAnalyses.spikingSpikingCohPositive.(areaSignal).(areaReference);
  elseif strcmpi(population, 'negative')
    compData = infraslowAnalyses.spikingSpikingCohNegative.(areaSignal).(areaReference);
  end

  % Combine recording phases
  nRecs = numel(compData);
  phase = [];
  for iRec = 1:nRecs
    recData = compData{iRec}.fullInterpCoherence;
    freqInds = ismember(recData.frequency(1,:), FOI);
    phase = [phase; recData.phase(:,freqInds)]; %#ok<*AGROW>
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