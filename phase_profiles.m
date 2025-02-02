% Script for producing phase profile figures. The script can only be
% executed on Matlab version R2021a or higher due to the use of functions
% with keyword-value argument pairs.

% Load parameters
params
population = 'all'; %'all', 'positive', 'negative'
areas = {'VPL_VPM_LG_PO_LP', 'SSp', 'RSP', 'CA_DG'};
%areas = {'DG', 'SSp', 'RSP', 'VPL_VPM_LG_PO_LP'};

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
  if strcmpi(population, 'all')
    compData = infraslowAnalyses.spikingSpikingCoh.(areas{iComp+1}).(areas{1});
  elseif strcmpi(population, 'positive')
    compData = infraslowAnalyses.spikingSpikingCoh.(areas{iComp+1}).(areas{1});
  elseif strcmpi(population, 'negative')
    compData = infraslowAnalyses.spikingSpikingCoh.(areas{iComp+1}).(areas{1});
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
  disp(phaseProfile);
end
ylim([-pi pi]);