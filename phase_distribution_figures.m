% Script for producing figures showing distributions of preferred unit
% firing phase wrt the pupila area size. The script can only be executed on
% Matlab version R2021a or higher due to the use of functions with
% keyword-value argument pairs.

% Load parameters
params

% Load data analysis results
analysisResultsFile = fullfile(processedDataFolder, 'bwmAnalysisResults.mat');
if ~exist('infraslowAnalyses', 'var') && exist(analysisResultsFile, 'file')
  load(analysisResultsFile);
elseif ~exist('infraslowAnalyses', 'var')
  error('Missing infraslow analyses file.')
end

% Extract data
nAreas = size(areaLabels,1);
areaPhase = cell(nAreas,1);
areaPhaseOI = cell(size(areasOI));
for iArea = 1:numel(areaLabels)
  areaName = areaLabels{iArea,2};
  areaPhase{iArea} = [];
  nRecs = numel(infraslowAnalyses.spikingPupilCoh.(areaName));
  for iRec = 1:nRecs
    phase = infraslowAnalyses.spikingPupilCoh.(areaName){iRec}.fullInterpCoherence.phase;
    areaFrequencies = infraslowAnalyses.spikingPupilCoh.(areaName){iRec}.fullInterpCoherence.frequency(1,:);
    fInds = ismember(areaFrequencies, FOI);

    % Agregate data
    areaPhase{iArea} = [areaPhase{iArea}; phase(:,fInds)];
  end
end

% Agregate data for areas of interest
% Select areas of interest
nAreas = numel(areasOI);
areaOIPhase = cell(nAreas,1);
for iArea = 1:nAreas
  areaInds = getAreaInds(areasOI{iArea}, infraslowAnalyses.areaSummaries.areaTable);
  areaOIPhase{iArea} = concatenateCells(areaPhase(areaInds));
end

% Bin phase values to histograms
binCounts = cell(size(areaPhase));
binLocs = cell(size(areaPhase));
totalCounts = cell(size(areaPhase));
phaseMeans = cell(size(areaPhase));
phaseSDs = cell(size(areaPhase));
significantFractions = cell(size(areaPhase));
for iArea = 1:numel(areaPhase)
  [binCounts{iArea}, binLocs{iArea}, totalCounts{iArea}, ...
    significantFractions{iArea}, phaseMeans{iArea}, phaseSDs{iArea}] = ...
    phaseHistrogram(areaPhase{iArea}, centre=0, nBins=10);
end

% Plot phase histograms
if ~exist(figFolder, 'dir')
  mkdir(figFolder);
end
for iArea = 1:numel(areaPhase)
  for iFreq = 1:numel(FOI)
    nValues = round(totalCounts{iArea}(iFreq)/significantFractions{iArea}(iFreq));
    if isnan(nValues)
      nValues = 0;
    end
    text2display = {[num2str(totalCounts{iArea}(iFreq)) '/' num2str(nValues)]};
    figTitle = ['Phase histogram: ' areaLabels{iArea} ' ' num2str(FOI(iFreq)) 'Hz'];
    fH = phaseHistogramPlot(binCounts{iArea}(:,iFreq), ...
      dataMean=phaseMeans{iArea}(iFreq), figText=text2display, ...
      figTitle=figTitle, figPath=figFolder);
  end
end