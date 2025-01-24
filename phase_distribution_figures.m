% Script for producing figures showing distributions of preferred unit
% firing phase wrt the pupila area size. The script can only be executed on
% Matlab version R2021a or higher due to the use of functions with
% keyword-value argument pairs.

% Load parameters
params
FOI = fliplr(FOI);
selectAreaHistos = true;
allAreaHistos = false;

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
for iArea = 1:size(areaLabels,1)
  areaName = strrep(areaLabels{iArea,2}, ' ', '_');
  areaName = strrep(areaName, '-', '_');
  areaPhase{iArea} = [];
  if ~strcmpi(areaName, 'root') && ~strcmpi(areaName, 'void')
    nRecs = numel(infraslowAnalyses.spikingPupilCoh.(areaName));
    for iRec = 1:nRecs
      if isfield(infraslowAnalyses.spikingPupilCoh.(areaName){iRec}, 'fullInterpCoherence')
        phase = infraslowAnalyses.spikingPupilCoh.(areaName){iRec}.fullInterpCoherence.phase;
        areaFrequencies = infraslowAnalyses.spikingPupilCoh.(areaName){iRec}.fullInterpCoherence.frequency(1,:);
        fInds = ismember(areaFrequencies, FOI);

        % Agregate data
        areaPhase{iArea} = [areaPhase{iArea}; phase(:,fInds)];
      end
    end
  end
end

% Agregate data for areas of interest
% Select areas of interest
nAreas = numel(areasOI);
areaOIPhase = cell(nAreas,1);
for iArea = 1:nAreas
  areaInds = getAreaInds(areasOI{iArea}, infraslowAnalyses.areaSummaries.areaTable);
  areaOIPhase{iArea} = concatenateCells(areaPhase(areaInds));
  %disp([areasOI{iArea} ': ' num2str(sum(areaInds))])
end

if selectAreaHistos
  % Bin phase values to histograms for areas of interest
  binCounts = cell(size(areaOIPhase)); %#ok<*UNRCH>
  binLocs = cell(size(areaOIPhase));
  totalCounts = cell(size(areaOIPhase));
  phaseMeans = cell(size(areaOIPhase));
  phaseSDs = cell(size(areaOIPhase));
  significantFractions = cell(size(areaOIPhase));
  for iArea = 1:nAreas
    if ~isempty(areaOIPhase{iArea})
      [binCounts{iArea}, binLocs{iArea}, totalCounts{iArea}, ...
        significantFractions{iArea}, phaseMeans{iArea}, phaseSDs{iArea}] = ...
        phaseHistrogram(areaOIPhase{iArea}, centre=0, nBins=10);
    end
  end

  % Plot phase histograms for areas of interest
  if ~exist(figFolder, 'dir')
    mkdir(figFolder);
  end
  for iArea = 1:nAreas
    for iFreq = 1:numel(FOI)
      if ~isempty(binCounts{iArea})
        nValues = round(totalCounts{iArea}(iFreq)/significantFractions{iArea}(iFreq));
        if isnan(nValues)
          nValues = 0;
        end
        text2display = {[num2str(totalCounts{iArea}(iFreq)) '/' num2str(nValues)]};
        figTitle = ['Phase histogram: ' areasOI{iArea} ' ' num2str(FOI(iFreq)) 'Hz'];
        fH = phaseHistogramPlot(binCounts{iArea}(:,iFreq), ...
          dataMean=phaseMeans{iArea}(iFreq), figText=text2display, ...
          figTitle=figTitle, figPath=figFolder);
        close(fH);
      end
    end
  end
end

if allAreaHistos
  % Bin phase values to histograms for all areas
  nAreas = size(areaLabels,1);
  binCounts = cell(size(areaPhase));
  binLocs = cell(size(areaPhase));
  totalCounts = cell(size(areaPhase));
  phaseMeans = cell(size(areaPhase));
  phaseSDs = cell(size(areaPhase));
  significantFractions = cell(size(areaPhase));
  for iArea = 1:nAreas
    if ~isempty(areaPhase{iArea})
      [binCounts{iArea}, binLocs{iArea}, totalCounts{iArea}, ...
        significantFractions{iArea}, phaseMeans{iArea}, phaseSDs{iArea}] = ...
        phaseHistrogram(areaPhase{iArea}, centre=0, nBins=10);
    end
  end

  % Plot phase histograms for all areas
  if ~exist(figFolder, 'dir')
    mkdir(figFolder);
  end
  for iArea = 1:nAreas
    for iFreq = 1:numel(FOI)
      if ~isempty(binCounts{iArea})
        nValues = round(totalCounts{iArea}(iFreq)/significantFractions{iArea}(iFreq));
        if isnan(nValues)
          nValues = 0;
        end
        text2display = {[num2str(totalCounts{iArea}(iFreq)) '/' num2str(nValues)]};
        figTitle = ['Phase histogram: ' areaLabels{iArea,2} ' ' num2str(FOI(iFreq)) 'Hz'];
        fH = phaseHistogramPlot(binCounts{iArea}(:,iFreq), ...
          dataMean=phaseMeans{iArea}(iFreq), figText=text2display, ...
          figTitle=figTitle, figPath=figFolder);
        close(fH);
      end
    end
  end
end