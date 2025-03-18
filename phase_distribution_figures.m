% Script for producing figures showing distributions of preferred unit
% firing phase wrt the pupila area size. The script can only be executed on
% Matlab version R2021a or higher due to the use of functions with
% keyword-value argument pairs.

% Load parameters
params
FOI = fliplr(FOI);
selectAreaHistos = false;
allAreaHistos = false;
oldAreaHistos = true;

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
  [binCounts, binLocs, totalCounts, significantFractions, ...
    phaseMeans, phaseSDs] = phaseHistrogramWrap(areaOIPhase)

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
  [binCounts, binLocs, totalCounts, significantFractions, ...
    phaseMeans, phaseSDs] = phaseHistrogramWrap(areaPhase)

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

if oldAreaHistos
  %[binCounts, binLocs, totalCounts, significantFractions, ...
  %  phaseMeans, phaseSDs] = phaseHistrogramWrap(areaOIPhase);

  options = struct();
  options.mainFolder = figFolder;
  options.histosSubfolder = 'phaseHistos';
  options.mapsSubfolder = 'phaseMaps';
  options.figSize = 15;
  options.figTitle = 'PUPIL';
  options.freqLim = [10e-3 2];
  phaseShift = pi/2;
  options.phaseLimHisto = phaseLim + phaseShift;
  options.phaseLimMap = phaseLim - phaseShift;
  options.xLabelHist = '# units';
  options.limExpansion = pi/2;
  options.mask = {[options.freqLim(1) options.phaseLimMap(end)+options.limExpansion; 2 options.phaseLimMap(end)+options.limExpansion;...
    options.freqLim(1) pi/2; 2 options.phaseLimMap(end)+options.limExpansion];...
    [0.05 options.phaseLimMap(1)-options.limExpansion; options.freqLim(end) -pi/2;...
    0.05 options.phaseLimMap(1)-options.limExpansion; options.freqLim(end) options.phaseLimMap(1)-options.limExpansion]};
  options.iAreasOI = 1:numel(areasOI);
  [phaseHistos, distributionStats] = phaseHistosPlotMaster([false false true], ...
    areasOI, {'awake'}, FOI, {areaOIPhase}, edges+phaseShift, options);
end



%% Local functions
function [binCounts, binLocs, totalCounts, significantFractions, ...
  phaseMeans, phaseSDs] = phaseHistrogramWrap(areaPhases) %#ok<*DEFNU>
% [binCounts, binLocs, totalCounts, significantFractions, ...
%   phaseMeans, phaseSDs] = phaseHistrogramWrap(areaPhases)
%
% Function calculates phase histogram profiles. It's a helper function
% within phase_distribution_figures.m script.

% Bin phase values to histograms for areas of interest
nAreas = numel(areaPhases);
binCounts = cell(size(areaPhases)); %#ok<*UNRCH>
binLocs = cell(size(areaPhases));
totalCounts = cell(size(areaPhases));
phaseMeans = cell(size(areaPhases));
phaseSDs = cell(size(areaPhases));
significantFractions = cell(size(areaPhases));
for iArea = 1:nAreas
  if ~isempty(areaPhases{iArea})
    [binCounts{iArea}, binLocs{iArea}, totalCounts{iArea}, ...
      significantFractions{iArea}, phaseMeans{iArea}, phaseSDs{iArea}] = ...
      phaseHistrogram(areaPhases{iArea}, centre=0, nBins=10);
  end
end
end