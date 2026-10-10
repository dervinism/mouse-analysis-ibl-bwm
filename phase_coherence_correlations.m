% Script for producing phase-coherence correlations. The script can only be
% executed on Matlab version R2021a or higher due to the use of functions
% with keyword-value argument pairs.

% Load parameters
params
FOI = FOI(12:end);
population = 'all'; %'all', 'positive', 'negative'
significantOnly = false;
%areas = {'DG'};
areas = areasOI;
freqOI = 0.3;
saveUnitTable = true;
drawFigs = false;

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

% Go through area comparisons and produce figures
nComps = numel(areas);
nFreqs = numel(freqOI);

% Loop over areas
clear unitTable
for iComp = 1:nComps

  % Get the data
  areaInds = find(getAreaInds(strrep(areas{iComp}, '_', '-'), infraslowAnalyses.areaSummaries.areaTable));
  compData = {};
  for iArea = 1:numel(areaInds)
    if strcmpi(population, 'all')
      if significantOnly
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohSignificant); %#ok<*UNRCH>
        ind = min([numel(structFieldNames) areaInds(iArea)]);
        compData = [compData, infraslowAnalyses.spikingPupilCohSignificant.(structFieldNames{ind})];
      else
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCoh);
        ind = min([numel(structFieldNames) areaInds(iArea)]);
        compData = [compData, infraslowAnalyses.spikingPupilCoh.(structFieldNames{ind})];
      end
    elseif strcmpi(population, 'positive')
      if significantOnly
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohPositiveSignificant);
        ind = min([numel(structFieldNames) areaInds(iArea)]);
        compData = [compData, infraslowAnalyses.spikingPupilCohPositiveSignificant.(structFieldNames{ind})];
      else
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohPositive);
        ind = min([numel(structFieldNames) areaInds(iArea)]);
        compData = [compData, infraslowAnalyses.spikingPupilCohPositive.(structFieldNames{ind})];
      end
    elseif strcmpi(population, 'negative')
      if significantOnly
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohNegativeSignificant);
        ind = min([numel(structFieldNames) areaInds(iArea)]);
        compData = [compData, infraslowAnalyses.spikingPupilCohNegativeSignificant.(structFieldNames{ind})];
      else
        structFieldNames = fieldnames(infraslowAnalyses.spikingPupilCohNegative);
        ind = min([numel(structFieldNames) areaInds(iArea)]);
        compData = [compData, infraslowAnalyses.spikingPupilCohNegative.(structFieldNames{ind})];
      end
    end
  end

  % Combine recording phases and coherences
  nRecs = numel(compData);
  phase = [];
  coherence = [];
  firingRate = [];
  rSpearman = [];
  pvalSpearman = [];
  for iRec = 1:nRecs
    if isfield(compData{iRec}, 'fullInterpCoherence')
      recData = compData{iRec}.fullInterpCoherence;
      %freqInds = ismember(recData.frequency(1,:), freqOI);
      freqInds = ismember(recData.frequency(1,:), FOI);
      phase = [phase; recData.phase(:,freqInds)]; %#ok<*AGROW>
      coherence = [coherence; recData.rateAdjustedCoherence(:,freqInds)];
      if saveUnitTable
        firingRate = [firingRate; compData{iRec}.firingRates];
        rSpearman = [rSpearman; compData{iRec}.rSpearman];
        pvalSpearman = [pvalSpearman; compData{iRec}.pvalSpearman];
      end
    end
  end
  freq = recData.frequency(1, freqInds);

  % Draw the figure
  if drawFigs
    fH = figure;
    plot(phase, coherence, '.', 'MarkerSize',5)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(phase-2*pi, coherence, '.', 'MarkerSize',5)
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(phase+2*pi, coherence, '.', 'MarkerSize',5)
    hold off
    xlim([-2*pi 2*pi]);
    xlabel('Phase (rad)');
    ylabel('Coherence');
    legend('0.01 Hz', '0.02 Hz', '0.03 Hz', '0.05 Hz', '0.07 Hz', '0.10 Hz', ...
      '0.20 Hz', '0.30 Hz', '0.50 Hz', '0.70 Hz', '1.00 Hz');
    if strcmpi(population, 'all')
      %tH = title(['Phase_v_coherence__' areas{iComp} '_' num2str(freqOI) 'Hz__all_units']);
      tT = ['Phase_v_coherence__' areas{iComp} '_wideband__all_units'];
      tH = title(tT, 'Interpreter','none');
    elseif strcmpi(population, 'positive')
      tT = ['Phase_v_coherence__' areas{iComp} '_' num2str(freqOI) 'Hz__positive_units'];
      tH = title(tT, 'Interpreter','none');
    elseif strcmpi(population, 'negative')
      tT = ['Phase_v_coherence__' areas{iComp} '_' num2str(freqOI) 'Hz__negative_units'];
      tH = title(tT, 'Interpreter','none');
    end

    % Save the figure
    savefig(fH, [tT '.fig'], 'compact');
    close(fH);
  end

  % Create the unit table
  if saveUnitTable
    unitAreas = repmat(areas(iComp), numel(rSpearman), 1);
    compUnitTable = table(unitAreas, firingRate, rSpearman, pvalSpearman, ...
      'VariableNames', {'Area', 'Rate', 'rSpearman', 'pvalSpearman'});
    for iF = 1:numel(FOI)
      label = ['Coh_' strrep(num2str(round(FOI(end-iF+1),4)), '.', 'p')];
      compUnitTable.(label) = coherence(:,iF);
    end
    for iF = 1:numel(FOI)
      label = ['Phi_' strrep(num2str(round(FOI(end-iF+1),4)), '.', 'p')];
      compUnitTable.(label) = phase(:,iF);
    end
    if ~exist('unitTable','var')
      unitTable = compUnitTable;
    else
      unitTable = [unitTable; compUnitTable]; %#ok<*AGROW>
    end
  end
end

% Save the unit table
if saveUnitTable
  save('phase_firingRate_table_IBL.mat', 'unitTable', '-v7.3');
end