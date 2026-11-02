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
%areas = {'Hp', 'DG'};
%areas = {'sensory_Th', 'association_Th', 'sensory_nCx', 'association_nCx', 'CA_DG'};
freqOI = 0.3; %[0.03 0.05 0.07 0.1 0.2 0.3 0.5 0.7 1]; %[0.03 0.3];

drawProfiles = true;
drawSummaries = true;

% Load preprocessed data
if ~exist('infraslowData', 'var')
  preprocessedDataFile = fullfile(processedDataFolder, 'bwmPreprocessedData2.mat');
  load(preprocessedDataFile);
end

% Load data analysis results
if ~exist('infraslowAnalyses', 'var')
  analysisResultsFile = fullfile(processedDataFolder, ['bwmAnalysisResults' resultsTag '.mat']);
  load(analysisResultsFile);
end

if drawProfiles

  % Go through area comparisons and produce figures
  nComps = numel(areas);
  nFreqs = numel(freqOI);

  % Confidence intervals
  for iComp = 1:nComps

    % Get the data
    areaInds = find(getAreaInds(strrep(areas{iComp}, '_', '-'), infraslowAnalyses.areaSummaries.areaTable));
    compData = {};
    for iArea = 1:numel(areaInds)
      if strcmpi(population, 'all')
        if significantOnly
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohSignificant, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        else
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCoh, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        end
      elseif strcmpi(population, 'positive')
        if significantOnly
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohPositiveSignificant, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        else
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohPositive, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        end
      elseif strcmpi(population, 'negative')
        if significantOnly
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohNegativeSignificant, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        else
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohNegative, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
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
    [phaseProfile, phaseProfileCI] = datamean(phase, 'circularNP');
    freq = recData.frequency(1, freqInds);

    % Draw CIs
    if iComp == 1
      fH = figure;
      colour = 'g';
      p = zeros(nComps, 1);
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
    if iComp == nComps
      hold off
      legend(p, areas);
    end
  end
  ylim([-pi pi]);
end


% Summary figures
if drawSummaries

  % Draw ordered phases
  [sortedAreasOI, areaOrderOI] = sort( ...
    infraslowAnalyses.spikingPupilCorr.areasOI.positiveSpearmanFractionsMeans(:,1), 'descend');
  areaOrderOI = areaOrderOI(~isnan(sortedAreasOI));

  % Go through area comparisons and compile them
  nComps = numel(areasOI);
  phaseProfiles = cell(nComps,1);
  phaseProfilesCI = cell(nComps,1);
  freq = cell(nComps,1);
  phaseProfilesOI = cell(nComps,1);
  phaseProfilesOICI = cell(nComps,1);
  phaseFOI = cell(nComps,1);
  for iComp = 1:nComps

    % Get the data
    areaInds = find(getAreaInds(strrep(areasOI{iComp}, '_', '-'), infraslowAnalyses.areaSummaries.areaTable));
    compData = {};
    for iArea = 1:numel(areaInds)
      if strcmpi(population, 'all')
        if significantOnly
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohSignificant, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        else
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCoh, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        end
      elseif strcmpi(population, 'positive')
        if significantOnly
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohPositiveSignificant, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        else
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohPositive, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        end
      elseif strcmpi(population, 'negative')
        if significantOnly
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohNegativeSignificant, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        else
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohNegative, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
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
        phase = [phase; datamean(recData.phase(:,freqInds), 'circularNP')]; %#ok<*SAGROW,*AGROW>
      end
    end
    disp([areasOI{iComp} ' ' num2str(numel(phase))]);
    [phaseProfiles{iComp}, phaseProfilesCI{iComp}] = datamean(phase, 'circularNP');
    freq{iComp} = recData.frequency(1, freqInds);
    freqInds = ismember(freq{iComp}, freqOI);
    phaseProfilesOI{iComp} = phaseProfiles{iComp}(freqInds);
    phaseProfilesOICI{iComp} = phaseProfilesCI{iComp}(:,freqInds);
    for iFreq = 1:nFreqs
      freqInd = ismember(freq{iComp}, freqOI(iFreq));
      phaseFOI{iFreq}{iComp} = phase(:,freqInd);
    end
  end

  % Draw figures
  fontSize = 18;
  for iF = 1:nFreqs

    % Compile phase for areas of interest only
    phase = zeros(1,nComps);
    phaseCI = zeros(2,nComps);
    for iComp = 1:nComps
      phase(iComp) = phaseProfilesOI{iComp}(iF);
      phaseCI(:,iComp) = phaseProfilesOICI{iComp}(:,iF);
    end
    phaseThr = 0;
    phase(phase < phaseThr) = phase(phase < phaseThr) + 2*pi;
    %phaseThr = 4.5;
    %phase(phase > phaseThr) = phase(phase > phaseThr) - 2*pi;

    fH = figure;
    includeAreaInds = [1:14 17:18, 24 29:32 36:47 51];
    dataPointCount = 0;
    for iComp = 1:nComps
      if ismember(areaOrderOI(iComp), includeAreaInds)
        if contains(areasOI(areaOrderOI(iComp)), {'Th', 'VP', 'PO', 'LG', 'LP'}) && ~strcmpi(areasOI(areaOrderOI(iComp)), 'STh')
          lineColour = [0.620 0.855 0.957]; %'b';
        elseif contains(areasOI(areaOrderOI(iComp)), {'Cx', 'SS', 'RSP', 'VIS'})
          lineColour = [0.843 0.788 0.361]; %'g';
        elseif contains(areasOI(areaOrderOI(iComp)), {'Hp', 'CA', 'DG'})
          lineColour = [0.867 0.561 0.891]; %'m';
        else
          lightGrey = 211;
          lineColour = [lightGrey lightGrey lightGrey]./256;
        end
        dataPointCount = dataPointCount + 1;
        x = [dataPointCount dataPointCount];
        y = phase(areaOrderOI(iComp)) + phaseCI(:,areaOrderOI(iComp));
        plot(x, y, '-', 'color',lineColour, 'LineWidth',0.75);
        if iComp == 1
          hold on;
        end
        plot(dataPointCount, phase(areaOrderOI(iComp)), 'o', ...
          'MarkerEdgeColor',lineColour, 'MarkerFaceColor',lineColour);
      end
    end
    hold off
    xticks(1:dataPointCount);
    xticklabels(areasOI_full(areaOrderOI(ismember(areaOrderOI, includeAreaInds))));
    xtickangle(45)
    if strcmpi(population, 'negative')
      yticks([0 pi/2 pi 3*pi/2]);
      yticklabels({'0', '\pi/2', '\pi', '3\pi/2'});
      ylim([0 3*pi/2]);
    elseif strcmpi(population, 'positive')
      yticks([-pi -pi/2 0 pi/4 pi/2 3*pi/4 pi 3*pi/2 2*pi]);
      yticklabels({'-\pi', '-\pi/2', '0', '\pi/4', '\pi/2' '3\pi/4', '\pi', '3\pi/2', '2\pi'});
    elseif strcmpi(population, 'all')
      yticks([-pi -pi/2 0 pi/4 pi/2 3*pi/4 pi 3*pi/2 2*pi]);
      yticklabels({'-\pi', '-\pi/2', '0', '\pi/4', '\pi/2' '3\pi/4', '\pi', '3\pi/2', '2\pi'});
    end
    yLim = ylim;
    ylabel('Phase (rad)')

    % Tidy the figure
    set(fH, 'Color', 'white');
    ax = gca;
    set(ax, 'box', 'off');
    set(ax, 'TickDir', 'out');
    ax.FontSize = fontSize - 4;
    %set(get(ax, 'XAxis'), 'FontWeight', 'bold');
    set(get(ax, 'YAxis'), 'FontWeight', 'bold');
    xlim([0.5 dataPointCount+0.5])
    if strcmpi(population, 'negative')
      ylim([0 3*pi/2]);
    elseif strcmpi(population, 'positive')
      ylim([pi/4 3*pi/4] - pi/16)
    elseif strcmpi(population, 'all')
      ylim([pi/4 pi] - pi/16)
    end
    ax.TickLength = [0.0025 0.0025];
  end

  % Perform stats
  for iFreq = 1:nFreqs
    stats{iFreq} = meanTest2(phaseFOI{iFreq}(areaOrderOI), 'wwtest');
    for iComp = 1:numel(stats{iFreq})
      stats{iFreq}(iComp).area1 = areasOI{areaOrderOI(stats{iFreq}(iComp).iCol1)};
      stats{iFreq}(iComp).area2 = areasOI{areaOrderOI(stats{iFreq}(iComp).iCol2)};
    end
    stats_np{iFreq} = meanTest2(phaseFOI{iFreq}(areaOrderOI), 'cmtest');
    for iComp = 1:numel(stats_np{iFreq})
      stats_np{iFreq}(iComp).area1 = areasOI{areaOrderOI(stats_np{iFreq}(iComp).iCol1)};
      stats_np{iFreq}(iComp).area2 = areasOI{areaOrderOI(stats_np{iFreq}(iComp).iCol2)};
    end
  end


  % Draw phase and positive fraction scatter
  % Go through area comparisons and compile them
  incMask = [false false false false false true true false true true true ...
    true true true false false false false false false false false false ...
    false false false false false false false false false false false false ...
    true true true true true true true true true true true true false false ...
    false true false];
  nComps = numel(areasOI);
  for iComp = 1:nComps

    % Get the data
    areaInds = find(getAreaInds(strrep(areasOI{iComp}, '_', '-'), infraslowAnalyses.areaSummaries.areaTable));
    compData = {};
    for iArea = 1:numel(areaInds)
      if strcmpi(population, 'all')
        if significantOnly
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohSignificant, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        else
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCoh, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        end
      elseif strcmpi(population, 'positive')
        if significantOnly
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohPositiveSignificant, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        else
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohPositive, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        end
      elseif strcmpi(population, 'negative')
        if significantOnly
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohNegativeSignificant, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
        else
          compData = [compData, getAreaCohData(infraslowAnalyses.spikingPupilCohNegative, infraslowAnalyses.areaSummaries.areaTable, areaInds(iArea))];
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
        phase = [phase; recData.phase(:,freqInds)]; %#ok<*SAGROW,*AGROW>
      end
    end
    disp([areasOI{iComp} ' ' num2str(numel(phase))]);
    if iComp == 1
      phaseProfiles = nan(nComps, size(phase,2));
      freq = recData.frequency(1, freqInds);
    end
    if ~isempty(phase)
      phaseProfiles(iComp,:) = datamean(phase, 'circularNP');
    end
  end

  % Correlate phase and positive fraction
  nFreqs = sum(freqInds);
  positiveSpearmanFractionsMeans = repmat( ...
    infraslowAnalyses.spikingPupilCorr.areasOI.positiveSpearmanFractionsMeans(:,1), 1, nFreqs);
  [rSpearman, pvalSpearman] = corrMulti( ...
    phaseProfiles(incMask,:)', positiveSpearmanFractionsMeans(incMask,:)', 'circlinearnp');
  phaseProfiles(phaseProfiles < 0) = phaseProfiles(phaseProfiles < 0) + 2*pi;

  % Fit the data
  nFreqs = numel(freqOI);
  slope = cell(nFreqs,1);
  coefficients = cell(nFreqs,1);
  for iF = 1:nFreqs
    freqInd = ismember(freq, freqOI(iF));
    [~, slope{iF}, coefficients{iF}] = fitLine( ...
      positiveSpearmanFractionsMeans(incMask,freqInd), ...
      phaseProfiles(incMask,freqInd), type='linear-circular-fma');
  end

  % Draw figures
  nFreqs = numel(freqOI);
  fontSize = 18;
  for iF = 1:nFreqs
    freqInd = ismember(freq, freqOI(iF));
    fH = figure;
    plot(positiveSpearmanFractionsMeans(incMask,freqInd), phaseProfiles(incMask,freqInd), 'b.', 'MarkerSize',15);
    xlabel('Positive fraction');
    yticks([-pi -pi/2 0 pi/4 pi/2 3*pi/4 pi 3*pi/2 2*pi])
    yticklabels({'-\pi', '-\pi/2', '0', '\pi/4', '\pi/2' '3\pi/4', '\pi', '3\pi/2', '2\pi'});
    ylabel('Phase (rad)')

    % Draw the line fit
    xLim = xlim;
    yFit = coefficients{iF}(1)*xLim + coefficients{iF}(2) ;%+ 2*pi;
    hold on; plot(xLim, yFit, ':k'); hold off

    % Add text
    yLim = ylim;
    xAxisSize = diff(xLim);
    yAxisSize = diff(yLim);
    xPos = xLim(1) + 0.05*xAxisSize;
    yPos = yLim(1) + 0.95*yAxisSize;
    txtStr = ['rho = ' num2str(round(rSpearman(freqInd), 1,'decimals')) ...
      '  p = ' num2str(round(pvalSpearman(freqInd), 1,'significant'))];
    text(xPos, yPos, txtStr);
    xlim(xLim);
    ylim(yLim);

    % Tidy the figure
    set(fH, 'Color', 'white');
    ax = gca;
    set(ax, 'box', 'off');
    set(ax, 'TickDir', 'out');
    ax.FontSize = fontSize - 4;
    %set(get(ax, 'XAxis'), 'FontWeight', 'bold');
    set(get(ax, 'YAxis'), 'FontWeight', 'bold');
  end
end