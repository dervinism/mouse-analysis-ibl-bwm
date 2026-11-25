% extract_uol_allen_fig5.m - per-unit pupil-state firing rates and pupil correlations of the UoL and Allen areas
% used in Figure 5 (distributions and rate change vs positive fraction), taken from the UoL/Allen analysis outputs
% (mouse-analysis globalFRPhaseAnalysis.m) instead of the values hard-coded in firing_rate_figures.m (2026-11-08).
% Sources (awake condition; data root = F: or E: \Datasets\infraslow-dynamics, whichever is mounted):
%   UoL:   <root>\04_data_analysis_3\noRun\rSpearmanBased_uol_mean\area2itself_comparisons_uol\units\globalUnits.mat
%   Allen: <root>\04_data_analysis\006_analysis_results\noRun\area2itself_comparisons_allensdk\units\globalUnits.mat
%          (Aug/Sep 2025 run, 24 sessions; user's choice 2026-11-08)
% Rates: areaFRPercentileIndividual_12p5{awake}{area}(:,[1 end]) = rate in the lowest / highest pupil octile (Hz).
% Positive fractions (positiveFractionTable, used by manuscript_figures.m / manuscript_numbers.m Fig. 5D-E):
%   UoL   = Table 2 per-recording means (stats file behind Figure 3A). 2026-11-10: the stats file's columns are in the
%           merged violin order Th, S1, Cx, RSC, Hp, CA, DG; the previous mapping assumed the published Table 2 row
%           order, which had CA and Hp swapped, so CA UoL and CA-DG UoL carried each other's fraction (fixed).
%   Allen = Spearman sign rule (user decision 2026-11-10): per-recording mean of the fraction of units with rho >= 0,
%           rho = arearSpearman5secIndividual of the Aug/Sep 2025 run, series with >= 10 units with a rho (as the
%           pupilCorrFractionPlot.m unit threshold). Was the phase-based split of the Feb-2022 allData_allensdk.mat
%           (published Table 3 up to MS1_v17); that value is kept in positiveFractionPhaseTable.
% Writes dependencies\uolAllen_fig5.mat.
cd(fileparts(mfilename('fullpath')));
roots = {'F:\Datasets\infraslow-dynamics', 'E:\Datasets\infraslow-dynamics'};
root = roots{find(cellfun(@isfolder, roots), 1)};
src = struct( ...
  'uol', fullfile(root, '04_data_analysis_3\noRun\rSpearmanBased_uol_mean\area2itself_comparisons_uol\units\globalUnits.mat'), ...
  'allen', fullfile(root, '04_data_analysis\006_analysis_results\noRun\area2itself_comparisons_allensdk\units\globalUnits.mat'));
unitThr = 10;   % minimum number of units with a rho per recording series (pupilCorrFractionPlot.m)
sets = {'uol', {'CA', 'DG', 'RSC', 'S1', 'Th', 'Cx', 'Hp'}, {'CA UoL', 'DG UoL', 'RSP UoL', 'S1 UoL', 'Th UoL', 'neoCx UoL', 'CA-DG UoL'};
        'allen', {'CA', 'DG', 'Th', 'V1', 'V2', 'VIS', 'Hp'}, {'CA Allen', 'DG Allen', 'Th Allen', 'V1 Allen', 'V2 Allen', 'Cx Allen', 'CA-DG Allen'}};
D = struct('label', {}, 'dataset', {}, 'area', {}, 'rates', {}, 'logRates', {}, 'rho', {}, 'recID', {}, 'source', {}, ...
  'positiveFractionPooled', {}, 'positiveFractionPerRec', {}, 'nRecs', {}, 'positiveFractionSpearmanThr', {}, ...
  'perRecFractionsSpearmanThr', {}, 'nRecsSpearmanThr', {});
for s = 1:size(sets, 1)
  f = src.(sets{s, 1});
  S = load(f, 'areas', 'conditions', 'areaFRPercentileIndividual_12p5', 'arearSpearman5secIndividual', ...
    'areaRecIDSpearman5secIndividual');
  iC = find(strcmp(S.conditions, 'awake'));
  for a = 1:numel(sets{s, 2})
    iA = find(strcmp(S.areas, sets{s, 2}{a}));
    X = S.areaFRPercentileIndividual_12p5{iC}{iA};
    r = S.arearSpearman5secIndividual{iC}{iA}; r = r(:);
    rec = S.areaRecIDSpearman5secIndividual{iC}{iA};
    if size(X, 1) ~= numel(r)
      error('extract:rows', '%s %s: %d rate rows but %d correlations', sets{s, 1}, sets{s, 2}{a}, size(X, 1), numel(r));
    end
    % keep only units with a pupil correlation (classified units, as counted in the Results and Tables 2/3; the Aug/Sep
    % 2025 Allen file also holds state rates for units without rho, e.g. 5 whole sessions) - 2026-11-08
    keepU = ~isnan(r);
    X = X(keepU, :); r = r(keepU); rec = rec(keepU);
    rates = X(:, [1 end]);
    L = log10(rates); L(isinf(L)) = -3;            % silent in a state -> -3, as in getLog.m
    if iscell(rec), rec = string(rec); end
    rec = rec(:);
    ok = ~isnan(r);
    ur = unique(rec(ok)); fr = nan(numel(ur), 1); nu = zeros(numel(ur), 1);
    for k = 1:numel(ur), m = ok & rec == ur(k); fr(k) = mean(r(m) >= 0); nu(k) = sum(m); end
    frThr = fr(nu >= unitThr);
    D(end + 1) = struct('label', sets{s, 3}{a}, 'dataset', sets{s, 1}, 'area', sets{s, 2}{a}, 'rates', rates, ...
      'logRates', L, 'rho', r, 'recID', rec, 'source', f, 'positiveFractionPooled', mean(r(ok) >= 0), ...
      'positiveFractionPerRec', mean(fr), 'nRecs', numel(ur), 'positiveFractionSpearmanThr', mean(frThr), ...
      'perRecFractionsSpearmanThr', frThr, 'nRecsSpearmanThr', numel(frThr)); %#ok<SAGROW>
    fprintf('%-12s n %5d (NaN rho %d), recs %2d | mean %.6f -> %.6f (change %.6f) | median %.4f -> %.4f | pos frac pooled %.6f per-rec %.6f | per-rec (>=%d units) %.6f (%d recs)\n', ...
      sets{s, 3}{a}, size(X, 1), sum(~ok), numel(ur), mean(L(:, 1), 'omitnan'), mean(L(:, 2), 'omitnan'), ...
      mean(L(:, 2), 'omitnan') - mean(L(:, 1), 'omitnan'), median(L(:, 1), 'omitnan'), median(L(:, 2), 'omitnan'), ...
      mean(r(ok) >= 0), mean(fr), unitThr, mean(frThr), numel(frThr));
  end
end
% Positive fractions of the stats files behind Figure 3A/B (per-recording means; merged violin column order).
tabSrc = struct('uol', 'C:\Users\m329786\Work\Leicester\Paper_1\Figures_v1\Figure3\statsFractionsUnitsPositiveMergedUOL.mat', ...
  'allen', 'C:\Users\m329786\Work\Leicester\Paper_1\Figures_v1\Figure3\statsFractionsUnitsPositiveMergedAllensdk.mat');
tabRows = struct('uol', {{'Th UoL', 'S1 UoL', 'neoCx UoL', 'RSP UoL', 'CA-DG UoL', 'CA UoL', 'DG UoL'}}, ...
  'allen', {{'Th Allen', 'V2 Allen', 'Cx Allen', 'V1 Allen', 'CA Allen', 'CA-DG Allen', 'DG Allen'}});
tabVals = struct('uol', [0.86 0.70 0.67 0.55 0.46 0.45 0.41], 'allen', [0.78 0.47 0.47 0.45 0.40 0.38 0.24]);
[D.positiveFractionTable] = deal(NaN); [D.positiveFractionTableSource] = deal(''); [D.positiveFractionPhaseTable] = deal(NaN);
for ds = {'uol', 'allen'}
  T = load(tabSrc.(ds{1})); st = T.statsFractionsUnitsPositiveMerged;
  colMean = nan(1, 7);
  for k = 1:numel(st), colMean(st(k).iCol1) = st(k).mean1; colMean(st(k).iCol2) = st(k).mean2; end
  assert(all(abs(round(colMean, 2) - tabVals.(ds{1})) < 1e-9), 'table check failed for %s: %s', ds{1}, mat2str(colMean, 4));
  if strcmp(ds{1}, 'uol')   % Hp (lHp1 + lCA3 series) 0.4576 and CA (lCA1 + lCA3) 0.4470 in area_properties_uol\all\pupilCorrFractions.mat
    assert(abs(colMean(5) - 0.4576) < 5e-4 && abs(colMean(6) - 0.4470) < 5e-4, 'UoL Hp/CA column check failed');
  end
  for k = 1:7
    i = find(strcmp({D.label}, tabRows.(ds{1}){k}));
    if strcmp(ds{1}, 'uol')
      D(i).positiveFractionTable = colMean(k); D(i).positiveFractionTableSource = tabSrc.(ds{1});
    else
      D(i).positiveFractionPhaseTable = colMean(k);
      D(i).positiveFractionTable = D(i).positiveFractionSpearmanThr;
      D(i).positiveFractionTableSource = sprintf('Spearman sign, %s (series with >= %d units with rho)', src.allen, unitThr);
    end
    fprintf('%-12s table fraction %.6f%s\n', tabRows.(ds{1}){k}, D(i).positiveFractionTable, ...
      ifelse(strcmp(ds{1}, 'allen'), sprintf(' (phase-based %.6f)', colMean(k)), ''));
  end
end
uolAllenFig5 = D; %#ok<NASGU>
save(fullfile('dependencies', 'uolAllen_fig5.mat'), 'uolAllenFig5', 'src');
fprintf('saved dependencies\\uolAllen_fig5.mat\n');


function out = ifelse(c, a, b)
if c, out = a; else, out = b; end
end
