function manuscript_figures(outDir, parts)
% manuscript_figures(outDir)
%
% IBL panels of the manuscript figures, redrawn from the saved analysis
% results after the probe01 fix, with the geometry, fonts and export
% settings of the Figures_v3 panels so they drop into the same layouts.
% Each panel is saved as .fig and .png under outDir\Figure<N>\:
%   Figure3\FigureC_IBL            per-recording positive fraction per area
%                                  (correlate_pupil.m 'lines'), with the
%                                  t-test brackets of the Results drawn in
%   Figure5\Firing_rates_DG,
%           Split_firing_rates_DG  copied from figures\firing_rate_distros
%   Figure5\firing_rate_change_collection, positive_fraction_v_firing_rate_change
%                                  panels D and E (firing_rate_figures.m)
%   Figure8\Figure8_mean_<pooled|separate>
%                                  panels A-C (phase_profiles_pupil.m
%                                  summaries), brackets = Table 7 tests,
%                                  Holm-corrected within the table version
%   Figure8\Figure8_mean_<all|positive|negative>_DEF
%                                  panels D-F (signed circular-linear rho)
%
% Args:
%   outDir (char, optional): default Work\Leicester\Paper_1\Figures_v4.
%   parts (string array, optional): figures to draw, default ["3" "5" "8"].

if nargin < 2 || isempty(parts), parts = ["3" "5" "8"]; end
if nargin < 1 || isempty(outDir)
  outDir = 'C:\Users\m329786\Work\Leicester\Paper_1\Figures_v4';
end
cd(fileparts(mfilename('fullpath')));
qualityFilter = false; %#ok<NASGU> read by params
params
load(fullfile(processedDataFolder, ['bwmAnalysisResults' resultsTag '.mat']), 'infraslowAnalyses');
A = infraslowAnalyses; areaTable = A.areaSummaries.areaTable; corrOI = A.spikingPupilCorr.areasOI;
for d = ["Figure3" "Figure5" "Figure8"]
  if ~isfolder(fullfile(outDir, d)), mkdir(fullfile(outDir, d)); end
end
TH = [0.620 0.855 0.957]; CX = [0.843 0.788 0.361]; HP = [0.867 0.561 0.891]; GR = [211 211 211]/256;
includeAreaInds = [1:14 17:18, 24 29:32 36:47 51];      % the 34 area groups of Fig. 3C / Fig. 8A-C (as the scripts)
[sorted, order] = sort(corrOI.positiveSpearmanFractionsMeans(:,1), 'descend');
order = order(~isnan(sorted)); order = order(ismember(order, includeAreaInds));
labels = areasOI_full(order); keys = areasOI(order);
xOf = @(k) find(strcmp(keys, k));

%% Figure 3C
if any(parts == "3")
% axes 31.02 cm wide with the width:height ratio of the v3 panel as placed in Figure3.pptx (3.553), 1.2 cm headroom
axH = 31.02 / 3.553; figH = 4.074 + axH + 1.6;
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [1 1 40.03 figH]);
ax = axes(fig, 'Position', [0.13 4.074/figH 0.775 axH/figH]); hold(ax, 'on');
perRec = cell(numel(order), 1);
for i = 1:numel(order)
  x = corrOI.positiveSpearmanFractionsPerRec{order(i)}(:,1); perRec{i} = x;
  [m, ci] = datamean(x);
  c = local_colour(keys{i}, TH, CX, HP, GR);
  plot(ax, [i i], [min(x) max(x)], ':', 'Color', c, 'LineWidth', 0.75);
  plot(ax, [i i], [max(min(x), m + ci(1)) min(max(x), m + ci(2))], 'Color', c, 'LineWidth', 3.5);
  plot(ax, i, m, '.', 'Color', 'k', 'MarkerSize', 5);
end
xticks(ax, 1:numel(order)); xticklabels(ax, strcat(labels, '+')); xtickangle(ax, 60);
ylabel(ax, 'Fraction'); set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 14, 'XLim', [0.5 numel(order) + 0.5], 'YLim', [0.35 0.85], ...
  'YTick', [0.4 0.6 0.8]);   % ticks as the hand-placed labels of Figure3.pptx
ax.YAxis.FontWeight = 'bold';
% Results t-tests (recordings as units): Th vs nCx, nCx vs Hp (level 2); sensory vs motor + association (level 1)
p1 = local_t(perRec{xOf('Th')}, perRec{xOf('nCx')}); p2 = local_t(perRec{xOf('nCx')}, perRec{xOf('Hp')});
p3 = local_t(perRec{xOf('sensory')}, [perRec{xOf('motor')}; perRec{xOf('association')}]);
fprintf('Fig 3C t-tests: Th-nCx p=%.3g, nCx-Hp p=%.3g, sensory-(motor+assoc) p=%.3g\n', p1, p2, p3);
yl = 0.85; R = 0.5; lw = 3.5; fs = 26;   % as the hand-drawn brackets of panels A and B once placed in the layout
local_join(ax, xOf('sensory'), [xOf('motor') xOf('association')], yl + 0.015*R, yl + 0.035*R, 0.015*R, local_stars(p3), lw, fs);
local_bracket(ax, xOf('Th'), xOf('nCx'), yl + 0.09*R, 0.015*R, local_stars(p1), lw, fs);
local_bracket(ax, xOf('nCx'), xOf('Hp'), yl + 0.09*R, 0.015*R, local_stars(p2), lw, fs);
local_save(fig, fullfile(outDir, 'Figure3', 'FigureC_IBL'), [40.03 figH], 100);
end

%% Figure 8A-C, both versions, and D-F
if any(parts == "8")
T7 = readtable(fullfile('tables', 'phase_ordering_stats_qc0.csv'), 'TextType', 'string');
pops = {'all', 'spikingPupilCoh', 'all units'; 'positive', 'spikingPupilCohPositive', 'positive units'; ...
  'negative', 'spikingPupilCohNegative', 'negative units'};
ph = cell(3, 1); phCI = cell(3, 1);
for iPop = 1:3   % per-recording circular mean phase at 0.3 Hz, then circular mean +- CI over recordings (as phase_profiles_pupil.m)
  coh = A.(pops{iPop, 2}); ph{iPop} = nan(numel(order), 1); phCI{iPop} = nan(numel(order), 2);
  for i = 1:numel(order)
    inds = find(getAreaInds(keys{i}, areaTable)); recPh = [];
    for k = 1:numel(inds)
      d = getAreaCohData(coh, areaTable, inds(k));
      for r = 1:numel(d)
        if isstruct(d{r}) && isfield(d{r}, 'fullInterpCoherence')
          fi = abs(d{r}.fullInterpCoherence.frequency(1,:) - 0.3) < 1e-6;
          recPh = [recPh; datamean(d{r}.fullInterpCoherence.phase(:, fi), 'circularNP')]; %#ok<AGROW>
        end
      end
    end
    [m, ci] = datamean(recPh, 'circularNP');
    if m < 0, m = m + 2*pi; end
    ph{iPop}(i) = m; phCI{iPop}(i, :) = ci(:)';
  end
end
families = {'pooled', {'Hp', 'nCx'; 'nCx', 'Th'; 'association-motor', 'sensory'; 'association', 'motor'}; ...
  'separate', {'Hp', 'nCx'; 'nCx', 'Th'; 'association', 'sensory'; 'motor', 'sensory'; 'association', 'motor'}};
for iv = 1:2
  comps = families{iv, 2};
  sel = false(height(T7), 1);
  for c = 1:size(comps, 1), sel = sel | (T7.group1 == comps{c, 1} & T7.group2 == comps{c, 2}); end
  pH = nan(height(T7), 1); pH(sel) = local_holm(T7.p(sel));
  P = @(g1, g2, pop) pH(T7.group1 == g1 & T7.group2 == g2 & T7.population == pop);
  fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [1 1 40.26 23.14]);
  axPos = {[0.1 0.687 0.87 0.213], [0.1 0.423 0.87 0.213], [0.1 0.16 0.87 0.213]};
  for iPop = 1:3
    ax = axes(fig, 'Position', axPos{iPop}); hold(ax, 'on');
    y = ph{iPop}; ci = phCI{iPop};
    for i = 1:numel(order)
      c = local_colour(keys{i}, TH, CX, HP, GR);
      plot(ax, [i i], y(i) + ci(i, :), '-', 'Color', c, 'LineWidth', 2.5);
      plot(ax, i, y(i), 'o', 'MarkerEdgeColor', c, 'MarkerFaceColor', c, 'MarkerSize', 7);
    end
    lo = min(y + ci(:, 1)); hi = max(y + ci(:, 2)); R = hi - lo; pop = pops{iPop, 1};
    lv = @(n) hi + (0.10 + 0.13 * (n - 1)) * R; tk = 0.035 * R; lw = 1.3; fs = 15;
    local_bracket(ax, xOf('Th'), xOf('nCx'), lv(1), tk, local_stars(P('nCx', 'Th', pop)), lw, fs);
    local_bracket(ax, xOf('nCx'), xOf('Hp'), lv(1), tk, local_stars(P('Hp', 'nCx', pop)), lw, fs);
    if iv == 1
      local_join(ax, xOf('sensory'), [xOf('motor') xOf('association')], lv(2), lv(3), tk, ...
        local_stars(P('association-motor', 'sensory', pop)), lw, fs);
      local_bracket(ax, xOf('motor'), xOf('association'), lv(4), tk, local_stars(P('association', 'motor', pop)), lw, fs);
      top = lv(4);
    else
      local_bracket(ax, xOf('motor'), xOf('association'), lv(2), tk, local_stars(P('association', 'motor', pop)), lw, fs);
      local_bracket(ax, xOf('sensory'), xOf('motor'), lv(3), tk, local_stars(P('motor', 'sensory', pop)), lw, fs);
      local_bracket(ax, xOf('sensory'), xOf('association'), lv(4), tk, local_stars(P('association', 'sensory', pop)), lw, fs);
      top = lv(4);
    end
    yl = [lo - 0.08 * R, top + 0.16 * R];
    tks = (ceil(yl(1) / (pi/4)):floor((hi + 0.05 * R) / (pi/4))) * pi/4;   % ticks over the data only, not the bracket area
    set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 20, 'XLim', [0.5 numel(order) + 0.5], 'YLim', yl, ...
      'YTick', tks, 'YTickLabel', arrayfun(@local_pilabel, tks / (pi/4), 'UniformOutput', false), 'TickLength', [0.0025 0.0025]);
    title(ax, sprintf('%s (mean, 0.3 Hz)', pops{iPop, 3}), 'FontSize', 15); ylabel(ax, 'Phase (rad)');
    xticks(ax, 1:numel(order));
    if iPop == 3, xticklabels(ax, labels); xtickangle(ax, 60); else, xticklabels(ax, {}); end
  end
  local_save(fig, fullfile(outDir, 'Figure8', ['Figure8_mean_' families{iv, 1}]), [43.18 44.45], 300);
end
% D-F: unit-level circular mean per area group (as the scatter part of phase_profiles_pupil.m), 21 groups
groups = {'paCx', 'pCx', 'sensory-Th', 'association-Th', 'motor-Th', 'sensory-nCx', 'association-nCx', 'motor-nCx', ...
  'CA', 'DG', 'AA', 'BG', 'Hyp', 'MBr', 'Cereb', 'CLA', 'BS', 'EPT', 'BF', 'STh', 'piCx'};
gi = cellfun(@(g) find(strcmp(areasOI, g)), groups);
frac = corrOI.positiveSpearmanFractionsMeans(gi, 1);
for iPop = 1:3
  coh = A.(pops{iPop, 2}); y = nan(numel(gi), 1);
  for j = 1:numel(gi)
    inds = find(getAreaInds(areasOI{gi(j)}, areaTable)); phase = [];
    for k = 1:numel(inds)
      d = getAreaCohData(coh, areaTable, inds(k));
      for r = 1:numel(d)
        if isstruct(d{r}) && isfield(d{r}, 'fullInterpCoherence')
          fi = abs(d{r}.fullInterpCoherence.frequency(1,:) - 0.3) < 1e-6;
          phase = [phase; d{r}.fullInterpCoherence.phase(:, fi)]; %#ok<AGROW>
        end
      end
    end
    y(j) = datamean(phase, 'circularNP');
  end
  [rho, pv] = corrMulti(y', frac', 'circlinearnp');
  y(y < 0) = y(y < 0) + 2*pi;
  [~, slope, coef] = fitLine(frac, y, type='linear-circular-fma');
  fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [1 1 21.17 17.2]);
  ax = axes(fig, 'Position', [0.135 0.133 0.77 0.792]); hold(ax, 'on');
  plot(ax, frac, y, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'LineWidth', 0.5);
  xl = [floor(min(frac) * 20) / 20, ceil(max(frac) * 20) / 20];
  k = round(mean(y - (coef(1) * frac + coef(2))) / (2*pi));   % circular fit: put the line on the same cycle as the points
  plot(ax, xl, coef(1) * xl + coef(2) + 2*pi*k, ':k', 'LineWidth', 1.2);
  R = max(y) - min(y); yl = [min(y) - 0.08 * R, max(y) + 0.18 * R];
  tks = (ceil(yl(1) / (pi/8)):floor(yl(2) / (pi/8))) * pi/8;
  set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 20, 'XLim', xl, 'YLim', yl, 'XTick', xl(1):0.1:xl(2), ...
    'YTick', tks, 'YTickLabel', arrayfun(@local_pilabel8, round(tks / (pi/8)), 'UniformOutput', false));
  xlabel(ax, 'Positive fraction'); if iPop == 1, ylabel(ax, 'Phase (rad)'); end
  text(ax, xl(2), yl(2), sprintf('\\rho = %s%.2f  p = %.2g', char(8722 * (slope < 0) + 43 * (slope >= 0)), rho, pv), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 20, 'FontWeight', 'bold');
  panelLetters = 'DEF';
  fprintf('Fig 8%s (%s): rho = %s%.2f, p = %.3g, %d groups\n', panelLetters(iPop), pops{iPop, 1}, char(45 * (slope < 0) + 43 * (slope >= 0)), rho, pv, numel(gi));
  local_save(fig, fullfile(outDir, 'Figure8', ['Figure8_mean_' pops{iPop, 1} '_DEF']), [21.17 17.2], 300);
end
end

%% Figure 5: DG distributions (already regenerated) and panels D, E
if ~any(parts == "5"), return, end
for f = ["Firing_rates_DG" "Split_firing_rates_DG"]
  for e = [".fig" ".png"]
    copyfile(fullfile(figFolder, 'firing_rate_distros', f + e), fullfile(outDir, 'Figure5', f + e));
  end
end
load(fullfile(processedDataFolder, 'bwmPreprocessedData2.mat'), 'infraslowData');
firingRateThr = 50/3600; nAreas = size(areaLabels, 1); firingRates = cell(nAreas, 1);
for iArea = 1:nAreas   % log10 rate in the most constricted / dilated 12.5% (as firing_rate_figures.m)
  areaInds = A.areaSummaries.groupedUnitInds{iArea}; fr = nan(size(areaInds, 1), 2);
  for iExp = unique(areaInds(:, 1))'
    unitMask = areaInds(:, 1) == iExp; e = infraslowData.experimentData{iExp};
    if ~isempty(e.leftPupilAreaSize) && ~isempty(e.rightPupilAreaSize), pupil = mean([e.leftPupilAreaSize; e.rightPupilAreaSize]);
    elseif ~isempty(e.leftPupilAreaSize), pupil = e.leftPupilAreaSize;
    elseif ~isempty(e.rightPupilAreaSize), pupil = e.rightPupilAreaSize;
    else, continue
    end
    counts = e.spikeCounts(areaInds(unitMask, 2), :); dur = size(counts, 2) / e.samplingRate;
    active = sum(counts, 2) / dur > firingRateThr; pv = prctile(pupil, [12.5 87.5]);
    lo = zeros(sum(unitMask), 1); hi = lo;
    lo(active) = sum(counts(active, pupil <= pv(1)), 2) / dur; hi(active) = sum(counts(active, pupil > pv(2)), 2) / dur;
    lo(~lo) = nan; hi(~hi) = nan; fr(unitMask, :) = [log10(lo) log10(hi)];
  end
  firingRates{iArea} = fr;
end
clear infraslowData
uolAreas = {'CA UoL', 'neoCx UoL', 'DG UoL', 'CA-DG UoL', 'RSP UoL', 'S1 UoL', 'Th UoL'};   % typed-in values of firing_rate_figures.m
uolChange = [0.145666 - 0.298024; 0.341518 - 0.259733; 0.289202 - 0.458551; 0.249652 - 0.424188; 0.434660 - 0.418259; ...
  0.300893 - 0.190590; 0.770454 - 0.415658];
uolFrac = [0.475467; 0.674299; 0.400523; 0.461996; 0.607605; 0.69995; 0.874251]; uolInc = [true false true false true true true];
allenAreas = {'CA Allen', 'DG Allen', 'CA-DG Allen', 'Th Allen', 'V1 Allen', 'V2 Allen', 'Cx Allen'};
allenChange = [0.675289 - 0.635800; 0.539034 - 0.615952; 0.647850 - 0.631803; 0.933029 - 0.693471; 0.452149 - 0.481036; ...
  0.476366 - 0.500610; 0.471062 - 0.496322];
allenFrac = [0.400404; 0.235446; 0.384781; 0.781899; 0.446529; 0.473333; 0.467418]; allenInc = [true true false true true true false];
rateChange = nan(numel(areasOI), 1);
for i = 1:numel(areasOI)
  x = concatenateCells(firingRates(getAreaInds(areasOI{i}, areaTable)));
  if ~isempty(x), rateChange(i) = mean(x(:, 2), 'omitnan') - mean(x(:, 1), 'omitnan'); end
end
keysAll = [areasOI; uolAreas'; allenAreas']; labAll = [areasOI_full; uolAreas'; allenAreas'];
labAll(strcmp(keysAll, 'CA-DG UoL')) = {'Hp UoL'}; labAll(strcmp(keysAll, 'CA-DG Allen')) = {'Hp Allen'};
fracAll = [corrOI.positiveSpearmanFractionsMeans(:, 1); uolFrac; allenFrac]; changeAll = [rateChange; uolChange; allenChange];
incD = [includeAreaInds, numel(areasOI) + (1:14)];                 % panel D: the 34 IBL groups + all UoL / Allen entries
[s, o] = sort(fracAll, 'descend'); o = o(~isnan(s)); o = o(ismember(o, incD));
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [1 1 40.47 11.11]);
ax = axes(fig, 'Position', [0.13 0.284 0.775 0.641]); hold(ax, 'on');
for i = 1:numel(o)
  plot(ax, [i i], [0 changeAll(o(i))], '-', 'Color', local_colour(keysAll{o(i)}, TH, CX, HP, GR), 'LineWidth', 1.5);
end
plot(ax, [0.5 numel(o) + 0.5], [0 0], ':');
yl = max(0.4, ceil(max(abs(changeAll(o))) * 10) / 10);
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 14, 'XLim', [0.5 numel(o) + 0.5], 'YLim', [-yl yl], 'XTick', 1:numel(o), ...
  'XTickLabel', labAll(o)); xtickangle(ax, 45); ylabel(ax, 'Log_{10}(Firing rate) change');
local_save(fig, fullfile(outDir, 'Figure5', 'firing_rate_change_collection'), [40.47 11.11], 100);
incE = [ismember(1:numel(areasOI), cellfun(@(g) find(strcmp(areasOI, g)), groups)), uolInc, allenInc];   % panel E: 21 + 5 + 5
[rhoE, pE] = corrMulti(fracAll(incE)', changeAll(incE)', 'Spearman');
[~, ~, cE] = fitLine(fracAll(incE), changeAll(incE), type='linear-linear');
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [1 1 14.82 11.11]);
ax = axes(fig, 'Position', [0.146 0.146 0.759 0.779]); hold(ax, 'on');
plot(ax, fracAll(incE), changeAll(incE), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
xl = [0.2 0.9]; plot(ax, xl, cE(1) * xl + cE(2), ':k');
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 14, 'XLim', xl, 'YLim', [-yl yl]);
xlabel(ax, 'Positive fraction'); ylabel(ax, 'Log_{10}(Firing rate) change');
text(ax, xl(1) + 0.05 * diff(xl), 0.95 * yl, sprintf('rho = %.2f  p = %.2g', rhoE, pE), 'FontSize', 10);
fprintf('Fig 5E: rho = %.2f, p = %.3g, %d groups\n', rhoE, pE, sum(incE));
local_save(fig, fullfile(outDir, 'Figure5', 'positive_fraction_v_firing_rate_change'), [14.82 11.11], 100);
end


function c = local_colour(key, TH, CX, HP, GR)
% thalamus / cortex / hippocampus / other, as in the figure scripts
if contains(key, {'Th', 'VP', 'PO', 'LG', 'LP'}) && ~strcmpi(key, 'STh'), c = TH;
elseif contains(key, {'Cx', 'SS', 'RSP', 'VIS', 'V2', 'V1', 'S1'}), c = CX;
elseif contains(key, {'Hp', 'CA', 'DG'}), c = HP;
else, c = GR;
end
end

function p = local_t(a, b)
[~, p] = ttest2(a, b);
end

function s = local_stars(p)
if p < 1e-4, s = '****'; elseif p < 1e-3, s = '***'; elseif p < 0.01, s = '**'; elseif p < 0.05, s = '*'; else, s = 'n.s.'; end
end

function adj = local_holm(p)
% Holm-Bonferroni adjusted p values (monotone)
[ps, o] = sort(p); m = numel(p); a = min(1, (m:-1:1)' .* ps(:)); a = cummax(a); adj = zeros(size(p)); adj(o) = a;
end

function local_bracket(ax, x1, x2, y, tk, label, lw, fs)
% horizontal bracket from x1 to x2 at height y with down ticks of length tk, label above the middle
plot(ax, [x1 x1 x2 x2], [y - tk y y y - tk], 'k-', 'LineWidth', lw, 'Clipping', 'off');
text(ax, mean([x1 x2]), y, label, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fs, ...
  'FontWeight', 'bold', 'Clipping', 'off');
end

function local_join(ax, x1, xs, ySub, y, tk, label, lw, fs)
% bracket from x1 to the middle of a sub-bracket joining the two positions xs (one group vs two pooled groups)
xs = sort(xs); xm = mean(xs);
plot(ax, [xs(1) xs(1) xs(2) xs(2)], [ySub - tk ySub ySub ySub - tk], 'k-', 'LineWidth', lw, 'Clipping', 'off');
plot(ax, [xm xm], [ySub y], 'k-', 'LineWidth', lw, 'Clipping', 'off');
local_bracket(ax, x1, xm, y, tk, label, lw, fs);
end

function s = local_pilabel(k)
% label of k*pi/4
switch k
  case 0, s = '0';
  otherwise
    g = gcd(abs(k), 4); n = k / g; d = 4 / g;
    if n == 1, num = '\pi'; elseif n == -1, num = '-\pi'; else, num = sprintf('%d\\pi', n); end
    if d == 1, s = num; else, s = sprintf('%s/%d', num, d); end
end
end

function s = local_pilabel8(k)
% label of k*pi/8
if k == 0, s = '0'; return, end
g = gcd(abs(k), 8); n = k / g; d = 8 / g;
if n == 1, num = '\pi'; elseif n == -1, num = '-\pi'; else, num = sprintf('%d\\pi', n); end
if d == 1, s = num; else, s = sprintf('%s/%d', num, d); end
end

function local_save(fig, base, sizeCm, dpi)
% .fig + .png at a fixed page size (no tight cropping, so the pixel geometry is known)
set(fig, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 sizeCm], 'PaperSize', sizeCm);
savefig(fig, [base '.fig']);
print(fig, [base '.png'], '-dpng', sprintf('-r%d', dpi));
close(fig);
end
