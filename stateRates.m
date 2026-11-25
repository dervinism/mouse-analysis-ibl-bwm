function [rates, logRates] = stateRates(counts, pupil, samplingRate, percentiles)
% [rates, logRates] = stateRates(counts, pupil, samplingRate, percentiles)
%
% Firing rates of units in the most constricted and most dilated pupil
% states, computed as in the UoL/Allen analysis (mouse-analysis
% globalFRPhaseAnalysis2.m / getLog.m): for each centile c the constricted
% state is the time with pupil area <= its c-th percentile and the dilated
% state the time with pupil area > its (100-c)-th percentile; a unit's rate
% in a state is its spike count in that state divided by the TIME SPENT IN
% THAT STATE (not by the whole recording). Units keep a rate in both states
% (zero when silent), so constricted vs dilated comparisons are paired;
% log10 of a zero rate is set to -3, as getLog.m does for UoL/Allen.
%
% Args:
%   counts (numeric, required, positional): shape-(nUnits, nSamples) spike
%     counts, binned at samplingRate.
%   pupil (numeric, required, positional): shape-(1, nSamples) pupil area.
%   samplingRate (numeric, required, positional): bins per second.
%   percentiles (numeric, required, positional): shape-(1, nCentiles)
%     centiles c (e.g. 12.5 for the most constricted / dilated 12.5%).
%
% Returns:
%   rates (numeric): shape-(nUnits, nCentiles, 2) rates in spikes/s,
%     (:,:,1) constricted, (:,:,2) dilated.
%   logRates (numeric): log10(rates) with log10(0) set to -3.

nUnits = size(counts, 1); nCentiles = numel(percentiles);
rates = zeros(nUnits, nCentiles, 2);
for iCent = 1:nCentiles
  lo = pupil <= prctile(pupil, percentiles(iCent));
  hi = pupil > prctile(pupil, 100 - percentiles(iCent));
  rates(:, iCent, 1) = sum(counts(:, lo), 2) / (sum(lo) / samplingRate);
  rates(:, iCent, 2) = sum(counts(:, hi), 2) / (sum(hi) / samplingRate);
end
logRates = log10(rates);
logRates(isinf(logRates)) = -3;
end
