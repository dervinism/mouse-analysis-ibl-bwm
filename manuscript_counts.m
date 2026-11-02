function manuscript_counts()
% IBL dataset counts for the Methods/Results: sessions, sessions with pupil, probe insertions, mice and units,
% overall and for the sessions with pupil data (the ones that enter the analyses), plus the units with a
% pupil correlation (classified). Writes tables/manuscript_counts.txt.
cd(fileparts(mfilename('fullpath')));
qualityFilter = false; %#ok<NASGU>
params
load(fullfile(processedDataFolder, 'bwmPreprocessedData2.mat'), 'infraslowData');
E = infraslowData.experimentData;
has = @(e, f) isfield(e, f) && ~isempty(e.(f));   % some sessions lack fields
hasP = cellfun(@(e) has(e, 'leftPupilAreaSize') || has(e, 'rightPupilAreaSize'), E);
subj = strings(size(E)); nU = zeros(size(E)); ins = zeros(size(E));
for i = 1:numel(E)
  e = E{i};
  if has(e, 'subjectID'), subj(i) = string(e.subjectID); end
  if has(e, 'spikeCounts'), nU(i) = size(e.spikeCounts, 1); end
  if has(e, 'probeIDs'), ins(i) = numel(unique(string(e.probeIDs))); end
end
noField = sum(~cellfun(@(e) isfield(e, 'leftPupilAreaSize'), E));
fid = fopen(fullfile('tables', 'manuscript_counts.txt'), 'w');
for f = [1 fid]
  fprintf(f, 'Sessions without pupil fields: %d; without subject ID: %d\n', noField, sum(subj == ""));
  fprintf(f, 'All sessions: %d sessions, %d insertions, %d mice, %d units\n', numel(E), sum(ins), numel(unique(subj)), sum(nU));
  fprintf(f, 'With pupil:   %d sessions, %d insertions, %d mice, %d units\n', sum(hasP), sum(ins(hasP)), ...
    numel(unique(subj(hasP))), sum(nU(hasP)));
end
A = load(fullfile(processedDataFolder, 'bwmAnalysisResults.mat'), 'infraslowAnalyses');
r = A.infraslowAnalyses.spikingPupilCorr.recordings.rSpearman;
if iscell(r), n = sum(cellfun(@(x) sum(isfinite(x(:))), r)); else, n = sum(isfinite(r(:))); end
for f = [1 fid], fprintf(f, 'Units with a pupil correlation (classified): %d\n', n); end
fclose(fid);
end
