function saveAnalysisResults(analysisResultsFile, infraslowAnalyses)
% saveAnalysisResults(analysisResultsFile, infraslowAnalyses)
%
% Saves the analysis results structure so that an interrupted save cannot
% corrupt the results file: the structure is written to a temporary file
% first, which then replaces the results file.
%
% Arguments:
%   analysisResultsFile (char, required, positional): full path of the
%     analysis results file.
%   infraslowAnalyses (struct, required, positional): the analysis results.

tempFile = [analysisResultsFile '.tmp'];
save(tempFile, 'infraslowAnalyses', '-v7.3');
movefile(tempFile, analysisResultsFile, 'f');
end
