function areaData = getAreaCohData(cohStruct, areaTable, areaInd)
% areaData = getAreaCohData(cohStruct, areaTable, areaInd)
%
% Returns the coherence results of one brain area, looked up by its acronym.
% Fields of the coherence structures (e.g., infraslowAnalyses.spikingPupilCoh)
% are named after area acronyms with spaces and hyphens replaced by
% underscores (see coherence_pupil.m). Areas without results (e.g., the
% unassigned '???' areas) return an empty cell.
%
% Arguments:
%   cohStruct (struct, required, positional): coherence results structure
%     with one field per brain area.
%   areaTable (table, required, positional): the brain area table
%     infraslowAnalyses.areaSummaries.areaTable.
%   areaInd (numeric, required, positional): row index of the area in
%     areaTable.
%
% Returns:
%   areaData (cell): the area's per-recording coherence results.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

fieldName = strrep(strrep(areaTable.Brain_area_acronym{areaInd}, ' ', '_'), '-', '_');
if isfield(cohStruct, fieldName)
  areaData = cohStruct.(fieldName);
else
  areaData = {};
end
end
