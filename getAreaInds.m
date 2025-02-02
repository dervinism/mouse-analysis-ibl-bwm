function areaInds = getAreaInds(areaOI, areaTable)
% areaInds = getAreaInds(areaOI, areaTable)
%
% Identifies indices of combined and non-combined brain areas of interest
% within the brain area table.
%
% Arguments:
%   areaOI (char, required, positional): a shape-(1, N) character array
%     containing the area name of interest. Areas of interest are defined
%     in the parameters.m script within the areasOI variable.
%   areaTable (table, required, positional): a shape-(M, 6) table
%     corresponding to the variable
%     infraslowAnalyses.areaSummaries.areaTable.
%
% Returns:
%   areaInds (logical): a shape-(244, 1) logical array of sub-area indices
%     comprising the area of interest.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  areaOI (1,:) %{mustBeVector,mustBeText}
  areaTable (244,6) %{mustBeA(areaTable,'table')}
end

if ismember(areaOI, areaTable.Brain_area_type)
  areaInds = ismember(areaTable.Brain_area_type, areaOI);
elseif ismember(areaOI, areaTable.Brain_area_group)
  areaInds = ismember(areaTable.Brain_area_group, areaOI);
elseif strcmpi(areaOI, 'sensory-Th')
  areaInds = ismember(areaTable.Brain_area_type, 'sensory') & ...
    ismember(areaTable.Brain_area_group, 'Th');
elseif strcmpi(areaOI, 'association-Th')
  areaInds = ismember(areaTable.Brain_area_type, 'association') & ...
    ismember(areaTable.Brain_area_group, 'Th');
elseif strcmpi(areaOI, 'motor-Th')
  areaInds = ismember(areaTable.Brain_area_type, 'motor') & ...
    ismember(areaTable.Brain_area_group, 'Th');
elseif strcmpi(areaOI, 'sensory-nCx')
  areaInds = ismember(areaTable.Brain_area_type, 'sensory') & ...
    ismember(areaTable.Brain_area_group, 'nCx');
elseif strcmpi(areaOI, 'association-nCx')
  areaInds = ismember(areaTable.Brain_area_type, 'association') & ...
    ismember(areaTable.Brain_area_group, 'nCx');
elseif strcmpi(areaOI, 'motor-nCx')
  areaInds = ismember(areaTable.Brain_area_type, 'motor') & ...
    ismember(areaTable.Brain_area_group, 'nCx');
elseif strcmpi(areaOI, 'association-paCx')
  areaInds = ismember(areaTable.Brain_area_type, 'association') & ...
    ismember(areaTable.Brain_area_group, 'paCx');
elseif strcmpi(areaOI, 'sensory-pCx')
  areaInds = ismember(areaTable.Brain_area_type, 'sensory') & ...
    ismember(areaTable.Brain_area_group, 'pCx');
elseif strcmpi(areaOI, 'association-pCx')
  areaInds = ismember(areaTable.Brain_area_type, 'association') & ...
    ismember(areaTable.Brain_area_group, 'pCx');
elseif strcmpi(areaOI, 'VPL-VPM')
  areaInds = ismember(areaTable.Brain_area_acronym, 'VPL') | ...
    ismember(areaTable.Brain_area_acronym, 'VPM');
elseif strcmpi(areaOI, 'LG')
  areaInds = ismember(areaTable.Brain_area_acronym, 'LGd') | ...
    ismember(areaTable.Brain_area_acronym, 'LGv');
elseif strcmpi(areaOI, 'VPL-VPM-PO')
  areaInds = ismember(areaTable.Brain_area_acronym, 'VPL') | ...
    ismember(areaTable.Brain_area_acronym, 'VPM')| ...
    ismember(areaTable.Brain_area_acronym, 'PO');
elseif strcmpi(areaOI, 'VPL-VPM-PO-LP')
  areaInds = ismember(areaTable.Brain_area_acronym, 'VPL') | ...
    ismember(areaTable.Brain_area_acronym, 'VPM') | ...
    ismember(areaTable.Brain_area_acronym, 'PO') | ...
    ismember(areaTable.Brain_area_acronym, 'LP');
elseif strcmpi(areaOI, 'VPL-VPM-LG')
  areaInds = ismember(areaTable.Brain_area_acronym, 'VPL') | ...
    ismember(areaTable.Brain_area_acronym, 'VPM') | ...
    ismember(areaTable.Brain_area_acronym, 'LGd') | ...
    ismember(areaTable.Brain_area_acronym, 'LGv');
elseif strcmpi(areaOI, 'VPL-VPM-LG-PO')
  areaInds = ismember(areaTable.Brain_area_acronym, 'VPL') | ...
    ismember(areaTable.Brain_area_acronym, 'VPM') | ...
    ismember(areaTable.Brain_area_acronym, 'LGd') | ...
    ismember(areaTable.Brain_area_acronym, 'LGv') | ...
    ismember(areaTable.Brain_area_acronym, 'PO');
elseif strcmpi(areaOI, 'VPL-VPM-LG-PO-LP')
  areaInds = ismember(areaTable.Brain_area_acronym, 'VPL') | ...
    ismember(areaTable.Brain_area_acronym, 'VPM') | ...
    ismember(areaTable.Brain_area_acronym, 'LGd') | ...
    ismember(areaTable.Brain_area_acronym, 'LGv') | ...
    ismember(areaTable.Brain_area_acronym, 'PO') | ...
    ismember(areaTable.Brain_area_acronym, 'LP');
elseif strcmpi(areaOI, 'SSp-body')
  areaInds = ismember(areaTable.Brain_area_acronym, 'SSp-ll') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-m') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-n') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-tr') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-ul') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-un');
elseif strcmpi(areaOI, 'SSp')
  areaInds = ismember(areaTable.Brain_area_acronym, 'SSp-bfd') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-ll') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-m') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-n') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-tr') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-ul') | ...
    ismember(areaTable.Brain_area_acronym, 'SSp-un');
elseif strcmpi(areaOI, 'RSPd-RSPv')
  areaInds = ismember(areaTable.Brain_area_acronym, 'RSPd') | ...
    ismember(areaTable.Brain_area_acronym, 'RSPv');
elseif strcmpi(areaOI, 'RSP')
  areaInds = ismember(areaTable.Brain_area_acronym, 'RSPagl') | ...
    ismember(areaTable.Brain_area_acronym, 'RSPd') | ...
    ismember(areaTable.Brain_area_acronym, 'RSPv');
elseif strcmpi(areaOI, 'CA')
  areaInds = ismember(areaTable.Brain_area_acronym, 'CA1') | ...
    ismember(areaTable.Brain_area_acronym, 'CA2') | ...
    ismember(areaTable.Brain_area_acronym, 'CA3');
elseif strcmpi(areaOI, 'CA-DG')
  areaInds = ismember(areaTable.Brain_area_acronym, 'CA1') | ...
    ismember(areaTable.Brain_area_acronym, 'CA2') | ...
    ismember(areaTable.Brain_area_acronym, 'CA3') | ...
    ismember(areaTable.Brain_area_acronym, 'DG');
else
  areaInds = ismember(areaTable.Brain_area_acronym, areaOI);
end