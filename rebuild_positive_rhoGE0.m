function infraslowAnalyses = rebuild_positive_rhoGE0(infraslowAnalyses)
% Rebuild infraslowAnalyses.spikingPupilCohPositive with the rho >= 0 rule (as in the UoL/Allen code)
% instead of rho > 0, without recomputing coherence (2026-11-08).
% Unit-to-pupil coherence is computed unit by unit (coherence_pupil.m), so the positive-population entry of
% a recording is the 'all' entry restricted to units with rho >= 0. The function first checks, for every
% recording, that the stored positive entry equals the 'all' entry restricted to rho > 0; it errors otherwise.
% Exception: fullInterpCoherence.phaseConf is one 2 x frequency array per recording (not per unit) and is not
% used by any analysis; it is not checked and is set to NaN in the rebuilt entries.

all = infraslowAnalyses.spikingPupilCoh;
pos = infraslowAnalyses.spikingPupilCohPositive;
areas = fieldnames(all);
nChecked = 0; nAdded = 0;
for iArea = 1:numel(areas)
  a = areas{iArea};
  if ~isfield(pos, a), pos.(a) = {}; end
  for iRec = 1:numel(all.(a))
    A = all.(a){iRec};
    if ~isstruct(A) || ~isfield(A, 'rSpearman'), continue; end
    r = A.rSpearman(:);
    if ~isfield(A, 'fullInterpCoherence')
      continue  % no coherence for this recording (e.g. no pupil); nothing to rebuild
    end
    nAll = numel(r);
    keepOld = r > 0; keepNew = r >= 0;
    hadOld = iRec <= numel(pos.(a)) && isstruct(pos.(a){iRec}) && isfield(pos.(a){iRec}, 'fullInterpCoherence');
    if any(keepOld) ~= hadOld
      error('rebuild:mismatch', '%s rec %d: positive entry presence does not match rho > 0 units', a, iRec);
    end
    if hadOld
      P = pos.(a){iRec};
      checkSubset(A.fullInterpCoherence, P.fullInterpCoherence, keepOld, nAll, sprintf('%s rec %d', a, iRec));
      if ~isequal(A.experimentIndex, P.experimentIndex)
        error('rebuild:mismatch', '%s rec %d: experiment index differs', a, iRec);
      end
      nChecked = nChecked + 1;
    end
    if ~any(keepNew)
      continue
    end
    N = struct();
    N.fullInterpCoherence = subsetRows(A.fullInterpCoherence, keepNew, nAll);
    N.spikeCounts = A.spikeCounts(keepNew);
    N.firingRates = A.firingRates(keepNew);
    um = A.unitMask; idx = find(um); um(idx(~keepNew)) = false; N.unitMask = um;
    N.experimentIndex = A.experimentIndex;
    N.timeOfCompletion = datetime;
    N.rebuiltFrom = 'spikingPupilCoh, rho >= 0 (rebuild_positive_rhoGE0.m, 2026-11-08)';
    if any(keepNew & ~keepOld), nAdded = nAdded + sum(keepNew & ~keepOld); end
    pos.(a){iRec} = N;
  end
end
infraslowAnalyses.spikingPupilCohPositive = pos;
fprintf('rebuild_positive_rhoGE0: %d recording entries checked identical to the rho > 0 subset; %d rho = 0 unit entries added\n', ...
  nChecked, nAdded);
end


function S = subsetRows(C, keep, nAll)
S = C;
f = fieldnames(C);
for k = 1:numel(f)
  v = C.(f{k});
  if strcmp(f{k}, 'phaseConf'), S.phaseConf = nan(size(v)); continue; end
  if (isnumeric(v) || islogical(v)) && size(v, 1) == nAll && nAll > 1
    S.(f{k}) = v(keep, :, :);
  elseif (isnumeric(v) || islogical(v)) && nAll == 1 && size(v, 1) == 1
    S.(f{k}) = v;  % one unit: keep is true here (otherwise not called)
  end
end
end


function checkSubset(C, P, keep, nAll, label)
f = fieldnames(C);
for k = 1:numel(f)
  v = C.(f{k});
  if strcmp(f{k}, 'phaseConf'), continue; end
  if ~isfield(P, f{k}), error('rebuild:mismatch', '%s: field %s missing in positive entry', label, f{k}); end
  w = P.(f{k});
  if (isnumeric(v) || islogical(v)) && size(v, 1) == nAll && nAll > 1
    vs = v(keep, :, :);
  else
    vs = v;
  end
  if ~isequaln(vs, w)
    % allow tiny floating-point differences only
    if isnumeric(vs) && isequal(size(vs), size(w)) && max(abs(double(vs(:)) - double(w(:))), [], 'omitnan') < 1e-9
      continue
    end
    error('rebuild:mismatch', '%s: field %s differs from the rho > 0 subset of the all-units entry', label, f{k});
  end
end
end
