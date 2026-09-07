function tieMap = isozymeTieMap(rxnIds,kcat0,varargin)
% isozymeTieMap  Index of the parameter each isozyme copy follows.
%
% A reaction catalysed by several isozymes gets one ec.rxns row per
% isozyme, each with its own kcat. Where the kcat assignment could tell
% them apart, those are genuinely different parameters; where it could
% not, every copy shares the same value from the same source, and an
% untied search is then free to invent a difference between them that
% means nothing biologically. Copies group together (and share one
% parameter) when they share a base reaction (the part of the rxn id
% before its "_EXP_" isozyme suffix), a prior value, and -- when sources
% is given -- a source.
%
% Parameters
% ----------
% rxnIds : cell
%     ec.rxns entries.
% kcat0 : double
%     prior kcat per row, same length as rxnIds.
%
% Name-Value Arguments
% --------------------
% sources : cell
%     ec.source entries, same length as rxnIds (default: not used, so
%     grouping ignores source).
% relTol : double
%     relative tolerance for treating two prior values as equal (default
%     1e-9).
%
% Returns
% -------
% tieMap : double
%     one entry per position. tieMap(i) == i for a free parameter,
%     otherwise the (lower) index of the representative it follows.
%
% See also
% --------
% cmaesKcatTuning, screenKcatLeverage
p = parseGECKOargs(varargin, { ...
    'sources', []; ...
    'relTol', 1e-9});
sources = p.sources;
relTol  = p.relTol;

n = numel(kcat0);
if numel(rxnIds) ~= n
    error('rxnIds has %d entries; kcat0 has %d.', numel(rxnIds), n)
end
tieMap = (1:n)';

keys = cell(n,1);
for i = 1:n
    base = baseReaction(rxnIds{i});
    if kcat0(i) > 0
        logKey = sprintf('%d', round(log(kcat0(i)) / relTol));
    else
        logKey = 'NONE';
    end
    if isempty(sources)
        srcKey = '';
    else
        srcKey = char(sources{i});
    end
    keys{i} = strjoin({base, logKey, srcKey}, char(31));
end

[~, ~, groupIdx] = unique(keys);
for g = 1:max(groupIdx)
    members = find(groupIdx == g);
    if numel(members) > 1
        tieMap(members) = members(1);
    end
end
end

function base = baseReaction(rxnId)
% The reaction an ec.rxns entry belongs to, isozyme suffix removed.
parts = strsplit(rxnId, '_EXP_');
base = parts{1};
end
