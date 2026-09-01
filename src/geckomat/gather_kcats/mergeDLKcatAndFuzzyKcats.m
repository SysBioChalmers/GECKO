function mergedKcatList = mergeDLKcatAndFuzzyKcats(kcatListDLKcat, kcatListFuzzy, varargin)
% mergeDLKcatAndFuzzyKcats  Merge DLKcat and BRENDA fuzzy matching results.
%
% Deprecated thin wrapper around mergeKcats, kept for the common DLKcat +
% fuzzy-BRENDA case. New code should call mergeKcats directly, which also
% handles a single list that mixes several sources (e.g. an
% OpenKineticsPredictor result with BRENDA/Sabio-RK/CataPro rows).
%
% Merges the results from DLKcat and fuzzy matching to the BRENDA database.
% Order of preference:
%
% 1. BRENDA match with correct E.C. number, with origin (see below) not
%    lower than the specified topOriginLimit.
% 2. DLKcat match.
% 3. BRENDA match with correct E.C. number, with origin below topOriginLimit
%    but not lower than the bottomOriginLimit.
% 4. BRENDA match with wildcards in the E.C. number, with not more wildcards
%    than wildcardLimit, and origin not lower than the bottomOriginLimit.
%
% Parameters
% ----------
% kcatListDLKcat : struct
%     kcatList derived from readDLKcatOutput.
% kcatListFuzzy : struct
%     kcatList derived from fuzzyKcatMatching.
%
% Name-Value Arguments
% --------------------
% topOriginLimit : double
%     origin limit for prioritized BRENDA matches. Origin is explained in
%     more detail below (default 6).
% bottomOriginLimit : double
%     origin limit for low priority BRENDA matches. Origin is explained in
%     more detail below (default 6).
% wildcardLimit : double
%     maximum number of wildcards in E.C. number of BRENDA matches
%     (default 3).
%
% Returns
% -------
% mergedKcatList : struct
%     merged kcatList with the same field layout produced by mergeKcats
%     (source, kcatSource, rxns, genes, substrates, kcats, eccodes,
%     wildcardLvl, origin), keeping only the rows selected by the order of
%     preference above.
%
% Notes
% -----
% The origin parameter:
%
% - 1 : correct organism, correct substrate, kcat.
% - 2 : any organism, correct substrate, kcat.
% - 3 : correct organism, any substrate, kcat.
% - 4 : any organism, any substrate, kcat.
% - 5 : correct organism, specific activity.
% - 6 : any organism, specific activity.
%
% Example of wildcards in E.C. number:
%
% - 0 : 1.1.1.3   glycerol-3-phosphate dehydrogenase (NAD+).
% - 1 : 1.1.1.-   oxidoreductase, acting on the CH-OH group of donors, with
%   NAD+ or NADP+ as acceptor.
% - 2 : 1.1.-.-   oxidoreductase, acting on the CH-OH group of donors.
% - 3 : 1.-.-.-   oxidoreductase.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     mergedKcatList = mergeDLKcatAndFuzzyKcats(kcatListDLKcat, kcatListFuzzy);
%     mergedKcatList = mergeDLKcatAndFuzzyKcats(kcatListDLKcat, kcatListFuzzy, 'wildcardLimit', 2);
%
% See also
% --------
% mergeKcats, readDLKcatOutput, fuzzyKcatMatching, assignKcatValues

p = parseGECKOargs(varargin, { ...
    'topOriginLimit',    6; ...
    'bottomOriginLimit', 6; ...
    'wildcardLimit',     3});
topOriginLimit    = p.topOriginLimit;
bottomOriginLimit = p.bottomOriginLimit;
wildcardLimit     = p.wildcardLimit;

warning(['mergeDLKcatAndFuzzyKcats is deprecated; use mergeKcats instead. ' ...
    'The old name will be removed in a future release.']);

% Fuzzy list first, then DLKcat, so surviving rows keep the legacy order
% (fuzzy block followed by DLKcat block).
mergedKcatList = mergeKcats({kcatListFuzzy, kcatListDLKcat}, ...
    {'database_top', 'dlkcat', 'database_bottom'}, ...
    'topOriginLimit', topOriginLimit, ...
    'bottomOriginLimit', bottomOriginLimit, ...
    'wildcardLimit', wildcardLimit);
end
