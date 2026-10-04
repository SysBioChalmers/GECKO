%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kcatAggregation = resolveKcatAggregation(kcatAggregation, modelAdapter)
% resolveKcatAggregation  Settle which BRENDA aggregate to read.
%
% The BRENDA snapshot records both a maximum and a median per
% (EC, substrate, organism) triple. Which one a run uses can be given per
% call, or once on the model adapter as params.kcatAggregation; an
% adapter that sets neither gets 'max', which is what GECKO has always
% used.
%
% Input Arguments
% ---------------
% kcatAggregation : char
%     the per-call choice, or empty to fall back to the adapter.
% modelAdapter : ModelAdapter
%     a loaded model adapter.
%
% Returns
% -------
% kcatAggregation : char
%     'max' or 'median'.
%
% See also
% --------
% loadBRENDAdata, fuzzyKcatMatching

if isempty(kcatAggregation)
    params = modelAdapter.getParameters();
    if isfield(params,'kcatAggregation') && ~isempty(params.kcatAggregation)
        kcatAggregation = params.kcatAggregation;
    else
        kcatAggregation = 'max';
    end
end

if ~(ischar(kcatAggregation) || isstring(kcatAggregation)) || ...
        ~any(strcmpi(kcatAggregation,{'max','median'}))
    error('resolveKcatAggregation:invalidKcatAggregation', ...
        'kcatAggregation must be ''max'' or ''median''.')
end
kcatAggregation = lower(char(kcatAggregation));
end
