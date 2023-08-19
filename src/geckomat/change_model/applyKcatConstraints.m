function model = applyKcatConstraints(model,updateRxns)
% applyKcatConstraints
%   Applies kcat-derived enzyme constraints to an ecModel. Existing enzyme
%   constraints are first removed (unless updateRxns is provided), and new
%   constraints are defined based on the content of model.ec.kcat.
%
% Input:
%   model       an ecModel in GECKO 3 format (with ecModel.ec structure)
%   updateRxns  if not all enzyme constraints should be updated, this can
%               be given as either a logical vector of length
%               model.ec.rxns, a vector of model.ec.rxns indices, or a
%               (cell array of) string(s) with model.ec.rxns identifiers.
%               For light models, these reactions should match model.rxns.
%
% Output:
%   model       ecModel where reactions are constrained by enzyme usage
%               if a kcat value was provided for the reaction-enzyme pair
%               in model.ec.kcat
%
% Usage: model = applyKcatConstraints(model,updateRxns);

%these lines are for the nargin lines below only
if (model.ec.geckoLight)
    rxns = model.rxns;
else
    rxns = model.ec.rxns;
end

if nargin<2
    updateRxns = true(numel(rxns),1);
elseif isnumeric(updateRxns)
    updateRxnsLog = false(numel(rxns),1);
    updateRxnsLog(updateRxns) = true;
    updateRxns = updateRxnsLog;
elseif iscellstr(updateRxns) || ischar(updateRxns) || isstring(updateRxns)
    updateRxnsIds = convertCharArray(updateRxns);
    updateRxns = ismember(rxns,updateRxnsIds);
end
 
if isempty(find(updateRxns, 1)) || isempty(updateRxns)
     error('No reaction to update or updateRxns is logical but without any true value')
end

if ~isfield(model,'ec')
    error(['No model.ec structure could be found: the provided model is'...
           ' not a valid GECKO3 ecModel. First run makeEcModel(model).'])
end
if all(model.ec.kcat==0)
    printOrange('WARNING: No kcat values are provided in model.ec.kcat, model remains unchanged.\n')
    return
end

%Clear existing incorporation of enzyme usage
if ~model.ec.geckoLight
    protMetIdx = startsWith(model.mets,'prot_') & ~strcmp(model.mets,'prot_pool');
    metabolRxn = unique(model.ec.rxns(updateRxns));
    metabolRxn = ismember(model.rxns,metabolRxn);
    model.S(protMetIdx,metabolRxn) = 0;
end %no clearing needed for GECKO Light, will be overwritten below
%For normal GECKO formulation (full model), where each enzyme is explicitly considered
if ~model.ec.geckoLight 
    %Column 1 = rxn idx; 2 = enzyme idx; 3 = subunit copies; 4 = kcat; 5 = MW
    newKcats=zeros(numel(updateRxns)*10,5);
    updateRxns=find(updateRxns);
    kcatFirst=0;
    for i=1:numel(updateRxns)
        j=updateRxns(i);
        enzymes   = find(model.ec.rxnEnzMat(j,:));
        kcatLast  = kcatFirst+numel(enzymes);
        kcatFirst = kcatFirst+1;
        newKcats(kcatFirst:kcatLast,1) = j;
        newKcats(kcatFirst:kcatLast,2) = enzymes;
        newKcats(kcatFirst:kcatLast,3) = model.ec.rxnEnzMat(j,enzymes);
        newKcats(kcatFirst:kcatLast,4) = model.ec.kcat(j);
        newKcats(kcatFirst:kcatLast,5) = model.ec.mw(enzymes);
        kcatFirst = kcatLast;
    end
    newKcats(kcatLast+1:end,:)=[];
    
    sel = newKcats(:,4) > 0; %Only apply to non-zero kcat
    newKcats(sel,4) = newKcats(sel,4) * 3600; %per second -> per hour
    newKcats(sel,4) = newKcats(sel,5) ./ newKcats(sel,4); %MW/kcat
    newKcats(sel,4) = newKcats(sel,3) .* newKcats(sel,4); %Multicopy subunits.
    newKcats(~sel,4) = 0; %Results in zero-cost
    
    %Replace rxns and enzymes with their location in model
    [~,newKcats(:,1)] = ismember(model.ec.rxns(newKcats(:,1)),model.rxns);
    [~,newKcats(:,2)] = ismember(strcat('prot_',model.ec.enzymes(newKcats(:,2))),model.mets);
    linearIndices     = sub2ind(size(model.S),newKcats(:,2),newKcats(:,1));
    model.S(linearIndices) = -newKcats(:,4); %Substrate = negative
    
else %GECKO light formulation, where prot_pool represents all usages
    prot_pool_idx = find(ismember(model.mets,'prot_pool'));
    %first remove the prefix of all rxns
    modRxns     = extractAfter(model.ec.rxns,4);
    % Map ec-reactions to model.rxns
    [hasEc,~]   = ismember(model.rxns,modRxns);
    hasEc       = find(hasEc & updateRxns);
    [~,rxnIdx]   = ismember(modRxns,model.rxns);
    for i = 1:numel(hasEc)
        % Get all isozymes per reaction
        ecIdx = find(rxnIdx == hasEc(i));
        % ecIdx = strcmp(model.rxns(hasEc(i)),modRxns);
        % Multiply enzymes with their MW (they are then automatically
        % summed per reaction), and divide by their kcat, to get a vector
        % of MW/kcat values.
        MWkcat = (model.ec.rxnEnzMat(ecIdx,:) * model.ec.mw) ./ model.ec.kcat(ecIdx);
        % If kcat was zero, MWkcat is Inf. If no enzyme info was found,
        % MWkcat is NaN. Correct both back to zero
        MWkcat(isinf(MWkcat) | isnan(MWkcat)) = 0;
        % Select the lowest MW/kcat (= most efficient), and convert to /hour
        model.S(prot_pool_idx, hasEc(i)) = -min(MWkcat/3600);
    end
end
end

