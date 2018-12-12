function ecModel = getSubset_ecModel(model,generalModel)
% getSubset_ecModel
%  
% 	Get a general ecModel and a context-specific metabolic model 
%   (strain/cell-line/tissue) and maps all of its reactions in the general 
%   ecModel in order to get an enzyme constrained version of the context-
%   specific model.
%
%   generalModel    Enzyme-constrained version of  general model of
%                   metabolism for a given species
%   model           Reduced model (subset of the general model) for a
%                   specific strain (microbes) or cell-line/tissues (mammals)
%
%   ecModel         Enzyme-constrained version of the context-specific model
%
% usage: ecModel = getSubset_ecModel(model,generalModel)
% 
% Ivan Domenzain.      Last edited: 2018-12-11
%

current = pwd;
%If model is in COBRA format convert to RAVEN-type 
if isfield(model,'rules')
    model = ravenCobraWrapper(model);
end
irrevModel = convertToIrrev(model);
ecModel    = irrevModel;
%Map each of the irrevModel rxns into the general model
enzymes = {};
for i=1:length(irrevModel.rxns)
    rxn      = irrevModel.rxns{i};
    rxnIndxs = find(contains(generalModel.rxns,rxn)); 
    %Exclude the original rxn index ('avoids adding it again)
    rxnIndxs = rxnIndxs(~strcmpi(generalModel.rxns(rxnIndxs),rxn));
    %As the model as been converted to irreversible format, correspondent
    %additional reactions for forward and backward reactions are added
    %separately
    if ~contains(rxn,'_REV')
        rxnIndxs = rxnIndxs(~contains(generalModel.rxns(rxnIndxs),'_REV'));
    else
        rxnIndxs = rxnIndxs(contains(generalModel.rxns(rxnIndxs),'_REV')); 
    end
    %Just keep the new reactions (isoenzymes, arm reactions)
    rxnIndxs    = rxnIndxs(startsWith(generalModel.rxns(rxnIndxs),[rxn 'No'])|...
                           strcmpi(generalModel.rxns(rxnIndxs),['arm_' rxn]));
    ecReactions = generalModel.rxns(rxnIndxs);
    
    if ~isempty(rxnIndxs)
        nR = length(rxnIndxs);
        grRules = cell(1,nR);
        %Create rxnsToAdd structure
        rxnsToAdd.rxns     = ecReactions;
        rxnsToAdd.rxnNames = generalModel.rxnNames(rxnIndxs);
        rxnsToAdd.lb       = generalModel.lb(rxnIndxs);
        rxnsToAdd.ub       = generalModel.ub(rxnIndxs);
        rxnsToAdd.c        = generalModel.lb(rxnIndxs);
        grRules(:)         = irrevModel.grRules(i);
        rxnsToAdd.grRules  = grRules;
        %Preallocate fields concerning the metabolites involved in each rxn
        mets           = cell(1,nR);
        StoichCoeffs   = cell(1,nR);
        rxnsToAdd.mets = cell(1,nR);
        %Loop through each of the mapped Indexes to get metabolites information
        for j=1:length(rxnIndxs)
            index = rxnIndxs(j);
            metsIndxs         = find(generalModel.S(:,index));
            StoichCoeffs{j}   = full(generalModel.S(metsIndxs,index));
            mets{j}           = generalModel.mets(metsIndxs);
            metNames          = generalModel.metNames(metsIndxs);
            compartments      = generalModel.metComps(metsIndxs);
            compartments      = generalModel.comps(compartments);
            bVector           = generalModel.b(metsIndxs);
            enzNames          = mets{j}(contains(mets{j},'prot_'));
            rxnsToAdd.mets{j} = mets{j};
            %Avoid ambiguities checking for proteins also in the general  
            %ecModel enzymes field
            enzNames  = strrep(enzNames,'prot_','');
            enzNames  = intersect(enzNames,generalModel.enzymes);
            enzymes   = [enzymes;enzNames(:)];
            %Add new metabolites for each reaction
            IA      = ismember(mets{j},ecModel.mets);
            mets{j} = mets{j}(~IA);
            if ~isempty(mets{j})
                metsToAdd.mets         = mets{j};
                metsToAdd.metNames     = metNames(~IA);
                metsToAdd.compartments = compartments(~IA);
                metsToAdd.bVector      = bVector(~IA);
                %Use RAVEN function for adding new metabolites
                ecModel = addMets(ecModel,metsToAdd);
            end
        end
        rxnsToAdd.stoichCoeffs = StoichCoeffs;
        %Use RAVEN function for adding new reactions
        ecModel = addRxns(ecModel,rxnsToAdd,1,[],false,true); 
        disp(['Added additional rxns for: ' rxn])
        if any(contains(rxnsToAdd.rxns,[rxn 'No']))
            ecModel = removeReactions(ecModel,rxn);
        end
    end
end
cd ../change_model
enzymes = unique(enzymes);
%Create additional fields in ecModel:
ecModel.enzymes  = cell(0,1);
ecModel.enzGenes  = cell(0,1);
ecModel.enzNames  = cell(0,1);
ecModel.MWs       = zeros(0,1);
ecModel.sequences = cell(0,1);
ecModel.pathways  = cell(0,1);
%Load protein databases for enzyme-genes matching
load ('../../Databases/ProtDatabase.mat')
%Add protein usage reactions
for i = 1:length(enzymes)
    disp(enzymes{i})
    ecModel = addProtein(ecModel,enzymes{i},kegg,swissprot);
end
cd (current)
end
