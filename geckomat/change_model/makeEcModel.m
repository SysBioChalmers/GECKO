function model = makeEcModel(model,geckoLight)
% makeEcModel
%   Expands a conventional genome-scale model (in RAVEN format) with enzyme
%   information and prepares the reactions for integration of enzyme usage
%   coefficients. This function contains all the steps that need to be done
%   to get a basic ec-model, without incorporating any kcat values or
%   constraints yet. This function should only have to be run once for a
%   model.
%
% Input:
%   model       a model in RAVEN format
%   geckoLight  true if a simplified GECKO light model should be generated.
%               Optional, default is false.
%
% Ouput:
%   model       a model with a model.ec structure where enzyme and kcat
%               information are stored. Protein pseudometabolites and their
%               draw reactions are added to the model, but their usage is
%               not yet implemented (due to absent kcat values at this
%               stage).
%
% The function goes through the following steps:
%   1.  Remove gene associations from pseudoreactions.
%   2.  Invert irreversible backwards reactions.
%   3.  Correct 'rev' vector to match lb and ub vectors.
%   4.  Convert to irreversible model (splits reversible reactions).
%   5.  [Skipped with geckoLight:] Expand model to split reactions with
%       'OR' in grRules (each reaction is then catalyzed by one enzyme
%       (complex).
%   6.  Sort identifiers (so that split reactions remain close to each
%       other, not real function, just makes it tidier.
%   7.  Make empty model.ec structure, that will contain enzyme and kcat
%       information. One entry per reaction, where isoenzymes have multiple
%       entries [not in geckoLight]. This model.ec structure will later be
%       populated with kcat values.
%   8.  Add enzyme information fields to model.ec structure: MW, sequence.
%   9.  Populate model.ec structure (from step 8) with information from
%       each reaction.
%   10. [Skipped with geckoLight:] Add proteins as pseudometabolites.
%   11. Add prot_pool pseudometabolite.
%   12. [Skipped with geckoLight:] Add draw reactions for the protein
%       pseudometabolites.
%   13. Add protein pool reaction, without upper bound.
%
%   Note that while protein pseudometabolites, draw & pool reactions might
%   be added to the model, the enzyme usage is not yet incorporated in each
%   metabolic reaction, so enzymes will not be used. applyKcatConstraints
%   incorporates protein pseudometabolites in reactions as enzyme usages by
%   applying the specified kcats as constraints.

if nargin<2
    geckoLight=false;
elseif ~islogical(geckoLight)
    error('geckoLight should be either true or false')
end

if geckoLight
    ec.geckoLight=true;
else
    ec.geckoLight=false;
end

[geckoPath, prevDir] = findGECKOroot();
uniprotDB = loadDatabases;

%1: Remove gene rules from pseudoreactions (if any):
for i = 1:length(model.rxns)
    if endsWith(model.rxnNames{i},' pseudoreaction')
        model.grRules{i}      = '';
        model.rxnGeneMat(i,:) = zeros(1,length(model.genes));
    end
end

%2: Swap direction of reactions that are defined to only carry negative flux
to_swap=model.lb < 0 & model.ub == 0;
model.S(:,to_swap)=-model.S(:,to_swap);
model.ub(to_swap)=-model.lb(to_swap);
model.lb(to_swap)=0;

%Delete blocked rxns (LB = UB = 0). Best not to do, as you cannot unblock
%these reactions later once removed. Better to do this once you run
%analysis.
% to_remove = logical((model.lb == 0).*(model.ub == 0));
% model     = removeReactions(model,model.rxns(to_remove),true,true,true);

%3: Correct rev vector: true if LB < 0 & UB > 0, or it is an exchange reaction:
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || sum(model.S(:,i) ~= 0) == 1
        model.rev(i) = true;
    end
end

%4: Make irreversible model (appends _REV to reaction IDs to indicate reverse
%reactions)
model=convertToIrrev(model);

%5: Expand model, to separate isoenzymes (appends _EXP_* to reaction IDs to
%indicate duplication)
if ~geckoLight
    model=expandModel(model);
end

%6: Sort reactions, so that reversible and isoenzymic reactions are kept near
model=sortIdentifiers(model);

%7: Make ec-extension structure, one for each kcat value
emptyCell    = cell(numel(model.rxns)*5,1); %Prepare plenty of space
emptyVect    = zeros(numel(model.rxns)*5,1);
ec.rxns      = emptyCell;
ec.eccodes   = emptyCell; % Can also be parsed from model. Only needed for classic GECKO matching
ec.rxnEnzMat = zeros(numel(model.rxns)*5,numel(model.genes)); % Non-zeros will indicate the number of subunits
%ec.reverse   = logical(emptyVect); % Whether reaction is in reverse. Might be unnecessary, can also be gathered from the reaction ID ('_REV')
ec.kcat      = emptyVect;
ec.source    = emptyCell; % Strings, like 'dlkcat', 'manual', 'kcatdb'
% (referring to BRENDA and/or SABIO-RK). Should be standardized for
% applyKcatConstraints to function. Can alternatively be numeric (with each
% number having a dedicated meaning), to avoid spelling mistakes and users
% defining other types of sources?
ec.notes     = emptyCell;

%8: Gather enzyme information via UniprotDB
[Lia,Locb]      = ismember(model.genes,uniprotDB.genes);
ec.genes        = model.genes(Lia); %Will often be duplicate of model.genes, but is done here to prevent issues when it is not.
ec.enzymes      = uniprotDB.ID(Locb);
ec.mw           = uniprotDB.MW(Locb);
ec.sequence     = uniprotDB.seq(Locb);
%Additional info
ec.concs        = nan(numel(ec.genes),1); % To be filled with proteomics data when available
%TODO: load Uniprot IDs from model annotation instead of from uniprotDB?
%To offer a choice, should then still be matched to a uniprotDB to obtain
%mw and sequence.
% if isfield(model,'geneMiriams')
%     uniprotDB = extractMiriam(model.geneMiriams,'uniprot');
% end


%Extract model provided eccodes
%TODO: loading of external ec-code database, possibly from uniprot, which
%should directly be parsed to the model, so that the eccodeDB entries match
%model.rxns. This is problematic, because Uniprot can contain multiple
%eccodes, while eccodes might not exactly match the reactions that are
%catalyzed. Only needed for classic GECKO matching, might just copy how it
%was dealt with there
%ec.eccodes      = uniprotDB.eccodes(Locb);

%9: Only parse rxns associated to genes
rxnWithGene = find(sum(model.rxnGeneMat,2));
ecCount=1;
for r=1:numel(rxnWithGene)
    ec.rxns(r) = model.rxns(rxnWithGene(r));
    rxnGenes   = model.genes(find(model.rxnGeneMat(rxnWithGene(r),:)));
    [~,locEnz] = ismember(rxnGenes,ec.genes); % Could also parse directly from rxnGeneMat, but some genes might be missing from Uniprot DB
    ec.rxnEnzMat(r,locEnz) = 1; %Assume 1 copy per subunit or enzyme, can be modified later
end
%Remove empty fields
for f=fieldnames(ec)'
    ec.(f{:})(ecCount-1:end,:)=[];
end

%10: Add proteins as pseudometabolites
if ~geckoLight
    [proteinMets.mets, uniprotSortId] = sort(ec.enzymes);
    proteinMets.mets         = strcat('prot_',proteinMets.mets);
    proteinMets.metNames     = proteinMets.mets;
    proteinMets.compartments = 'c';
    proteinMets.metMiriams   = repmat({struct('name',{{'sbo'}},'value',{{'SBO:0000252'}})},numel(proteinMets.mets),1);
    proteinMets.metNotes     = repmat({'Enzyme-usage pseudometabolite'},numel(proteinMets.mets),1);
    model = addMets(model,proteinMets);
end
%11: Add protein pool pseudometabolite
pool.mets         = 'prot_pool';
pool.metNames     = pool.mets;
pool.compartments = 'c';
pool.metNotes     = 'Enzyme-usage protein pool';
model = addMets(model,pool);

%12: Add protein draw reactions.
if ~geckoLight
    %TODO: Separate function to switch to exchange reactions for proteomics
    %data
    drawRxns.rxns            = strcat('draw_',proteinMets.mets);
    drawRxns.mets            = cell(numel(drawRxns.rxns),1);
    drawRxns.stoichCoeffs    = cell(numel(drawRxns.rxns),1);
    for i=1:numel(drawRxns.mets)
        drawRxns.mets{i}         = {'prot_pool',proteinMets.mets{i}};
        drawRxns.stoichCoeffs{i} = [1,-(ec.mw(uniprotSortId(i)))/1000];
    end
    drawRxns.lb              = zeros(numel(drawRxns.rxns),1);
    drawRxns.grRules         = ec.gene(uniprotSortId);
    model = addRxns(model,drawRxns);
end

%13: Add protein pool reaction (with open UB)
poolRxn.rxns            = 'prot_pool_exchange';
poolRxn.mets            = {'prot_pool'};
poolRxn.stoichCoeffs    = {-1};
poolRxn.lb              = 0;
model = addRxns(model,poolRxn);

model.ec=ec;
end
