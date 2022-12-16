function model = makeEcModel(model,modelAdapter,geckoLight)
% makeEcModel
%   Expands a conventional genome-scale model (in RAVEN format) with enzyme
%   information and prepares the reactions for integration of enzyme usage
%   coefficients. This function contains all the steps that need to be done
%   to get a basic ec-model, without incorporating any kcat values or
%   constraints yet. This function should only have to be run once for a
%   model.
%
% Input:
%   model        a model in RAVEN format
%   modelAdapter A modelAdapter, an instance of a class inheriting the 
%                ModelAdapter class, e.g., YeastAdapter.
%   geckoLight   true if a simplified GECKO light model should be generated.
%                Optional, default is false.
%
% Ouput:
%   model        a model with a model.ec structure where enzyme and kcat
%                information are stored. Protein pseudometabolites and their
%                draw reactions are added to the model, but their usage is
%                not yet implemented (due to absent kcat values at this
%                stage).
%
% The function goes through the following steps:
%   1.  Remove gene associations from pseudoreactions.
%   2.  Invert irreversible backwards reactions.
%   3.  Correct 'rev' vector to match lb and ub vectors.
%   4.  Convert to irreversible model (splits reversible reactions).
%   5.  [Skipped with geckoLight:] Expand model to split reactions with
%       'OR' in grRules (each reaction is then catalyzed by one enzyme
%       (complex).
%   6.  [Skipped with geckoLight:] Sort identifiers (so that split reactions 
%       remain close to each other, not real function, just makes it tidier.
%   7.  Make empty model.ec structure, that will contain enzyme and kcat
%       information. One entry per reaction, where isoenzymes have multiple
%       entries. This model.ec structure will later be populated with kcat values.
%       For geckoLight the structure is different, where each reaction can have
%       multiple isozymes.
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

if nargin<3
    geckoLight=false;
elseif ~islogical(geckoLight) && ~(geckoLight == 0) && ~(geckoLight == 1)
    error('geckoLight should be either true or false')
end

if geckoLight
    ec.geckoLight=true;
else
    ec.geckoLight=false;
end

[geckoPath, prevDir] = findGECKOroot();
uniprotDB = loadDatabases(modelAdapter);


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
if ~geckoLight
    model=sortIdentifiers(model);
end

%7: Make ec-extension structure, one for gene-associated reaction.
%   The structure is different for light and full models
rxnWithGene  = find(sum(model.rxnGeneMat,2));
if ~geckoLight
    ec.rxns      = model.rxns(rxnWithGene);
    emptyCell    = cell(numel(rxnWithGene),1);
    emptyVect    = zeros(numel(rxnWithGene),1);
    if isfield(model,'eccodes')
        ec.eccodes = model.eccodes(rxnWithGene);
    else
        ec.eccodes = emptyCell;
        %TODO: loading of external ec-code database, possibly from uniprot, which
        %should directly be parsed to the model, so that the eccodeDB entries match
        %model.rxns. This is problematic, because Uniprot can contain multiple
        %eccodes, while eccodes might not exactly match the reactions that are
        %catalyzed. Only needed for classic GECKO matching, might just copy how it
        %was dealt with there
    end
    ec.kcat      = emptyVect;
    ec.source    = emptyCell; % Strings, like 'dlkcat', 'manual', 'brenda', etc.
    ec.notes     = emptyCell; % Additional comments
else
    %Different strategy for GECKO light: Each reaction can exist multiple times in 
    %ec.rxns and similar fields - one time per isozyme. The number of copies is
    %the number of ORs in the GPR + 1
    numOrs = count(model.grRules(rxnWithGene), ' or ');
    cpys = numOrs + 1;
    prevNumRxns = length(numOrs);
    cpyIndices = repelem(rxnWithGene, cpys);
    ec.rxns      = model.rxns(cpyIndices);
    emptyCell    = cell(numel(ec.rxns),1);
    emptyVect    = zeros(numel(ec.rxns),1);
    if isfield(model,'eccodes')
        ec.eccodes = model.eccodes(cpyIndices);
    else
        ec.eccodes = emptyCell;
        %TODO: loading of external ec-code database, possibly from uniprot, which
        %should directly be parsed to the model, so that the eccodeDB entries match
        %model.rxns. This is problematic, because Uniprot can contain multiple
        %eccodes, while eccodes might not exactly match the reactions that are
        %catalyzed. Only needed for classic GECKO matching, might just copy how it
        %was dealt with there
    end
    
    ec.kcat      = emptyVect;
    ec.source    = emptyCell; % Strings, like 'dlkcat', 'manual', 'brenda', etc.
    ec.notes     = emptyCell; % Additional comments
end
    
%8: Gather enzyme information via UniprotDB
uniprotCompatibleGenes = modelAdapter.getUniprotCompatibleGenes(model);
[Lia,Locb]      = ismember(uniprotCompatibleGenes,uniprotDB.genes);
ec.genes        = model.genes(Lia); %Will often be duplicate of model.genes, but is done here to prevent issues when it is not.
ec.enzymes      = uniprotDB.ID(Locb(Lia));
ec.mw           = uniprotDB.MW(Locb(Lia));
ec.sequence     = uniprotDB.seq(Locb(Lia));
%Additional info
ec.concs        = nan(numel(ec.genes),1); % To be filled with proteomics data when available
%TODO: load Uniprot IDs from model annotation instead of from uniprotDB?
%To offer a choice, should then still be matched to a uniprotDB to obtain
%mw and sequence.
% if isfield(model,'geneMiriams')
%     uniprotDB = extractMiriam(model.geneMiriams,'uniprot');
% end

%9: Only parse rxns associated to genes
if ~geckoLight
    ec.rxnEnzMat = zeros(numel(rxnWithGene),numel(ec.genes)); % Non-zeros will indicate the number of subunits
    for r=1:numel(rxnWithGene)
        ec.rxns(r) = model.rxns(rxnWithGene(r)); %why is this needed, already set above? Test to remove later
        rxnGenes   = model.genes(find(model.rxnGeneMat(rxnWithGene(r),:)));
        [~,locEnz] = ismember(rxnGenes,ec.genes); % Could also parse directly from rxnGeneMat, but some genes might be missing from Uniprot DB
        %Changed because this allows for a selection of zero genes - all genes are not always found in 
        %locEnz = ismember(ec.genes, rxnGenes); % Could also parse directly from rxnGeneMat, but some genes might be missing from Uniprot DB
        if locEnz ~= 0
            ec.rxnEnzMat(r,locEnz) = 1; %Assume 1 copy per subunit or enzyme, can be modified later
        end
    end
else
    %For light models, we need to split up all GPRs
    ec.rxnEnzMat = zeros(numel(ec.rxns),numel(ec.genes)); % Non-zeros will indicate the number of subunits
    nextIndex = 1;
    rxnToCheck={};
    for i=1:prevNumRxns
        %ind is the index in the model, not to confuse with the index in the ec struct (i),
        %which only contains reactions with GPRs.
        ind = rxnWithGene(i); 
        %Check that it doesn't contain nested 'and' and 'or' relations and
        %print a warning if it does
        %TODO: This warning is not right - should do something similar to the old GECKO, i.e., 
        %warn if you have another form of the GPR than (subunit1 and subunit2 and ...) or (subunit3 and subunit4 and ...) or ...
        %Same problem in expandModel
        if ~isempty(strfind(model.grRules{ind},' and '))
            rxnToCheck{end+1,1}=model.rxns{ind};
        end
        %Get rid of all '(' and ')' since I'm not looking at complex stuff
        %anyways
        geneString=model.grRules{ind};
        geneString=strrep(geneString,'(','');
        geneString=strrep(geneString,')','');
        geneString=strrep(geneString,' or ',';');
        
        if (numOrs(i) == 0)
            geneNames = {geneString};
        else
            %Split the string into gene names
            geneNames=regexp(geneString,';','split');
        end
        
        %Now loop through the isozymes and set the rxnGeneMat
        for j = 1:length(geneNames)
            %Find the gene in the gene list If ' and ' relationship, first
            %split the genes
            fnd = strfind(geneNames{j},' and ');
            if ~isempty(fnd)
                andGenes=regexp(geneNames{j},' and ','split');
                ec.rxnEnzMat(nextIndex,ismember(ec.genes,andGenes)) = 1; %should be subunit stoichoimetry
            else
                ec.rxnEnzMat(nextIndex,ismember(ec.genes,geneNames(j)))=1;%should be subunit stoichoimetry
            end
            nextIndex = nextIndex + 1;
        end
    end
    
    if ~isempty(rxnToCheck) 
        warning(['The following reactions contain nested and/or-relations, which might not have been interpreted correctly:%s\n' ...
            strjoin(rxnToCheck,', ')],'');
    end
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
    drawRxns.rxns            = strcat('draw_',proteinMets.mets);
    drawRxns.mets            = cell(numel(drawRxns.rxns),1);
    drawRxns.stoichCoeffs    = cell(numel(drawRxns.rxns),1);
    for i=1:numel(drawRxns.mets)
        drawRxns.mets{i}         = {'prot_pool',proteinMets.mets{i}};
        drawRxns.stoichCoeffs{i} = [1,-(ec.mw(uniprotSortId(i)))/1000];
    end
    drawRxns.lb              = zeros(numel(drawRxns.rxns),1);
    drawRxns.grRules         = ec.genes(uniprotSortId);
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
