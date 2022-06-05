function model = makeEcModel(model)
%preprocessModel
%
% Makes the model ec-ready.
%

%Remove gene rules from pseudoreactions (if any):
for i = 1:length(model.rxns)
    if endsWith(model.rxnNames{i},' pseudoreaction')
        model.grRules{i}      = '';
        model.rxnGeneMat(i,:) = zeros(1,length(model.genes));
    end
end

%Swap direction of reactions that are defined to only carry negative flux
to_swap=model.lb < 0 & model.ub == 0;
model.S(:,to_swap)=-model.S(:,to_swap);
model.ub(to_swap)=-model.lb(to_swap);
model.lb(to_swap)=0;

%Delete blocked rxns (LB = UB = 0):
to_remove = logical((model.lb == 0).*(model.ub == 0));
model     = removeReactions(model,model.rxns(to_remove),true,true,true);

%Correct rev vector: true if LB < 0 & UB > 0, or it is an exchange reaction:
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || sum(model.S(:,i) ~= 0) == 1
        model.rev(i) = true;
    end
end

%Make irreversible model (appends _REV to reaction IDs to indicate reverse
%reactions)
model=convertToIrrev(model);

%Expand model, to separate isoenzymes (appends _EXP_* to reaction IDs to
%indicate duplication)
model=expandModel(model);

%Sort reactions, so that reversible and isoenzymic reactions are kept near
model=sortIdentifiers(model);

%Remove simple substrates
reducedS     = model.S; %make a reduced S-matrix, with zero entries
metsToIgnore = contains(model.metNames,{'H2O','H+','sodium','phosphate',...
    'diphosphate','H2O','H+','iron(2+)','iron(3+)','hydrogen peroxide',...
    'ammonium','thiosulfate','hydrogen sulfide','potassium','sodium',...
    'sulphate','chloride','Mg(2+)','Mn(2+)','Zn(2+)','Ca(2+)','Cu2(+)'});
reducedS(metsToIgnore,:) = 0;

%Remove currency metabolites
fid         = fopen('currencyMets.tsv');
fileContent = textscan(fid,'%q %q','Delimiter','\t');
fclose(fid);
currMets.sub = fileContent{1};
currMets.pro = fileContent{2};
for i=1:numel(currMets.sub)
    idS = contains(model.metNames,currMets.sub{i});
    idP = contains(model.metNames,currMets.pro{i});
    [~,hasS] = find(reducedS(idS,:));
    [~,hasP] = find(reducedS(idP,:));
    hasBoth  = intersect(hasS,hasP);
    reducedS([idS,idP],hasBoth) = 0;
end



%Make ec-extension structure
emptyCell    = cell(numel(model.rxns)*5,1);
emptyVect    = zeros(numel(model.rxns)*5,1);
ec.rxns      = emptyCell;
ec.eccodes   = emptyCell;
ec.enzyme    = emptyCell; % The enzyme matching kcat value
ec.metNames  = emptyCell;
ec.metSmiles = emptyCell; %
ec.reverse   = emptyVect;
ec.kcat      = emptyVect;
ec.source    = emptyCell;
ec.notes     = emptyCell;
ec.sequence  = emptyCell; %
%Not duplicating unnecessary

ec.uniprot   = cell(numel(model.genes),1); % Could also be in metMiriams, but there it is more difficult to access
ec.mw        = zeros(numel(model.genes),1);
ec.gene      = model.genes; % Still good to keep it separate. Or should it contain all genes from the genome?
ec.sequence  = cell(numel(model.genes),1);
ec.metNames  = unique(ec.metNames);
ec.metSmiles = cell(numel(model.mets),1); % Could also be in model.metSmiles

%Extract model provided SMILES
%TODO: include a metSmiles field in RAVEN (same name as in COBRA), stored
%as annotation note in the SBML file
%TODO: loading of external SMILES database
if isfield(model,'metSmiles')
    metSmileDB   = model.metSmiles;
end

%Extract model provided Uniprot
%TODO: also add amino acid sequence and MW
% if isfield(model,'geneMiriams')
%     uniprotDB = extractMiriam(model.geneMiriams,'uniprot');
% end
%TODO: loading of external Uniprot database
%TODO: call downloadKEGG from updateDabases() to get the required data from
%KEGG. Currently manually load ProtDatabase.mat
uniprotDB.id    = kegg(:,1);
uniprotDB.gene  = kegg(:,3);
uniprotDB.mw    = cell2mat(kegg(:,5));
uniprotDB.seq   = kegg(:,7);

%Make protein metabolites
[Lia,Locb]      = ismember(model.genes,uniprotDB.gene);
ec.gene         = model.genes(Lia);
ec.uniprot      = uniprotDB.id(Locb);
ec.mw           = uniprotDB.mw(Locb);
ec.sequence     = uniprotDB.seq(Locb);
%Extract model provided eccodes
%TODO: loading of external ec-code database, possibly from uniprot, which
%should directly be parsed to the model, so that the eccodeDB entries match
%model.rxns.
if isfield(model,'eccodes')
    eccodeDB = model.eccodes;
end

%     %Leave out simple molecules (H+, H2O, Pi)
%     %Should be expanded to include any feasible option
%     simpleMetsDB = {'H+','hydrogen','H','H(+)','proton',...
%         'water','H2O',...
%         'Pi','PPi','phosphate'};



%Only parse rxns associated to genes
rxnWithGene = find(sum(model.rxnGeneMat,2));
ecCount=1;
for r=rxnWithGene' % model.rxns(r)
    %Parse through substrates
    substrates = find(reducedS(:,r)<0);
    substrateName = model.metNames(substrates);
    
    %Leave out simple molecules (H+, H2O, Pi)
    %Should be expanded to include any feasible option
%     simpleMets = matches(substrateName,simpleMetsDB); % matches() introduced R2019b, code alternative?
%     substrates(simpleMets)    = '';
%     substrateName(simpleMets) = '';
%     
    %TODO: Leave out currency metabolites if there is > 1 substrate present, and
    %the currency metabolite has its pair as product (not to affect e.g.
    %NAD biosynthesis)
%     if numel(substrates)>1
%         currMets1 = {'ATP','adenosine triphosphate','ADP','adenosine diphosphate',...
%             'NADH','NAD','NADPH','NAD'
%             'ADP','adenosine monophosphate','hydrogen','H','H(+)','proton','',...
%             'water','H2O',...
%             'Pi','phosphate'}); % matches() introduced R2019b, code alternative?
%         currMetsPairs = matches(substrateName,); % matches() introduced R2019b, code alternative?
%         substrates(currMets)    = '';
%         substrateName(currMets) = '';
%     end

    %Loop through all substrates
    for s=substrates' % s = model.mets(s)
        e.rxns                  = model.rxns(r);
        e.eccodes               = model.eccodes(r);
        if ~isempty(regexp(model.rxns{r},'_REV($|_)'))
            e.reverse           = 1;
        end
        e.metNames              = model.metNames(s);
        %e.metSmiles             = model.metSmiles(s);
                
        %Fill in the gene information
        rxnGenes = find(model.rxnGeneMat(r,:)); 
        for g=rxnGenes % g = model.genes(g)
            ec.rxns(ecCount)        = e.rxns;
            %ec.eccodes(ecCount)     = e.eccodes;
            ec.reverse(ecCount)     = e.reverse;
            ec.metNames(ecCount)    = e.metNames;
            %ec.metSmiles(ecCount)   = e.metSmiles;

            %ec.enzyme(ecCount)    = model.genes(g);
            %ec.uniprot(ecCount) = uniprotDB.id(g);
            %ec.mw(ecCount)      = uniprotDB.mw(g);
            %ec.sequence(ecCount)= uniprotDB.seq(g);

            ecCount=ecCount+1;
        end
    end
end
%Remove empty fields
for f=fieldnames(ec)'
    ec.(f{:})(ecCount-1:end)=[];
end
model.ec=ec;
end
