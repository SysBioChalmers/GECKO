function eModel = convertToEnzymeModel(irrevModel,Genes,uniprots,kcats)
%convertToEnzymeModel
%
% eModel = convertToEnzymeModel(model,Genes,uniprots,kcats)
% Converts standard GEM to GEM accounting for enzymes as pseudo
% metabolites, with -(1/kcat) as the corresponding stoich. coeffs.
%
% INPUT:
% model             The GEM structure (1x1 struct)
% Genes             Genes that were matched to uniprot codes for each Rxn
% uniprots          uniprot codes of each reaction
% kcats             kcats for each enzyme/reaction
%
% OUTPUT:
% eModel            Modified GEM structure (1x1 struct)
% 
% Cheng Zhang.          Last edited: 2018-05-24
% Benjamin J. Sanchez.  Last edited: 2018-11-11
% Ivan Domenzain.       Last edited: 2020-02-07

fprintf('Adding enzymes to reactions...')

%Load databases:
data      = load('../../databases/ProtDatabase.mat');
swissprot = data.swissprot;
kegg      = data.kegg;
[nrows,~] = size(kegg);
eModel  = irrevModel;
enzymes = cell(5000,1);
[m,n]   = size(uniprots);
y       = 0;

for i = 1:m
    rxnID = irrevModel.rxns{i};
    x     = 0;
    for j = 1:n
        %Update vector enzymes and count the number of isozymes (x):
        if ~isempty(uniprots{i,j}) && kcats(i,j) > 0 % if kcat=0 no mets will be added
            uniprots{i,j} = strsplit(uniprots{i,j},' ');
            for k = 1:length(uniprots{i,j})
                y          = y+1;
                enzymes{y} = uniprots{i,j}{k};
            end
            x = x+1;
        end
    end
        
    if x > 0
        %>1 enzyme: Will include an "arm reaction" for controlling the
        %total flux through the system of parallel rxns.
        if x > 1
            eModel  = addArmReaction(eModel,rxnID);
        end
        x = 0;
        %For each enzyme adds a new parallel rxn with that enzyme as pseudo metabolite.
        for j = 1:n
            if ~isempty(uniprots{i,j}) && kcats(i,j) > 0 % if kcat=0 no mets will be added
                x         = x+1;
                newID     = [rxnID 'No' num2str(x)];
                newName   = [irrevModel.rxnNames{i} ' (No' num2str(x) ')'];
                kvalues   = ones(size(uniprots{i,j}))*kcats(i,j);
                newMets   = cell(size(uniprots{i,j}));
                for k = 1:length(newMets)
                    newMets{k} = ['prot_' uniprots{i,j}{k}];
                end
                grRule = strrep(Genes{i,j},' ',' and ');
                eModel = addEnzymesToRxn(eModel,kvalues,rxnID,newMets,{newID,newName},grRule); 
            end
        end
        eModel = removeReactions(eModel,{rxnID});  %Remove the original rxn
    end
    if rem(i,100) == 0 || i == m 
        fprintf('.')
    end
end

%Eliminate repeated uniprots from enzymes vector:
enzymes(y+1:end) = [];
enzymes          = unique(enzymes)';

%Create additional fields in model:
eModel.enzymes   = cell(0,1);
eModel.enzGenes  = cell(0,1);
eModel.enzNames  = cell(0,1);
eModel.MWs       = zeros(0,1);
eModel.sequences = cell(0,1);
eModel.pathways  = cell(0,1);
for i = 1:length(enzymes)
    if nrows>1
        eModel = addProtein(eModel,enzymes{i},kegg,swissprot);
    else
        eModel = addProtein(eModel,enzymes{i},[],swissprot);
    end
end

fprintf(' Done!\n')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
