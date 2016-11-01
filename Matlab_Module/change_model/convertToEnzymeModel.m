%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eModel = convertToEnzymeModel(model,uniprots,kcats)
% Converts standard GEM to GEM accounting for enzymes as pseudo
% metabolites, with -(MW/kcat) as the corresponding stoich. coeffs.
%
% INPUT:
% model             The GEM structure (1x1 struct)
% uniprots          uniprot codes of each reaction
% kcats             kcats for each enzyme/reaction
%
% OUTPUT:
% eModel            Modified GEM structure (1x1 struct)
% 
% Cheng Zhang. Last edited: 2016-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eModel = convertToEnzymeModel(irrevModel,uniprots,kcats)

eModel  = irrevModel;
enzymes = cell(5000,1);
[m,n]   = size(uniprots);
y       = 0;
for i = 1:m
    rxnName = irrevModel.rxns{i};
    x = 0;
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
            eModel = addArmReaction(eModel,rxnName);
        end
        x = 0;
        %For each enzyme adds a new parallel rxn with that enzyme as pseudo metabolite.
        for j = 1:n
            if ~isempty(uniprots{i,j}) && kcats(i,j) > 0 % if kcat=0 no mets will be added
                x       = x+1;
                newName = [rxnName 'No' num2str(x)];
                kvalues = ones(size(uniprots{i,j}))*kcats(i,j);
                newMets = cell(size(uniprots{i,j}));
                for k = 1:length(newMets)
                    newMets{k} = ['prot_' uniprots{i,j}{k}];
                end 
                eModel = addEnzymesToRxn(eModel,kvalues,rxnName,newMets,newName); 
            end
        end
        eModel = removeRxns(eModel,{rxnName});  %Remove the original rxn
    end
end

%Eliminate repeated uniprots from enzymes vector:
enzymes(y+1:end) = [];
enzymes          = unique(enzymes)';

%Load databases:
cd ../../Databases
data      = load('ProtDatabase.mat');
swissprot = data.swissprot;
kegg      = data.kegg;
cd ../Matlab_Module/change_model

%Create additional fields in model:
eModel.enzymes   = cell(0,1);
eModel.genes2    = cell(0,1);
eModel.geneNames = cell(0,1);
eModel.MWs       = zeros(0,1);
eModel.sequence  = cell(0,1);
eModel.pathways  = cell(0,1);
for i = 1:length(enzymes)
    eModel = addProtein(eModel,enzymes{i},kegg,swissprot);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%