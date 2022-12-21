%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = removeIncorrectPathways(model)

fprintf('Removing incorrect pathways...')

%Construct vector of non repeated pathways:
for i = 1:length(model.enzymes)
    pathway_i = model.pathways{i};
    pathway_new = '';
    pos       = strfind(pathway_i,'sce0');
    for j = 1:length(pos)
        if j == length(pos)
            pos_end = length(pathway_i);
        else
            pos_end = pos(j+1)-2;
        end
        pathway = pathway_i(pos(j):pos_end);
        %Detect weird pathways:
        weird = strcmp(pathway,'sce00591  Linoleic acid metabolism') || ...
                strcmp(pathway,'sce00590  Arachidonic acid metabolism') || ...
                strcmp(pathway,'sce00592  alpha-Linolenic acid metabolism') || ...
                strcmp(pathway,'sce00565  Ether lipid metabolism') || ...
                strcmp(pathway,'sce00460  Cyanoamino acid metabolism') || ...
                strcmp(pathway,'sce00680  Methane metabolism');
            
        if ~weird
            pathway_new = [pathway_new pathway ' '];
        end
        
        model.pathways{i} = pathway_new(1:end-1);
    end
    
    if strcmp(model.pathways{i},'') || strcmp(model.pathways{i},'-')
        model.pathways{i} = 'sce00NNN  None';
    end
end

fprintf(' Done!\n')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%