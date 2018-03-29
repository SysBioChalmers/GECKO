%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function model = fitGAM(model,biomassRxn)
% 
% Gets a GEM mat structure and refits the growth associated manteinance
% cost in the biomass rxn for reaching a given treshold in growth rate on
% minimal glucose media.
%
% INPUT:    model       GEM file (.mat format)
%           biomassRxn  String containing the Biomass pseudoreaction name
%           Ptot        Experimentally measured total protein content in
%                       the cell.
% OUTPUTS:  model       The model with modified biomass pseudoreaction
%                       stoichiometry.
            
% 
% Ivan Domenzain. Last edited: 2017-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = fitGAM(model,biomassRxn,gR_exp,Ptot)
    %Biopos           = zeros(1,4);
    tempModel        = model;
    %Find biomass pseudoRxn in model
    pos              = find(contains(model.rxnNames,biomassRxn));
    disp(pos)
    %Set biomass as objective function
    tempModel.c      = zeros(length(model.c),1);
    tempModel.c(pos) = 1;
    solution         = solveLP(tempModel,1);
    GamValues        = [];
    errors           = [];
    %Gets substrate and products in the biomass pseudoRxn
    subsIndexes  = find(tempModel.S(:,pos)<0);
    prodsIndexes = find(tempModel.S(:,pos)>0);
    substrates   = tempModel.metNames(subsIndexes);
    products     = tempModel.metNames(prodsIndexes);
    % Correct protein requeriments according to experimentally measured
    % total protein content
    %if nargin > 3
        Protpos = find(contains(substrates,'protein'));
        if ~isempty(Protpos)
            model.S(subsIndexes(Protpos),pos) = -1*Ptot;
        end
    %end
    % Find original GAM and set a consistent biomass stoichiometry
    Biopos(1) = subsIndexes(find(contains(substrates,'ATP')));
    Biopos(2) = subsIndexes(find(contains(substrates,'H2O')));
    Biopos(3) = prodsIndexes(find(contains(products,'ADP')));
    Biopos(4) = prodsIndexes(find(contains(products,'phosphate')));
    GAM       = abs(model.S(Biopos(1),pos));
    disp(GAM)
    tempModel = setBiomassStoich(tempModel,pos,Biopos,GAM);
    
    for i=1:100
        newGAM    = GAM-((i-1)*GAM/100);
        tempModel = setBiomassStoich(tempModel,pos,Biopos,newGAM);
        solution  = solveLP(tempModel,1);
        if ~isempty(solution.f)
            error     = (gR_exp-solution.x(pos))*100/gR_exp;
            GamValues = [GamValues; newGAM];
            errors    = [errors; error]; 
        end
    end
    [minError, minIndx] = min(errors);
    figure
    plot(GamValues,errors,'LineWidth',5)
    title('GAM fitting for growth on glucose minimal media')
    xlabel('Growth associated manteinance [mmol/gDwh]')
    ylabel('Absolute relative error [%]')

end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function model = setBiomassStoich(model,RxnPos,metsPos,GAM)
    for i=1:length(metsPos)
        stoichSign                 = sign(model.S(metsPos(i),RxnPos));
        model.S(metsPos(i),RxnPos) = stoichSign*GAM;
    end
end
    
    
    
    
