%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eModel = readKcatData(model_data,kcats)
% Reads the output of the getBRENDAdata module, and creates the modified
% enzyme model, with additional metabolites (the enzymes) and reactions
% (for all isoenzymes and also the enzyme exchange reactions).
%
% INPUT:
% model_data        model and EC numbers and substrates/products from each
%                   reaction (output from "getECnumbers.m")
% kcats             kcats for each reaction/enzyme (output from
%                   "matchKcats.m")
%
% OUTPUT:
% eModel            modified model accounting for enzymes
% 
% Cheng Zhang. Last edited: 2015-12-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eModel = readKcatData(model_data,kcats)

%Get kcat value for both directions:
Fkcat = kcats.forw.kcats;
Bkcat = kcats.back.kcats;
rev   = boolean(model_data.model.rev);
kcats = [Fkcat;Bkcat(rev,:)];
%Get model:
model = model_data.model;
%Predefine ECnumber and uniprots for enzyme model:
ECnumbers = [model_data.EC_numbers ; model_data.EC_numbers(rev,:)];
uniprots  = [model_data.uniprots   ; model_data.uniprots(rev,:)  ];
%Convert to irreversible model with RAVEN function (will split in 2 any reversible rxn):
model = convertToIrrev(model);
%Convert original model to enzyme model according to uniprots and kcats:
eModel = convertToEnzymeModel(model,uniprots,kcats);

%Leave all UB = +Inf:
eModel.ub(eModel.ub == 1000) = +Inf;

%Add backup data to model:
eModel.rxnUniprots  = uniprots;
eModel.rxnECnumbers = ECnumbers;
eModel.rxnKcats     = kcats;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%