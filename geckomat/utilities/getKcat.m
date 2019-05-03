function [kcat,rxnIdx,rxnName,MW] = getKcat(ecModel,protein)
% getKcat
%
%   This function can quickly get protein related data from an ecModel,
%   which can be useful when curating an ecModel. By providing a protein
%   ID, the ecModel is queried to find reactions involving that protein.
%   The output includes kcat values (in sec-1, not corrected for the fact
%   that stoichiometric coefficients of the substrate can influence
%   assigning the kcat by matchKcats()), index and names of reactions
%   involving the protein, in addition to the proteins molecular weight.
%   This func
%
%   Input:
%   ecModel         enzyme-constrained model
%   protein         string specifying a UniProt protein ID
%
%   Output:
%   kcat            kcat values of all reactions involving the protein.
%                   Note: kcat values assigned by matchKcats() might have
%                   been corrected for the stoichiometric coefficients of
%                   the substrates in the reaction
%   rxnIdx          indices of all reactions involving the protein
%   rxnName         names of all reactions involving the protein
%   MW              molecular weight of the protein, in kDa
%
% usage: [kcat,rxnIdx,rxnName,MW] = getKcat(ecModel,protein)
%
% Eduard Kerkhoven  Last edited: 2018-12-19


pIdx    = strcmpi(ecModel.metNames,['prot_' protein]);
rxnIdx  = find(ecModel.S(pIdx,:) < 0);
kcat    = -1./ecModel.S(pIdx,rxnIdx)/3600;
rxnName = ecModel.rxnNames(rxnIdx);

MW = strcmpi(model.enzymes,protein);
MW = model.MWs(MW);
end
