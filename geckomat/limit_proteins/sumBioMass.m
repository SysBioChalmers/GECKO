%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,P,C,R,D,L] = sumBioMass(model)
% Calculates breakdown of biomass:
% X -> Biomass fraction without lipids [g/gDW]
% P -> Protein fraction [g/gDW]
% C -> Carbohydrate fraction [g/gDW]
% R -> RNA fraction [g/gDW]
% D -> DNA fraction [g/gDW]
% L -> Lipid fraction [g/gDW]
% Function adapted from SLIMEr: https://github.com/SysBioChalmers/SLIMEr
%
% Benjamin Sanchez. Last update: 2018-08-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,P,C,R,D,L] = sumBioMass(model)

%Components of biomass:
%        id         MW [g/mol]  class     name
comps = {'s_0404'	89.09       'P'     % A     Alanine         ala
         's_0542'	121.16      'P'     % C     Cysteine        cys
         's_0432'   133.11      'P'     % D     Aspartic acid   asp
         's_0748'   147.13      'P'     % E     Glutamic acid   glu
         's_1314'   165.19      'P'     % F     Phenylalanine   phe
         's_0757'   75.07       'P'     % G     Glycine         gly
         's_0832'   155.15      'P'     % H     Histidine       his
         's_0847'   131.17      'P'     % I     Isoleucine      ile
         's_1099'   146.19      'P'     % K     Lysine          lys
         's_1077'   131.17      'P'     % L     Leucine         leu
         's_1148'   149.21      'P'     % M     Methionine      met
         's_0430'   132.12      'P'     % N     Asparagine      asn
         's_1379'   115.13      'P'     % P     Proline         pro
         's_0747'   146.14      'P'     % Q     Glutamine       gln
         's_0428'   174.2       'P'     % R     Arginine        arg
         's_1428'   105.09      'P'     % S     Serine          ser
         's_1491'   119.12      'P'     % T     Threonine       thr
         's_1561'   117.15      'P'     % V     Valine          val
         's_1527'   204.23      'P'     % W     Tryptophan      trp
         's_1533'   181.19      'P'     % Y     Tyrosine        tyr
         's_0001'	180.16      'C'     % (1->3)-beta-D-glucan
         's_0004'	180.16      'C'     % (1->6)-beta-D-glucan
         's_0509'   221.21      'C'     % chitin
         's_0773'   180.16      'C'     % glycogen
         's_1107'   180.16      'C'     % mannan
         's_1520'   342.296 	'C'     % trehalose
         's_0423'   347.22      'R'     % AMP
         's_0526'   323.2       'R'     % CMP
         's_0782'   363.22      'R'     % GMP
         's_1545'   324.18      'R'     % UMP
         's_0584'   331.22      'D'     % dAMP
         's_0589'   307.2       'D'     % dCMP
         's_0615'   345.21      'D'     % dGMP
         's_0649'   322.21      'D'     % dTMP
         's_3714'   852.83      'N'     % heme a
         's_1405'   376.36      'N'     % riboflavin
         's_1467'   96.06       'N'};   % sulphate

%Get main fractions:
[P,X] = getFraction(model,comps,'P',0);
[C,X] = getFraction(model,comps,'C',X);
[R,X] = getFraction(model,comps,'R',X);
[D,X] = getFraction(model,comps,'D',X);
[L,X] = getFraction(model,comps,'L',X);

%Add up any remaining components:
bioPos = strcmp(model.rxns,'r_4041');
for i = 1:length(model.mets)
    pos = strcmp(comps(:,1),model.mets{i});
    if sum(pos) == 1
        abundance = -model.S(i,bioPos)*comps{pos,2}/1000;
        X         = X + abundance;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,X] = getFraction(model,comps,compType,X)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');
rxnName = strrep(rxnName,'L','lipid backbone');

%Add up fraction:
fractionPos = strcmp(model.rxnNames,rxnName);
if contains(rxnName,'lipid')
    subs = model.S(:,fractionPos) < 0;        %substrates in pseudo-rxn
    F    = -sum(model.S(subs,fractionPos));   %g/gDW
else
    comps = comps(strcmp(comps(:,3),compType),:);
    F = 0;
    %Add up all components:
    for i = 1:length(model.mets)
        pos = strcmp(comps(:,1),model.mets{i});
        if sum(pos) == 1
            abundance = -model.S(i,fractionPos)*comps{pos,2}/1000;
            F         = F + abundance;
        end
    end
end
X = X + F;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
