function [model,pos] = changeMedia(model,c_source,media,anox,flux)
%changeMedia
%
% Function that modifies the ecModel and makes it suitable for batch growth
% simulations on different carbon sources.
%
% INPUT:
%   - model:  An enzyme constrained model
%   - media:  Media type ('YEP' for complex, 
%                         'MAA' minimal with Aminoacids,
%                         'MIN' for minimal media,
%                         'MIN+His' for minimal media with his
%                         'MIN+Arg' for minimal media with arg
%                         'MIN+Citrate' for minimal media with Citrate)
%   - anox:   (optional) TRUE if anaerobic conditions are desired, DEFAULT=
%             FALSE
%   - flux:   (Optional) A cell array with measured uptake fluxes in mmol/gDwh
%
% OUTPUT:
%   - model: The ECmodel with the specified medium constraints
%
% Ivan Domenzain        2020-01-17
% Feiran Li             2020-11-10

if nargin<4
    anox = false;
    if nargin<3
        media = 'MIN';
    end
end
% Give the carbon source (c_source) input variable with the following
% format: c_source  = 'D-glucose exchange'
c_source = [c_source,' exchange'];
% First block all uptake reactions
[~,exchange]  = getExchangeRxns(model,'out');
% But make exception for protein pool exchange
[~,idx] = ismember('prot_pool_exchange',model.rxnNames);
exchange = setdiff(exchange,idx);
model = setParam(model,'lb',exchange,0);
% Also block uptake if reactions are defined in the opposite direction
[~,exchange]  = getExchangeRxns(model,'in');
model = setParam(model,'ub',exchange,0);

pos = getComponentIndexes(model,c_source);


% %For growth on fructose and mannose the transport takes place in a passive
% %way. [Boles & Hollenberg, 2006]
% if strcmp(c_source,'D-fructose exchange')
%     model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1134')) = 0;
%     model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1134')) = 0;
% elseif strcmp(c_source,'D-mannose exchange')
%     model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1139')) = 0;
%     model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1139')) = 0;
% end
%The media will define which rxns to fix:
switch media
    case 'YEP'
        N = 25;     %Aminoacids + Nucleotides
    case 'MAA'
        N = 21;     %Aminoacids
    case 'MIN'
        N = 1;      %Only the carbon source
    case 'MIN+His'
        N = 1;      %Only the carbon source
        model = setParam(model,'lb','r_1893',-0.08);	%Histidine exchange
    case 'MIN+Arg'
        N = 1;      %Only the carbon source
        model = setParam(model,'lb','r_1879',-0.08);	%L-arginine exchange
    case 'MIN+Citrate'
        N = 1;      %Only the carbon source
        model = setParam(model,'lb','r_1687',-0.08);	%citrate exchange   
end
%LB parameter (manually optimized for glucose on Min+AA):
b = -0.08;
%LB parameter (manually optimized for glucose complex media):
c = -2;
%Define fluxes in case of ec model:
if nargin < 5   %Limited protein
    if N>1
       flux    = b*ones(1,N);
       if N>21
           flux(22:25) = c;
       end
    end
    flux(1) = -1000;
end
%Fix values as LBs:
for i = 1:N
    model.lb(pos(i)) = flux(i);
end
%Allow uptake of essential components
model = setParam(model, 'lb', 'r_1654', -Inf); % 'ammonium exchange';
model = setParam(model, 'lb', 'r_2100', -Inf); % 'water exchange' ;
model = setParam(model, 'lb', 'r_1861', -Inf); % 'iron(2+) exchange';
model = setParam(model, 'lb', 'r_1992', -Inf); % 'oxygen exchange';
model = setParam(model, 'lb', 'r_2005', -Inf); % 'phosphate exchange';
model = setParam(model, 'lb', 'r_2060', -Inf); % 'sulphate exchange';
model = setParam(model, 'lb', 'r_1832', -Inf); % 'H+ exchange' ;
model = setParam(model, 'lb', 'r_4593', -Inf); % 'chloride exchange' ;
model = setParam(model, 'lb', 'r_4595', -Inf); % Mn(2+) exchange
model = setParam(model, 'lb', 'r_4596', -Inf); % Zn(2+ exchange
model = setParam(model, 'lb', 'r_4597', -Inf); % Mg(2+) exchange
model = setParam(model, 'lb', 'r_2049', -Inf); % sodium exchange
model = setParam(model, 'lb', 'r_4594', -Inf); % Cu(2+) exchange
model = setParam(model, 'lb', 'r_4600', -Inf); % Ca(2+) exchange
model = setParam(model, 'lb', 'r_2020', -Inf); % potassium exchange
%Block some production fluxes
model = setParam(model, 'ub', 'r_1663', 0); % bicarbonate exchange
model = setParam(model, 'ub', 'r_4062', 0); % lipid backbone exchange
model = setParam(model, 'ub', 'r_4064', 0); % lipid chain exchange

%Allow biomass production
model = setParam(model, 'ub', 'r_2111', Inf); % growth

% change aeerobic or anaerobic
if strcmp(anox,'anaerobic')
    model = anaerobicModel_GECKO(model);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getComponentIndexes(model,c_source)
    pos(1)  = find(strcmpi(model.rxnNames,c_source));
    pos(2)  = find(strcmpi(model.rxnNames,'L-alanine exchange'));
    pos(3)  = find(strcmpi(model.rxnNames,'L-arginine exchange'));
    pos(4)  = find(strcmpi(model.rxnNames,'L-asparagine exchange'));
    pos(5)  = find(strcmpi(model.rxnNames,'L-aspartate exchange'));
    pos(6)  = find(strcmpi(model.rxnNames,'L-cysteine exchange'));
    pos(7)  = find(strcmpi(model.rxnNames,'L-glutamine exchange'));
    pos(8)  = find(strcmpi(model.rxnNames,'L-glutamate exchange'));
    pos(9)  = find(strcmpi(model.rxnNames,'L-glycine exchange'));
    pos(10) = find(strcmpi(model.rxnNames,'L-histidine exchange'));
    pos(11) = find(strcmpi(model.rxnNames,'L-isoleucine exchange'));
    pos(12) = find(strcmpi(model.rxnNames,'L-leucine exchange'));
    pos(13) = find(strcmpi(model.rxnNames,'L-lysine exchange'));
    pos(14) = find(strcmpi(model.rxnNames,'L-methionine exchange'));
    pos(15) = find(strcmpi(model.rxnNames,'L-phenylalanine exchange'));
    pos(16) = find(strcmpi(model.rxnNames,'L-proline exchange'));
    pos(17) = find(strcmpi(model.rxnNames,'L-serine exchange'));
    pos(18) = find(strcmpi(model.rxnNames,'L-threonine exchange'));
    pos(19) = find(strcmpi(model.rxnNames,'L-tryptophan exchange'));
    pos(20) = find(strcmpi(model.rxnNames,'L-tyrosine exchange'));
    pos(21) = find(strcmpi(model.rxnNames,'L-valine exchange'));
    pos(22) = find(strcmpi(model.rxnNames,'2''-deoxyadenosine exchange'));
    pos(23) = find(strcmpi(model.rxnNames,'2''-deoxyguanosine exchange'));
    pos(24) = find(strcmpi(model.rxnNames,'thymidine exchange'));
    pos(25) = find(strcmpi(model.rxnNames,'deoxycytidine exchange'));
    pos(26) = find(strcmpi(model.rxnNames,'D-glucose exchange'));
    pos = transpose(pos);
end
