%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeMedia(model,media,flux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,pos] = changeMedia_batch(model,c_source,media,b,flux)

% Block O2 and glucose production (for avoiding multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;

%Find substrate production rxn and block it:
pos_rev = strcmpi(model.rxnNames,c_source(1:strfind(c_source,...
                                            ' (reversible)')-1));
model.ub(pos_rev) = 0;
%The media will define which rxns to fix:
if strcmpi(media,'YEP')
    N = 25;     %Aminoacids + Nucleotides
elseif strcmpi(media,'MAA')
    N = 21;     %Aminoacids
elseif strcmpi(media,'Min')
    N = 1;      %Only the carbon source
end
if nargin < 4   
    %UB parameter (manually optimized for glucose on Min+AA):
    b = 0.07;
    %UB parameter (manually optimized for glucose complex media):
    c = 2;
end

%Define fluxes in case of ec model:
if nargin < 5   %Limited protein    
    if N>1
       flux    = b*ones(1,N);
       if N>21
           flux(22:25) = c;
       end
    end
    flux(1) = 1000;
end

%Consumption positions:
pos(1)  = find(strcmpi(model.rxnNames,c_source));
pos(2)  = find(strcmpi(model.rxnNames,'L-alanine exchange (reversible)'));
pos(3)  = find(strcmpi(model.rxnNames,'L-arginine exchange (reversible)'));
pos(4)  = find(strcmpi(model.rxnNames,'L-asparagine exchange (reversible)'));
pos(5)  = find(strcmpi(model.rxnNames,'L-aspartate exchange (reversible)'));
pos(6)  = find(strcmpi(model.rxnNames,'L-cysteine exchange (reversible)'));
pos(7)  = find(strcmpi(model.rxnNames,'L-glutamine exchange (reversible)'));
pos(8)  = find(strcmpi(model.rxnNames,'L-glutamate exchange (reversible)'));
pos(9)  = find(strcmpi(model.rxnNames,'glycine exchange (reversible)'));
pos(10) = find(strcmpi(model.rxnNames,'L-histidine exchange (reversible)'));
pos(11) = find(strcmpi(model.rxnNames,'L-isoleucine exchange (reversible)'));
pos(12) = find(strcmpi(model.rxnNames,'L-leucine exchange (reversible)'));
pos(13) = find(strcmpi(model.rxnNames,'L-lysine exchange (reversible)'));
pos(14) = find(strcmpi(model.rxnNames,'L-methionine exchange (reversible)'));
pos(15) = find(strcmpi(model.rxnNames,'L-phenylalanine exchange (reversible)'));
pos(16) = find(strcmpi(model.rxnNames,'L-proline exchange (reversible)'));
pos(17) = find(strcmpi(model.rxnNames,'L-serine exchange (reversible)'));
pos(18) = find(strcmpi(model.rxnNames,'L-threonine exchange (reversible)'));
pos(19) = find(strcmpi(model.rxnNames,'L-tryptophan exchange (reversible)'));
pos(20) = find(strcmpi(model.rxnNames,'L-tyrosine exchange (reversible)'));
pos(21) = find(strcmpi(model.rxnNames,'L-valine exchange (reversible)'));
pos(22) = find(strcmpi(model.rxnNames,'2''-deoxyadenosine exchange (reversible)'));
pos(23) = find(strcmpi(model.rxnNames,'2''-deoxyguanosine exchange (reversible)'));
pos(24) = find(strcmpi(model.rxnNames,'thymidine exchange (reversible)'));
pos(25) = find(strcmpi(model.rxnNames,'deoxycytidine exchange (reversible)'));

%Fix values as UBs:
for i = 1:N
    model.ub(pos(i)) = flux(i);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%