%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,pIDs] = fixComplex(rxn,model,data,pIDs)

%Find subunits (name & stoichiometry) of the rxn's complex:
rxn_coeffs  = full(model.S(:,strcmpi(rxn,model.rxns)));
is_rxn      = rxn_coeffs ~= 0;
is_prot     = ~cellfun(@isempty,strfind(model.mets,'prot_'));
pos         = find(is_rxn.*is_prot);
pos_mets    = find(is_rxn.*~is_prot);
prot_names  = strrep(model.mets(pos),'prot_','');
prot_stoich = abs(rxn_coeffs(pos))/min(abs(rxn_coeffs(pos)));

%Go through all subunits and retrieve MW & experimental measurement:
prot_MWs    = zeros(size(prot_names));
prot_values = zeros(size(prot_names));
for i = 1:length(prot_names)
    %Warn if subunit is involved in another rxn in the model (not isozyme):
    rxn_pos = find(model.S(pos(i),:) < 0);
    if length(rxn_pos) > 1
        isozyme = true;
        for j = 1:length(rxn_pos)
            is_rxn_new   = full(model.S(:,rxn_pos(j))) ~= 0;
            pos_mets_new = find(is_rxn_new.*~is_prot);
            if ~isequal(pos_mets,pos_mets_new)
                isozyme = false;
            end
        end
        if ~isozyme && isempty(strfind(rxn,'r_1021'))
            disp(['WARNING: Protein ' prot_names{i} ' in ' num2str(length(rxn_pos)) ...
                  ' rxns: ' strjoin(model.rxns(rxn_pos)',' - ')])
        end
    end
    %Retrieve MW:
    enz_pos     = strcmpi(prot_names{i},model.enzymes);
    prot_MWs(i) = model.MWs(enz_pos);
    %Retrieve exp. measurement:
    data_pos    = strcmpi(prot_names{i},pIDs);
    if sum(data_pos) > 0
        prot_values(i) = data(data_pos);
    end
    %Display missing measurements:
    if prot_values(i) == 0
        disp(['WARNING: Protein ' prot_names{i} ' not present in ' ...
              'dataset (rxn ' rxn ')'])
    end
end

%Replace missing data to be proportional to the average measured subunit:
prot_prop = prot_values;
new_value = mean(prot_values(prot_values > 0)./prot_stoich(prot_values > 0));
for i = 1:length(prot_values)
    prot_prop(i) = prot_stoich(i)*new_value;
end

%Display increase in mass:
prev_total = sum(prot_values.*prot_MWs);
new_total  = sum(prot_prop.*prot_MWs);
disp(['Solving ' rxn ' issue: Inconsistent complex had ' num2str(prev_total) ...
      ' g/gDW -> new complex has ' num2str(new_total) ' g/gDW.'])
  
%Update exp. data:
for i = 1:length(prot_names)
    data_pos = strcmpi(prot_names{i},pIDs);
    if sum(data_pos) == 0
        data(length(data)+1) = prot_prop(i);
        pIDs{length(pIDs)+1} = prot_names{i};
    else
        data(data_pos) = prot_prop(i);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%