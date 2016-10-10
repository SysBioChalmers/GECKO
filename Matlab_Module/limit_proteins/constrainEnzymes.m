%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = constrainEnzymes(model,Ptotal,sigma,pIDs,data)
% 
%
% Benjamín J. Sánchez. Last edited: 2016-08-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = constrainEnzymes(model,Ptotal,sigma,pIDs,data)

%Current values:
f       = 0.4462; %Yeast 7.6 (all enzymes) [g(Pmodel)/g(Ptot)]
options = [2;3;0.7;1;2;2;0.4266];

if options(4) == 1
    Pbase = options(7);
else
    Pbase = Ptotal;
end

%No UB will be changed if no data is available -> pool = all enzymes(FBAwMC)
if nargin == 3
    pIDs = {};
    data = [];
end

%Assign concentrations as UBs [mmol/gDW]:
model.concs = nan(size(model.enzymes));      %OBS: min value is zero!!
disp('Matching data to enzymes in model...')
for i = 1:length(model.enzymes)
    match = false;
    for j = 1:length(pIDs)
        if strcmpi(pIDs{j},model.enzymes{i}) && ~match
            model.concs(i) = data(j)*model.MWs(i); %g/gDW
            rxn_name       = ['prot_' model.enzymes{i} '_exchange'];
            pos            = strcmpi(rxn_name,model.rxns);
            model.ub(pos)  = data(j);
            match          = true;
        end
    end
end

%Count mass of non-measured enzymes:
measured       = ~isnan(model.concs);
concs_measured = model.concs(measured);
Pmeasured      = sum(concs_measured);

if Pmeasured > 0
    %Calculate fraction of non measured proteins in model out of remaining mass:
    fm     = Pmeasured/Ptotal;
    [fn,~] = measureAbundance(model.enzymes(~measured),'prot_abundance.txt');
    f      = fn/(1-fm);
    %Discount measured mass from global constrain:
    fs = (Ptotal - Pmeasured)/Pbase*f*sigma;
else
    fs = f*sigma;
end


%Constrain the rest of enzymes with the pool assumption:
model = constrainPool(model,~measured,1000);
model = fixedModifications(model,options);
model = changeProtein(model,Ptotal,fs);

%Display some metrics:
disp(['Total protein amount measured = '     num2str(Pmeasured)              ' g/gDW'])
disp(['Total enzymes measured = '            num2str(sum(measured))          ' enzymes'])
disp(['Enzymes in model with 0 g/gDW = '     num2str(sum(concs_measured==0)) ' enzymes'])
disp(['Total protein amount not measured = ' num2str(Pbase - Pmeasured)      ' g/gDW'])
disp(['Total enzymes not measured = '        num2str(sum(~measured))         ' enzymes'])
disp(['Total protein in model = '            num2str(Pbase)                  ' g/gDW'])


%Plot histogram:
figure
hist(concs_measured*1e3,10.^(-3:0.5:3))
set(gca,'xscale','log')
xlim([1e-3,1e3])
xlabel('Protein amount [mg/gDW]');
ylabel('Frequency');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%