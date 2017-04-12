%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = manualModifications(model)
% 
%
% Benjamín J. Sánchez. Last edited: 2017-01-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = manualModifications(model)

%Read manual data:
cd ../../Databases
fID       = fopen('manual_data.txt');
data      = textscan(fID,'%s %s %s %s %f','delimiter','\t');
structure = data{2};
kcats     = data{5}.*3600;
data      = load('ProtDatabase.mat');
swissprot = data.swissprot;
kegg      = data.kegg;
fclose(fID);
cd ../Matlab_Module/change_model

%Construct curated complexes:
uniprots = cell(size(kcats));
stoich   = cell(size(kcats));
for i = 1:length(kcats)
    uniprots{i} = strsplit(structure{i},' + ');
    stoich{i}   = ones(size(uniprots{i}));
    %Separate complex strings in units and amount of each unit:
    for j = 1:length(uniprots{i})
        unit = uniprots{i}{j};
        pos  = strfind(unit,' ');
        if isempty(pos)
            stoich{i}(j)   = 1;
            uniprots{i}{j} = unit;
        else
            stoich{i}(j)   = str2double(unit(1:pos-1));
            uniprots{i}{j} = unit(pos+1:end);
        end
    end
end

for i = 1:length(model.rxns)
    %Find set of proteins present in rxn:
    S        = full(model.S);
    subs_pos = find(S(:,i) < 0);
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos);
    prot_set = cell(size(int_pos));
    MW_set   = 0;
    for j = 1:length(int_pos)
        met_name    = model.mets{int_pos(j)};
        prot_set{j} = met_name(6:end);
        MW_set      = MW_set + model.MWs(strcmp(model.enzymes,prot_set{j}));
    end
    %Find intersection with manual curated data:
    for j = 1:length(uniprots)
        int    = intersect(prot_set,uniprots{j});
        if length(int)/max(length(prot_set),length(uniprots{j})) > 0.50 % 50% match
            %Erase previous protein stoich. coeffs from rxn:
            for k = 1:length(prot_set)
                model.S(int_pos(k),i) = 0;
            end
            %If some proteins where not present previously, add them:
            newMets = uniprots{j};
            for k = 1:length(uniprots{j})
                if sum(strcmp(model.enzymes,uniprots{j}{k})) == 0
                    model = addProtein(model,uniprots{j}{k},kegg,swissprot);
                end
                newMets{k} = ['prot_' newMets{k}];
            end
            %Add new protein stoich. coeffs to rxn:
            kvalues = kcats(j)./stoich{j};
            rxnID   = model.rxns{i};
            rxnName = model.rxnNames{i};
            model   = addEnzymesToRxn(model,kvalues,rxnID,newMets,{rxnID,rxnName});
        end
    end
    %Update int_pos:
    S        = full(model.S);
    subs_pos = find(S(:,i) < 0);
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos)';
    
    %Individual Changes:
    for j = 1:length(int_pos)
        %INITIAL MANUAL CURATION:
        % Aconitase (P19414/EC 4.2.1.3): The rxn is represented as a two step rxn in Yeast
        % 7.5, so the kcat must be multiplied by 2 (2015-11-05)
        if ~isempty(strfind(model.rxns{i},'r_0280')) || ...  %mitochondrion
           ~isempty(strfind(model.rxns{i},'r_0302')) || ...  %mitochondrion
           ~isempty(strfind(model.rxns{i},'r_0303')) || ...  %cytoplasm
           ~isempty(strfind(model.rxns{i},'r_2305'))         %cytoplasm
            model.S(int_pos(j),i) = model.S(int_pos(j),i)/2;
        end
        %MANUAL CURATION FOR TOP USED ENZYMES:
        % FAS (P07149+P19097/EC2.3.1.86): No kcat available in BRENDA (it was using ECs of
        % subunits: EC3.1.2.14 and EC2.3.1.41).Value corrected with max. s.a. in S.cerevisiae
        % for NADPH in BRENDA (2015-08-25)
        if ~isempty(strfind(model.rxns{i},'r_2140'))
            model.S(int_pos(j),i) = -(3/14*60*1e3/1e3*MW_set)^-1;    %3/14 [umol/min/mg]
        elseif ~isempty(strfind(model.rxns{i},'r_2141'))
            model.S(int_pos(j),i) = -(3/16*60*1e3/1e3*MW_set)^-1;    %3/16 [umol/min/mg]
        end
        % Glycogen Synthase (P27472/EC2.4.1.11): No kcat available in BRENDA (it was using
        % the kcat from the subunit glycogenin glucosyltransferase (P47011/EC2.4.1.186) in
        % O.cuniculus). Value corrected with max. s.a. in S.cerevisiae (2015-08-27)
        if ~isempty(strfind(model.rxns{i},'r_0510'))
            model.S(int_pos(j),i) = -(90.5*60*1e3/1e3*MW_set)^-1;    %90.5 [umol/min/mg]
        end
        % Ketol-acid Reductoisomerase (P06168/EC1.1.1.86) - 2-acetyllactic acid: Substrate
        % name in BRENDA was a synonim as name in model. Changed manually (2015-08-28).
        if ~isempty(strfind(model.rxns{i},'r_0096'))
            model.S(int_pos(j),i) = -(18.3*3600)^-1;    %BRENDA: 2-acetolactate
        end
        % Ketol-acid Reductoisomerase (P06168/EC1.1.1.86) - (S)-2-acetyl-2-hydroxybutanoate:
        % Substrate name in BRENDA was a synonim as name in model. Changed manually
        % (2015-08-28).
        if ~isempty(strfind(model.rxns{i},'r_0669'))
            model.S(int_pos(j),i) = -(78.3*3600)^-1;    %BRENDA: 2-aceto-2-hydroxybutyrate
        end
        % Phosphoribosylformylglycinamidine synthase (P38972/EC6.3.5.3): Only kcat available
        % in BRENDA was for NH4 in E.coli. Value corrected with max. s.a. in E.coli
        % (2015-08-28).
        if ~isempty(strfind(model.rxns{i},'r_0079'))
            model.S(int_pos(j),i) = -(2*60*1e3/1e3*MW_set)^-1;   %2 [umol/min/mg]
        end
        % HMG-CoA reductase (P12683-P12684/EC1.1.1.34): Only kcat available in BRENDA was
        % for Rattus Norvegicus. Value corrected with max. s.a. in S.cerevisiae, checked in
        % original source (I.F.Durr & H.Rudney, J. Biol. Chem. 1960 235:2572-2578
        % http://www.jbc.org/content/235/9/2572) (2015-08-31).
        if ~isempty(strfind(model.rxns{i},'r_0558'))
            model.S(int_pos(j),i) = -(13.5*0.002*60*1e3/1e3*MW_set)^-1;     %13.5 units*0.002 umol = 0.027 [umol/min/mg]
        end
        % FPP synthase (P08524/EC2.5.1.1): The recommended EC is 2.5.1.10. Value corrected
        % with s.a. in S.cerevisiae from BRENDA (2015-08-31).
        if ~isempty(strfind(model.rxns{i},'r_0355')) || ...
           ~isempty(strfind(model.rxns{i},'r_0462'))
            model.S(int_pos(j),i) = -(2.33*60*1e3/1e3*MW_set)^-1;    %2.33 [umol/min/mg]
        end        
        % Amino-acid N-acetyltransferase (P40360-Q04728/EC2.3.1.1): Only kcat available
        % in BRENDA was for Mycobacterium tuberculosis. Value corrected with s.a. in
        % E.coli from BRENDA (2015-08-31).
        if ~isempty(strfind(model.rxns{i},'r_0761'))
            model.S(int_pos(j),i) = -(133*60*1e3/1e3*MW_set)^-1;    %133 [umol/min/mg]
        end        
        % Glutamate N-acetyltransferase (Q04728/EC2.3.1.1): Missanotated EC (should be
        % 2.3.1.35), and no kcats for the correct EC. Value corrected with s.a. in
        % S.cerevisiae from BRENDA (2015-08-31).
        if ~isempty(strfind(model.rxns{i},'r_0818'))
            model.S(int_pos(j),i) = -(22*60*1e3/1e3*MW_set)^-1;    %22 [umol/min/mg]
        end
        %MANUAL CURATION FOR CARBON SOURCES
        % alpha,alpha-trehalase (P32356/EC3.2.1.28): The available kcat values were not for
        % S.cerevisiae. Value corrected with s.a. in S.cerevisiae, checked in original source
        % (N.Biswas &  A.K.Ghosh, BBA. 1996 1290(1):95-100) (2015-09-02)
        if ~isempty(strfind(model.rxns{i},'r_0194'))
            model.S(int_pos(j),i) = -(22*60*1e3/1e3*MW_set)^-1;    %22 [umol/min/mg]
        end
        %MANUAL CURATION FOR STRESS DATA
        % 1. Ribose-phosphate pyrophosphokinase (Q12265/EC2.7.6.1): Substrate name
        % (ribose-5-phosphate)in BRENDA was a synonim as name in model so it was using ATP
        % as substrate instead (100 times lower). Changed manually (2015-10-05).
        if ~isempty(strfind(model.rxns{i},'r_0916'))
            model.S(int_pos(j),i) = -(60.68*3600)^-1;    %BRENDA: D-ribose 5-phosphate
        end
        % 2. Glutamine synthetase (P32288/EC6.3.1.2): No data in BRENDA for yeast or fungi.
        % S.A. retrieved from manual search (Mitchell & Magasanik, Journal of Bio.Chem. 
        % 1983, 258:119-124) (2015-10-05).
        if ~isempty(strfind(model.rxns{i},'r_0476'))
            model.S(int_pos(j),i) = -(236*60*1e3/1e3*MW_set)^-1;    %236 [umol/min/mg]
        end
        % 3. Chorismate synthase (P28777/EC4.2.3.5): Only available kcat was for N. crassa.
        % Replaced with s.a. of E.coli, from BRENDA (2015-10-05).
        if ~isempty(strfind(model.rxns{i},'r_0279'))
            model.S(int_pos(j),i) = -(14.8*60*1e3/1e3*MW_set)^-1;    %14.8 [umol/min/mg]
        end
        % 4. Homoaconitase, mitochondrial (P49367/EC4.2.1.36): Only available kcats are for
        % Methanocaldococcus jannaschii. Instead, we used a study in yeast (Strassman &
        % Ceci, Journal of Bio.Chem. 1966, 241:5401-5407) that shows that aconitase is
        % 0.062/0.005 = 12.4 times faster than homo-aconitase. Aconitase's kcat (143.3 1/s)
        % is taken from manual data and multiplied by 2 to represent both parts of the
        % reaction (2015-11-05).
        if ~isempty(strfind(model.rxns{i},'r_0027')) || ...
           ~isempty(strfind(model.rxns{i},'r_0542'))
            model.S(int_pos(j),i) = -(143.3*2/12.4*3600)^-1;    %Aconitase: 143.3 [1/s]
        end
        % 5. Formyltetrahydrofolate synthetase (P07245/EC6.3.4.3): kcat of S.cerevisiae in
        % BRENDA was hidden as 'additional information'. Fixed by going to the original
        % reference (Mejillano et al, Biochemistry 1989, 28:5136-5145) (2015-11-11).
        if ~isempty(strfind(model.rxns{i},'r_0446'))
            model.S(int_pos(j),i) = -(200*3600)^-1;             %200 [1/s]
        end
        % 5. Methenyltetrahydrofolate cyclohydrolase (P07245/EC6.3.4.3): Missanotated EC
        % (should be 3.5.4.9), value corrected with only kcat available (2015-11-11).
        if ~isempty(strfind(model.rxns{i},'r_0725'))
            model.S(int_pos(j),i) = -(134*3600)^-1;             %134 [1/s] - Homo sapiens
        end
        % 5. Methylenetetrahydrofolate dehydrogenase (P07245/EC6.3.4.3): Missanotated EC
        % (should be 1.5.1.5), and no kcats for the correct EC. Value corrected with s.a.
        % in S.cerevisiae from BRENDA (2015-11-11).
        if ~isempty(strfind(model.rxns{i},'r_0732'))
            model.S(int_pos(j),i) = -(259*60*1e3/1e3*MW_set)^-1;    %259 [umol/min/mg]
        end
        % 6. Phosphoserine transaminase (P33330/EC2.6.1.52): Only values were for E.coli
        % with fusion proteins. Value changed with s.a. found with manual search (Hirsch
        % & Greenberg, Journal of Bio.Chem. 1967, 242:2283-2287) (2015-11-11).
        if ~isempty(strfind(model.rxns{i},'r_0918'))
            model.S(int_pos(j),i) = -(78*60*1e3/1e3*MW_set)^-1;    %78 [umol/min/mg]
        end
        %ADITIONAL CURATION FOR CHEMOSTAT GROWTH
        % Succinate-semialdehyde dehydrogenase (P38067/EC1.2.1.16): Retrieved value was
        % from E.coli under extreme conditions. Value changed with s.a. in S.cerevisiae
        % from BRENDA (2017-02-13).
        if ~isempty(strfind(model.rxns{i},'r_1023'))
            model.S(int_pos(j),i) = -(0.66*60*1e3/1e3*MW_set)^-1;   %0.66 [umol/min/mg]
        end
        % 1,3-beta-glucan synthase component FKS1 (P38631/EC2.4.1.34): Retrieved value
        % was from Staphylococcus aureus. Value changed with s.a. in S.cerevisiae
        % from BRENDA (2017-03-05).
        if ~isempty(strfind(model.rxns{i},'r_0005'))
            model.S(int_pos(j),i) = -(4*60*1e3/1e3*MW_set)^-1;   %4 [umol/min/mg]
        end
    end
    disp(['Improving model with curated data: Ready with rxn #' num2str(i)])
end

% Other manual changes:

% Remove protein P14540 (Fructose-bisphosphate aldolase) from missanotated rxns
% (2017-01-16):
rxns_tochange = {'r_0322No1','r_0322_REVNo1','r_0990No1','r_0990_REVNo1'};
for i = 1:length(rxns_tochange)
    rxn_name = rxns_tochange{i};
    pos_rxn  = strcmp(model.rxns,rxn_name);
    model.S(strcmp(model.mets,'prot_P14540'),pos_rxn) = 0;
    model.rxns{pos_rxn} = rxn_name(1:end-3);
end

% Remove protein P40009 (Golgi apyrase) from missanotated rxns (2016-12-14):
model = removeRxns(model,'r_0227No3');

% Remove repeated reactions (2017-01-16):
rem_rxn = false(size(model.rxns));
for i = 1:length(model.rxns)-1
    for j = i+1:length(model.rxns)
        if isequal(model.S(:,i),model.S(:,j)) && model.lb(i) == model.lb(j) && ...
           model.ub(i) == model.ub(j) 
            rem_rxn(j) = true;
            disp(['Removing repeated rxn: ' model.rxns{i} ' & ' model.rxns{j}])
        end
    end
end
model = removeRxns(model,model.rxns(rem_rxn));

% Merge arm reactions to reactions with only one isozyme (2017-01-17):
arm_pos = zeros(size(model.rxns));
p       = 0;
for i = 1:length(model.rxns)
    rxn_name = model.rxns{i};
    if ~isempty(strfind(rxn_name,'arm_'))
        rxn_code  = rxn_name(5:end);
        k         = 0;
        for j = 1:length(model.rxns)
            if ~isempty(strfind(model.rxns{j},[rxn_code 'No']))
                k   = k + 1;
                pos = j;
            end
        end
        if k == 1
            %Condense both reactions in one:
            new_name   = model.rxns{pos};
            stoich     = model.S(:,i) + model.S(:,pos);
            model      = addReaction(model,new_name,model.mets,stoich,true,0,1000);
            p          = p + 1;
            arm_pos(p) = i;
            disp(['Merging reactions: ' model.rxns{i} ' & ' model.rxns{pos}])
        end
    end
end
% Remove saved arm reactions:
model = removeRxns(model,model.rxns(arm_pos(1:p)));

% Remove unused enzymes after manual curation (2017-01-16):
rem_enz = false(size(model.enzymes));
rem_met = false(size(model.mets));
for i = 1:length(model.enzymes)
    pos_met = strcmp(model.mets,['prot_' model.enzymes{i}]);
    if sum(model.S(pos_met,:)~=0) == 1
        rem_met(pos_met) = true;
        rem_enz(i)       = true;
        disp(['Removing unused protein: ' model.enzymes{i}])
    end
end
model.enzymes(rem_enz)   = [];
model.genes2(rem_enz)    = [];
model.geneNames(rem_enz) = [];
model.MWs(rem_enz)       = [];
model.sequence(rem_enz)  = [];
model.pathways(rem_enz)  = [];
model                    = removeMetabolites(model,model.mets(rem_met));

% Block O2 and glucose production (for avoiding multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;

% Remove incorrect pathways:
model = removeIncorrectPathways(model);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%