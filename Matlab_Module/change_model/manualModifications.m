%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = manualModifications(model)
% 
%
% Benjam?n J. S?nchez. Last edited: 2017-10-29
% Ivan Domenzain.      Last edited: 2018-01-25
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
    %Get the proteins that are part of the i-th rxn
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos)';
    
    %Individual Changes:
    for j = 1:length(int_pos)
%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL MANUAL CURATION: %%%%%%%%%%%%%%%%%%%%%%%%
         % Aconitase (P19414/EC 4.2.1.3): The rxn is represented as a two step rxn in Yeast
         % 7.5, so the kcat must be multiplied by 2 (2015-11-05)
         if strcmpi('prot_P19414',model.mets(int_pos(j))) && ...
            (~isempty(strfind(model.rxnNames{i},'cis-aconitate(3-) to isocitrate')) || ...
             ~isempty(strfind(model.rxnNames{i},'citrate to cis-aconitate(3-) (')) || ...
             ~isempty(strfind(model.rxnNames{i},'citrate to cis-aconitate(3-) (')) || ...
             ~isempty(strfind(model.rxnNames{i},'citrate to cis-aconitate(3-), cytoplasmic')))             
             model.S(int_pos(j),i) = model.S(int_pos(j),i)/2;
         end
%%%%%%%%%%%%%%%%%%%%%%%%%% MANUAL CURATION FOR TOP GROWTH LIMITING ENZYMES:       
           % [1.2.1.12] Kcat (Sce & natural substrate) 29 [1/s], S.A. 100 (Sce)
%          if (strcmpi('prot_P00358',model.mets(int_pos(j))) || ...
%              strcmpi('prot_P00359',model.mets(int_pos(j))) || ...
%              strcmpi('prot_P00360',model.mets(int_pos(j)))) && ...
          if strcmpi('prot_P00359',model.mets(int_pos(j))) &&...
             ~isempty(strfind(model.rxnNames{i},...
                               'glyceraldehyde-3-phosphate dehydrogenase'))
              %model.S(int_pos(j),i) = -(29*3600)^-1;
              %model.S(int_pos(j),i) = -(9.1*3600)^-1;
              model.S(int_pos(j),i) = -(100*1e3/1e3*MW_set)^-1;
          end
        % FAS (P07149+P19097/EC2.3.1.86): No kcat available in BRENDA (it was using ECs of
        % subunits: EC3.1.2.14 and EC2.3.1.41).Value corrected with max. s.a. in S.cerevisiae
        % for NADPH in BRENDA (2015-08-25)
          if (strcmpi('prot_P07149',model.mets(int_pos(j)))  || ...
              strcmpi('prot_P19097',model.mets(int_pos(j))))
            if ~isempty(strfind(model.rxnNames{i},...
                                   'fatty-acyl-CoA synthase (n-C16:0CoA)')) 
                model.S(int_pos(j),i) = -(3/14*60*1e3/1e3*MW_set)^-1;    %3/14 [umol/min/mg]
            elseif ~isempty(strfind(model.rxnNames{i},...
                                  'fatty-acyl-CoA synthase (n-C18:0CoA)'))
                model.S(int_pos(j),i) = -(3/16*60*1e3/1e3*MW_set)^-1;    %3/16 [umol/min/mg]
            end
          end
        % Ketol-acid Reductoisomerase (P06168/EC1.1.1.86) - 2-acetyllactic acid: Substrate
        % name in BRENDA was a synonim as name in model. Changed manually (2015-08-28).
          if strcmpi('prot_P06168',model.mets(int_pos(j))) 
              if (~isempty(strfind(model.rxnNames{i},...
                                  'acetohydroxy acid isomeroreductase (')))
                 model.S(int_pos(j),i) = -(18.3*3600)^-1;    %BRENDA: 2-acetolactate
              elseif (~isempty(strfind(model.rxnNames{i},...
                  'ketol-acid reductoisomerase (2-aceto-2-hydroxybutanoate) (')))
                 model.S(int_pos(j),i) = -(78.3*3600)^-1;    %BRENDA: 2-aceto-2-hydroxybutyrate
              end     
          end
        % Phosphoribosylformylglycinamidine synthase (P38972/EC6.3.5.3): Only kcat available
        % in BRENDA was for NH4 in E.coli. Value corrected with max. s.a. in E.coli
        % (2015-08-28).
          if strcmpi('prot_P38972',model.mets(int_pos(j))) && ...
             (~isempty(strfind(model.rxnNames{i},'phosphoribosylformyl glycinamidine synthetase')))
             model.S(int_pos(j),i) = -(2.15*141418*60/1000)^-1;    
          end
        % HMG-CoA reductase (P12683-P12684/EC1.1.1.34): Only kcat available in BRENDA was
        % for Rattus Norvegicus. Value corrected with max. s.a. in S.cerevisiae, checked in
        % original source (I.F.Durr & H.Rudney, J. Biol. Chem. 1960 235:2572-2578
        % http://www.jbc.org/content/235/9/2572) (2015-08-31).
          if (strcmpi('prot_P12683',model.mets(int_pos(j)))  || ...
              strcmpi('prot_P12684',model.mets(int_pos(j)))) && ...
              (~isempty(strfind(model.rxnNames{i},'hydroxymethylglutaryl')))
             model.S(int_pos(j),i) = -(13.5*0.001*60*1e3/1e3*MW_set)^-1;
             %model.S(int_pos(j),i) = -(3.0851*3600)^-1;
          end
          % amidophosphoribosyltransferase [P04046/EC2.4.2.14]
          % No Kcat reported on BRENDA, the maximum S.A. (E. coli is taken
          % instead)
          if strcmpi('prot_P04046',model.mets(int_pos(j))) && ...
             (~isempty(strfind(model.rxnNames{i},'phosphoribosylpyrophosphate amidotransferase')))
             model.S(int_pos(j),i) = -(17.2*194000*60/1000)^-1;
          end
         % Enolase (1&2) [4.2.1.11]
         % 71.4 (1/s) is the Kcat reported for 2-phospho-D-glycerate
         % 230 (1/s) is the Kcat reported for 2-phosphoglycerate
         % both measurements are for native enzymes
        %if (strcmpi('prot_P00924',model.mets(int_pos(j))) || ...
        %    strcmpi('prot_P00925',model.mets(int_pos(j)))) && ...
           if strcmpi('prot_P00924',model.mets(int_pos(j)))  && ...
              (~isempty(strfind(model.rxnNames{i},'enolase'))) 
                model.S(int_pos(j),i) = -1/(230*3600);
           end
           
        % [4.1.1.-, 4.1.1.43, 4.1.1.72, 4.1.1.74] Pyruvate decarboxylase 
        %  Resulted to be a growth limiting enzyme but the Kcat
        % value seems to be the best candidate for this reaction
        %if (strcmpi('prot_P06169',model.mets(int_pos(j))) || ...
        %     strcmpi('prot_P26263',model.mets(int_pos(j)))) && ...
         if strcmpi('prot_P06169',model.mets(int_pos(j))) && ...
            (~isempty(strfind(model.rxnNames{i},'pyruvate decarboxylase')))
            model.S(int_pos(j),i) = -1/(145*3600); 
         end
         
         % 1,3-beta-glucan synthase component FKS1 (P38631/EC2.4.1.34): Retrieved value
         % was from Staphylococcus aureus. Value changed with s.a. in S.cerevisiae
         % from BRENDA (2017-03-05).
          if strcmpi('prot_P38631',model.mets(int_pos(j))) && ...
             (~isempty(strfind(model.rxnNames{i},'1,3-beta-glucan synthase (')))
             model.S(int_pos(j),i) = -(4*60*1e3/1e3*MW_set)^-1;    %4 [umol/min/mg]
             disp(['P38631' ' modified'])
             
          end
        % [Q06817//EC6.1.1.14] glycyl-tRNA synthetase
        %  Resulted to be a growth limiting enzyme, maximum catalytic value
        % found by the automatic algorithm was 15.9 1/s
        % from BRENDA (2018-01-25).
         if strcmpi('prot_Q06817',model.mets(int_pos(j))) && ...
            (~isempty(strfind(model.rxnNames{i},'glycyl-tRNA synthetase')))
            model.S(int_pos(j),i) = -1/(15.9*3600); 
         end
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MANUAL CURATION FOR CARBON SOURCES
        % alpha,alpha-trehalase (P32356/EC3.2.1.28): The available kcat values were not for
        % S.cerevisiae. Value corrected with s.a. in S.cerevisiae, checked in original source
        % (N.Biswas &  A.K.Ghosh, BBA. 1996 1290(1):95-100) (2015-09-02)
        if ~isempty(strfind(model.rxns{i},'r_0194'))
            model.S(int_pos(j),i) = -(44*60*1e3/1e3*MW_set)^-1;    %44 [umol/min/mg]
        end
          if (strcmpi('prot_P32356',model.mets(int_pos(j)))  || ...
              strcmpi('prot_P48016',model.mets(int_pos(j)))) && ...
             (~isempty(strfind(model.rxnNames{i},'alpha,alpha-trehalase')))
             model.S(int_pos(j),i) = -(0.67*60*1e3/1e3*MW_set)^-1;    
          end        
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MANUAL CURATION FOR TOP USED ENZYMES:
          % atp synthase mitochondrial (P07251-P00830/EC3.6.3.14): No Kcat 
          % reported for S. cerevisiae, value for Geobacillus stearothermophilus
          % for ATP is used instead 217 [1/s]-
          % from BRENDA (2017-01-18).
          if (strcmpi('prot_P07251',model.mets(int_pos(j))) || ...
              strcmpi('prot_P00830',model.mets(int_pos(j))))&& ...
             (~isempty(strfind(model.rxnNames{i},'ATP synthase (No1)')))
             model.S(int_pos(j),i) = -(217*3600)^-1; 
          end
          % glyceraldehyde-3-phosphate dehydrogenase 1 (gapdh 1) 
          % (P00360/EC1.2.1.12): The maximum catalytic value available for
          % Sce (S.A. = 100 umol/min/mg) was taking ~10-12% of the
          % simulated proteome, so the highest Kcat value available for the
          % natural substrate is chosen instead (for H. sapiens
          % from BRENDA (2017-01-18).
          if strcmpi('prot_P00360',model.mets(int_pos(j))) &&...
             ~isempty(strfind(model.rxnNames{i},...
                               'glyceraldehyde-3-phosphate dehydrogenase'))
              model.S(int_pos(j),i) = -(199*3600)^-1;
          end
          % transaldolase (reversible) (P15019/EC2.2.1.2): 
          % The protein usage is represents around 10% of the used proteome
          % on several carbon sources (batch simulations). The highest S.A.
          % was used instead (60*75000*60/1000) [1/hr] for E. coli
          % from BRENDA (2017-01-18).
          if strcmpi('prot_P15019',model.mets(int_pos(j))) &&...
             ~isempty(strfind(model.rxnNames{i},...
                               'transaldolase (reversible)'))
              model.S(int_pos(j),i) = -(60*75000*60/1000)^-1;
          end
           % Glutamate N-acetyltransferase (Q04728/EC2.3.1.1): Missanotated 
           % EC (should be 2.3.1.35), and no kcats for the correct EC. 
           % Value corrected with s.a. in S.cerevisiae from BRENDA (2015-08-31).
          if strcmpi('prot_Q04728',model.mets(int_pos(j))) && ...
             (~isempty(strfind(model.rxnNames{i},'ornithine transacetylase (No1)')))
             model.S(int_pos(j),i) = -(22*60*1e3/1e3*MW_set)^-1;    %22 [umol/min/mg]
          end
           % P80235/EC 2.3.1.7 - carnitine O-acetyltransferase 
           % The assigned Kcat was for Mus musculus and resulted to take
           % 20% of the used proteome on batch simulation for several
           % carbon sources. Sce S.A. is used instead [200 umol/min/mg]
           % from BRENDA (2018-01-18)
          if strcmpi('prot_P802358',model.mets(int_pos(j))) && ...
             (~isempty(strfind(model.rxnNames{i},'carnitine O-acetyltransferase')))
             model.S(int_pos(j),i) = -(200*60*1e3/1e3*MW_set)^-1;    %22 [umol/min/mg]
          end
          
     end
    disp(['Improving model with curated data: Ready with rxn #' num2str(i)])
end

%%%%%%%%%%%%%%%%%%%%%%%%% Other manual changes: %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove protein P14540 (Fructose-bisphosphate aldolase) from missanotated rxns
% (2017-01-16):
 rxns_tochange = {'D-fructose 1-phosphate D-glyceraldehyde-3-phosphate-lyase (No1)';...
                  'sedoheptulose 1,7-bisphosphate D-glyceraldehyde-3-phosphate-lyase (No1)';...
                  'D-fructose 1-phosphate D-glyceraldehyde-3-phosphate-lyase (reversible) (No1)';...
                  'sedoheptulose 1,7-bisphosphate D-glyceraldehyde-3-phosphate-lyase (reversible) (No1)'};
 for i = 1:length(rxns_tochange)
     rxn_name = rxns_tochange{i};
     pos_rxn  = find(strcmpi(model.rxnNames,rxn_name));
     rxn_id   = model.rxns{pos_rxn};
     model.S(strcmpi(model.mets,'prot_P14540'),pos_rxn) = 0;
     model.rxns{pos_rxn} = rxn_id(1:end-3);
 end
% 
% % Remove protein P40009 (Golgi apyrase) from missanotated rxns (2016-12-14):
 %disp(class('r_0227No3'))
% model = removeRxns(model,{'r_0227No3'});
  pos_rxn = find(~cellfun(@isempty,strfind(model.rxnNames,'ATPase, cytosolic (No3)')));
  if ~isempty(pos_rxn) && pos_rxn~=0
     model = removeRxns(model,model.rxns(pos_rxn));
  end

% % Remove 2 proteins from missanotated rxns: Q12122 from cytosolic rxn (it's
% % only mitochondrial) & P48570 from mitochondrial rxn (it's only cytosolic).
% % Also rename r_0543No2 to r_0543No1 (for consistency) (2017-08-28):
% model = removeRxns(model,'r_0543No1');
% model = removeRxns(model,'r_1838No2');
% model.rxnNames{strcmp(model.rxns,'r_0543No2')} = 'homocitrate synthase (No1)';
% model.rxns{strcmp(model.rxns,'r_0543No2')}     = 'r_0543No1';

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
    rxn_id = model.rxns{i};
    if ~isempty(strfind(rxn_id,'arm_'))
        rxn_code  = rxn_id(5:end);
        k         = 0;
        for j = 1:length(model.rxns)
            if ~isempty(strfind(model.rxns{j},[rxn_code 'No']))
                k   = k + 1;
                pos = j;
            end
        end
        if k == 1
            %Condense both reactions in one:
            new_id     = model.rxns{pos};
            new_name   = model.rxnNames{pos};
            stoich     = model.S(:,i) + model.S(:,pos);
            model      = addReaction(model,{new_id,new_name},model.mets,stoich,true,0,1000);
            p          = p + 1;
            arm_pos(p) = i;
            disp(['Merging reactions: ' model.rxns{i} ' & ' model.rxns{pos}])
        end
    end
end
% Remove saved arm reactions:
model = removeRxns(model,model.rxns(arm_pos(1:p)));


% Possible incompatibillity with other MATLAB versions (needs to be
% corrected).
% %Change gene rules:
% for i = 1:length(model.rules)
%     if ~isempty(model.rules{i})
%         %Change gene ids:
%         model.rules{i} = strrep(model.rules{i},'x(','');
%         model.rules{i} = strrep(model.rules{i},')','');
%         model.rules{i} = model.genes{str2double(model.rules{i})};
%     end
% end

% Remove unused enzymes after manual curation (2017-01-16):
rem_enz = false(size(model.enzymes));
for i = 1:length(model.enzymes)
    pos_met = strcmp(model.mets,['prot_' model.enzymes{i}]);
    if sum(model.S(pos_met,:)~=0) == 1
        rem_enz(i) = true;
    end
end
rem_enz = model.enzymes(rem_enz);
for i = 1:length(rem_enz)
    model = deleteProtein(model,rem_enz{i});
    disp(['Removing unused protein: ' rem_enz{i}])
end

% Block O2 and glucose production (for avoiding multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;

% Remove incorrect pathways:
model = removeIncorrectPathways(model);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%