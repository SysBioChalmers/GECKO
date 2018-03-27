%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = manualModifications(model)
% 
% Benjamin J. Sanchez. Last edited: 2017-10-29
% Ivan Domenzain.      Last edited: 2018-01-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,modifications] = manualModificationsKma(model)

%Read manual data:
% cd ../../Databases
% fID           = fopen('manual_data.txt');
% data          = textscan(fID,'%s %s %s %s %f','delimiter','\t');
% structure     = data{2};
% kcats         = data{5}.*3600;
% data          = load('ProtDatabase.mat');
% swissprot     = data.swissprot;
% kegg          = data.kegg;
% fclose(fID);
% cd ../Matlab_Module/change_model
 modifications{1} = [];
 modifications{2} = [];
% 
% %Construct curated complexes:
% uniprots = cell(size(kcats));
% stoich   = cell(size(kcats));
% for i = 1:length(kcats)
%     uniprots{i} = strsplit(structure{i},' + ');
%     stoich{i}   = ones(size(uniprots{i}));
%     %Separate complex strings in units and amount of each unit:
%     for j = 1:length(uniprots{i})
%         unit = uniprots{i}{j};
%         pos  = strfind(unit,' ');
%         if isempty(pos)
%             stoich{i}(j)   = 1;
%             uniprots{i}{j} = unit;
%         else
%             stoich{i}(j)   = str2double(unit(1:pos-1));
%             uniprots{i}{j} = unit(pos+1:end);
%         end
%     end
% end

 for i = 1:length(model.rxns)
     reaction = model.rxnNames{i};
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
%     %Find intersection with manual curated data:
%     for j = 1:length(uniprots)
%         int    = intersect(prot_set,uniprots{j});
%         if length(int)/max(length(prot_set),length(uniprots{j})) > 0.50 % 50% match
%             %Erase previous protein stoich. coeffs from rxn:
%             for k = 1:length(prot_set)
%                 model.S(int_pos(k),i) = 0;
%             end
%             %If some proteins where not present previously, add them:
%             newMets = uniprots{j};
%             for k = 1:length(uniprots{j})
%                 if sum(strcmp(model.enzymes,uniprots{j}{k})) == 0
%                     model = addProtein(model,uniprots{j}{k},kegg,swissprot);
%                 end
%                 newMets{k}       = ['prot_' newMets{k}];
%             end
%             %Add new protein stoich. coeffs to rxn:
%             kvalues = kcats(j)./stoich{j};
%             rxnID   = model.rxns{i};
%             rxnName = model.rxnNames{i};
%             model   = addEnzymesToRxn(model,kvalues,rxnID,newMets,{rxnID,rxnName});
%         end
%     end
     %Update int_pos:
     S        = full(model.S);
     subs_pos = find(S(:,i) < 0);
     %Get the proteins that are part of the i-th rxn
     prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
     int_pos  = intersect(subs_pos,prot_pos)';
 %%%%%%%%%%%%%%%%%%%%%%%%%%%  Individual Changes:  %%%%%%%%%%%%%%%%%%%%%%%%
     for j = 1:length(int_pos)
         enzName = model.mets(int_pos(j));
         %%%%%%%%%%%%%%%%%% MANUAL CURATION FOR TOP GROWTH LIMITING ENZYMES:  
         [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications);
         if ~isempty(newValue)
             model.S(int_pos(j),i) = newValue;
         else
             %%%%%%%%%%%%%%%%%%%% MANUAL CURATION FOR CARBON SOURCES
             [newValue,modifications] = curation_topUsedEnz(reaction,enzName,MW_set,modifications);
             if ~isempty(newValue)
                 model.S(int_pos(j),i) = newValue;
             %else
                 %%%%%%%%%%%%%%% MANUAL CURATION FOR TOP USED ENZYMES:
%                  [newValue,modifications] = curation_carbonSources(reaction,enzName,MW_set,modifications);
%                  if ~isempty(newValue)
%                      model.S(int_pos(j),i) = newValue;
%                  end
             end
          end          
      end
     disp(['Improving model with curated data: Ready with rxn #' num2str(i)])
end
%%%%%%%%%%%%%%%%%%%%%%%%% Other manual changes: %%%%%%%%%%%%%%%%%%%%%%%%%%%
%model = otherChanges(model);
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
%Change gene rules:
for i = 1:length(model.rules)
    if ~isempty(model.rules{i})
        %Change gene ids:
        model.rules{i} = strrep(model.rules{i},'x(','');
        model.rules{i} = strrep(model.rules{i},')','');
        model.rules{i} = model.genes{str2double(model.rules{i})};
    end
end
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
%model         = removeIncorrectPathways(model);
% Map the index of the modified Kcat values to the new model (after rxns
% removals).
modifications = mapModifiedRxns(modifications,model);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the top growth limiting enzymes that were detected by the
% modifyKcats.m script in a preliminary run. 
function [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications)  
        newValue = [];
        reaction = string(reaction);
        
        % 3-hydroxy-3-methylglutaryl coenzyme A reductase (W0TEE1/EC1.1.1.34): 
        % Only kcat available in BRENDA was for Rattus Norvegicus. Value 
        % corrected with max. s.a. in Rattus norvegicus [0.03 umol/min/mg, Mw=226 kDa]
        % from BRENDA (2018-01-27)
          if (strcmpi('prot_W0TEE1',enzName)) && (~isempty(strfind(reaction,'(R)-Mevalonate:NADP+ oxidoreductase (CoA acylating)')))
             newValue         = -(0.03*226000*0.06)^-1;
            % Homarus americanus
             %newValue         = -(60*540000*0.06)^-1;                        
             modifications{1} = [modifications{1}; string('W0TEE1')];
             modifications{2} = [modifications{2}; reaction];
          end
          
        % ATP:adenosine 5'-phosphotransferase (W0TF72/EC2.7.1.20): 
        % No Kcat available for K. marxianus, a value for the same
        % substrate and for S. cerevisiae was found in BRENDA -> 25.5 1/s
        % [2018-03-23]
          if (strcmpi('prot_W0TF72',enzName)) && (~isempty(strfind(reaction,'ATP:adenosine 5')))
             newValue         = -(234*3600)^-1;
             modifications{1} = [modifications{1}; string('W0TF72')];
             modifications{2} = [modifications{2}; reaction];
          end
 
         % [1.2.1.12] Kcat (Sce & natural substrate) 29 [1/s] is a highly 
        % growthRate limiting value, S.A. 100 (Sce)
          %enzIDs = {'prot_P00359','prot_P00360','prot_P00358'};
          if ~isempty(strfind(reaction,'glyceraldehyde-3-phosphate'))
              if (strcmpi('prot_W0T9E0',enzName))
                  %newValue     = -(29*3600)^-1;
                 newValue = -(253*1e3/1e3*MW_set)^-1;
                 modifications{1} = [modifications{1}; string('W0T9E0')];
                 modifications{2} = [modifications{2}; reaction];
%               elseif strcmpi('prot_P00358',enzName)
%                   newValue      = -(16.7*3600)^-1;
%                   %newValue = -(100*1e3/1e3*MW_set)^-1;
%                   modifications{1} = [modifications{1}; string('P00358')];
%                   modifications{2} = [modifications{2}; reaction];
%               elseif strcmpi('prot_P00359',enzName)
%                   newValue      = -(9.1*3600)^-1;
%                   %newValue = -(100*1e3/1e3*MW_set)^-1;
%                   modifications{1} = [modifications{1}; string('P00359')];
%                   modifications{2} = [modifications{2}; reaction];
               end
          end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify those kcats involved in extreme misspredictions for growth on 
% several carbon sources. This values were obtained by specific searches on
% the involved pathways for the identification of the ec numbers and then
% its associated Kcat values were gotten from BRENDA.
function [newValue,modifications] = curation_carbonSources(reaction,enzName,MW_set,modifications)
        newValue = [];
        reaction = string(reaction);
       % alpha,alpha-trehalase (P32356/EC3.2.1.28): The growth on
       % alpha,alpha-trehalose was extremely highly overpredicted, the
       % minimum substrate specific kcat value (0.67 1/s) is chosen instead
       % (Metarhizium flavoviride) 
       % from BRENDA (2018-01-22).
         if (strcmpi('prot_P32356',enzName)  || strcmpi('prot_P48016',enzName)) && ...
                      (~isempty(strfind(reaction,'alpha,alpha-trehalase')))
             newValue      = -(0.67*60*1e3/1e3*MW_set)^-1; 
             modifications{1} = [modifications{1}; string('P48016')];
             modifications{2} = [modifications{2}; reaction];
         end
        % alcohol dehydrogenase (ethanol to acetaldehyde)[P00331/EC1.1.1.1] 
        % Growth on ethanol was extremely highly overpredicted, the
        % specific Kcat value for ethanol found in BRENDA is 143 (1/s) @pH
        % 8.0 and 25 C, and 432.3 @20?C and pH 9.0 
       % from BRENDA (2018-01-22).
         if strcmpi('prot_P00331',enzName) && (~isempty(strfind(reaction,...
                      'alcohol dehydrogenase (ethanol to acetaldehyde)')))
             newValue      = -(143*3600)^-1;
             modifications{1} = [modifications{1}; string('P00331')];
             modifications{2} = [modifications{2}; reaction];
         end
        % Glycerol dehydrogenase [P14065/EC1.1.1.1] 
        % Growth on glycerol was highly overpredicted, the minimum S.A. 
        % value for  glycerol found in BRENDA is 0.54 (umol/ming/mg) and
        % the maximal is 15.7, both for Schizosaccharomyces pombe. The MW
        % is 57 kDa. From BRENDA (2018-01-26).
         if strcmpi('prot_P14065',enzName) && (~isempty(strfind(reaction,...
                      'glycerol dehydrogenase (NADP-dependent)')))
             newValue      = -(0.54*57000*0.06)^-1; 
             modifications{1} = [modifications{1}; string('P14065')];
             modifications{2} = [modifications{2}; reaction];
         end
        % alpha-glucosidase [EC3.2.1.10/EC3.2.1.20] 
        % Growth on maltose is still underpredicted, for this reason all
        % the isoenzymes catalysing its conversion are set up to the
        % maximal glucosidase rate for maltose (3.2.1.20) <- 709 1/s for 
        % Schizosaccharomyces pombe, a value of 1278 is also available for 
        % Sulfolobus acidocaldarius.
        % From BRENDA (2018-01-28).
         if ~isempty(strfind(reaction,'alpha-glucosidase'))
              newValue = -(1278*3600)^-1; 
         end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After the growth limiting Kcats analysis and the curation for several
% carbon sources, a simulation for the model growing on minimal glucose
% media yielded a list of the top used enzymes (mass-wise), those that were
% taking more than 10% of the total proteome are chosen for manual curation
function [newValue,modifications] = curation_topUsedEnz(reaction,enzName,MW_set,modifications)
    newValue = [];
    reaction = string(reaction);
    % S-Adenosyl-L-methionine:zymosterol C-methyltransferase [W0TAZ1/EC2.1.1.41]
    % No Kcat reported on BRENDA, the maximum S.A. (E. coli is taken
    % instead)
    if strcmpi('prot_W0TAZ1',enzName) && ...
            (~isempty(strfind(reaction,'zymosterol C-methyltransferase')))
        newValue      = -(0.53*172000*60/1000)^-1;
        modifications{1} = [modifications{1}; string('W0TAZ1')];
        modifications{2} = [modifications{2}; reaction];
    end
    % Ketol-acid Reductoisomerase (W0TAV6/EC1.1.1.86)
    % Substrate name in BRENDA was a synonim as name in model. Changed manually
    % (2018-03-24).
    if strcmpi('prot_W0TAV6',enzName) && (~isempty(strfind(reaction,'zymosterol C-methyltransferase')))
        newValue         = -(78.3*3600)^-1;
        modifications{1} = [modifications{1}; string('W0TAV6')];
        modifications{2} = [modifications{2}; reaction];
    end
    % Amidophosphoribosyltransferase (W0T3E2/EC2.4.2.14)
    % No reported values for K. marxianus, E. coli was the maximal S.A.
    % value reported for a microbial organism.
    % (2018-03-24).
    if strcmpi('prot_W0T3E2',enzName) && (~isempty(strfind(reaction,'5-phosphoribosylamine:')))
        newValue         = -(17.2*194000*60/1000)^-1;
        modifications{1} = [modifications{1}; string('W0T3E2')];
        modifications{2} = [modifications{2}; reaction];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = otherChanges(model)
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
    % Remove protein P40009 (Golgi apyrase) from missanotated rxns (2016-12-14):
    % model = removeRxns(model,{'r_0227No3'});
      pos_rxn = find(~cellfun(@isempty,strfind(model.rxnNames,'ATPase, cytosolic (No3)')));
      if ~isempty(pos_rxn) && pos_rxn~=0
         model = removeRxns(model,model.rxns(pos_rxn));
      end
    % Remove 2 proteins from missanotated rxns: Q12122 from cytosolic rxn (it's
    % only mitochondrial) & P48570 from mitochondrial rxn (it's only cytosolic).
    % Also rename r_0543No2 to r_0543No1 (for consistency) (2017-08-28):
    model = removeRxns(model,'r_0543No1');
    model = removeRxns(model,'r_1838No2');
    index = find(strcmp(model.rxns,'r_0543No2'));
    if isempty(index)
        model.rxnNames{index} = 'homocitrate synthase (No1)';
        model.rxns{index}     = 'r_0543No1';
    end
    % Aconitase (P19414/EC 4.2.1.3): The rxn is represented as a two
    % step rxn in Yeast7.5, so the kcat must be multiplied by 2
    %(2015-11-05)
    index = find(strcmpi('prot_P19414',model.mets));
    if ~isempty(index)
        rxnIndxs = find(model.S(:,index));
        if ~isempty(rxnIndxs)
            rxnIndxs = rxnIndxs(1:end-2);
            model.S(rxnIndxs,index) = model.S(rxnIndxs,index)/2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modified = mapModifiedRxns(modifications,model)
    modified = [];
    for i=1:length(modifications{1})
        rxnIndex = find(strcmp(model.rxnNames,modifications{2}(i)),1);
        str      = {horzcat(modifications{1}{i},'_',num2str(rxnIndex))};
        modified = [modified; str];
    end
end
