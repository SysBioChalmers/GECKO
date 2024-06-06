function kcatList = fuzzyKcatMatching(model, ecRxns, modelAdapter, forceWClvl)
% fuzzyKcatMatching
%   Matchs the model EC numbers and substrates to the BRENDA database, to
%   return the corresponding kcats for each reaction. If no exact match is
%   found, less specific kcat values are found from (a) evolutionary
%   closely related organism; (b) different substrate; (c) calculated from
%   specific activities; (d) wildcards in the EC number. The model organism
%   is provided in the model adapter as obj.params.org_name, and
%   evolutionary distance to other organisms is determined via KEGG
%   phylogeny. If an organism name occurs multiple times in KEGG, the first
%   instance will be used when determining evolutionary distance.
%
% Input:
%   model        an ecModel in GECKO 3 format (with ecModel.ec structure)
%   ecRxns       for which reactions (from model.ec.rxns) kcat values should
%                be found, provided as logical vector with same length as
%                model.ec.rxns. (Opt, default is all reactions)
%   modelAdapter a loaded model adapter (Optional, will otherwise use the
%                default model adapter).
%   forceWClvl   force a minimum wildcard level (Optional, default 0). 
%
% Output:
%   kcatList    structure array with list of BRENDA derived kcat values,
%               with separate entries for each kcat value
%               source      'brenda'           
%               rxns        reaction identifiers
%               substrate   substrate names
%               kcat        proposed kcat value in /sec
%               eccodes     as used to query BRENDA
%               wildCardLvl which level of EC wild-card was necessary to
%                           find a match
%                           0: w.x.y.z
%                           1: w.x.y.-
%                           2: w.x.-.-
%                           3: w.-.-.-
%               origin      which level of specificity was necessary to
%                           find a match
%                           1: correct organism, correct substrate, kcat
%                           2: any organism, correct substrate, kcat
%                           3: correct organism, any substrate, kcat
%                           4: any organism, any substrate, kcat
%                           5: correct organism, specific activity
%                           6: any organism, specific activity
%
%   Note: If a wildcard is used, origin levels 1 and 2 are ignored. The
%   last digits in the E.C. number indicate the substrate specificity, so
%   if this should be ignored, then correct substrate matches should not be
%   prioritized.
%
% Usage:
%   kcatList = fuzzyKcatMatching(model, ecRxns, modelAdapter, forceWClvl)

if nargin<2 || isempty(ecRxns)
    ecRxns = true(numel(model.ec.rxns),1);
elseif isnumeric(ecRxns)
    ecRxnsVec = false(numel(model.ec.rxns),1);
    ecRxnsVec(ecRxns) = true;
    ecRxns = ecRxnsVec;
end
ecRxns=find(ecRxns); % Get indices instead of logical

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if nargin < 4 || isempty(forceWClvl)
    forceWClvl = 0;
end

if ~isfield(model.ec,'eccodes')
    error('No EC codes defined in model.ec.eccodes. First run getECfromGEM() and/or getECfromDatabase().')
end
eccodes      = model.ec.eccodes(ecRxns);
substrates   = cell(numel(ecRxns),1);
substrCoeffs = cell(numel(ecRxns),1);

%Need to remove the prefix of GECKO light rxn names in the ec structure
if ~model.ec.geckoLight
    rxnNames = model.ec.rxns;
else
    rxnNames = extractAfter(model.ec.rxns, 4);
end
[~,originalRxns] = ismember(rxnNames(ecRxns),model.rxns);
for i = 1:length(ecRxns)
    sel = find(model.S(:,originalRxns(i)) < 0);
    substrates{i}  = model.metNames(sel); 
    substrCoeffs{i} = -model.S(sel,originalRxns(i));
end

%Load BRENDA data:
[KCATcell, SAcell] = loadBRENDAdata(modelAdapter);

%Creates a Structure with KEGG codes for organisms, names and taxonomical
%distance matrix and extract the organism index in the KEGG struct
phylDistStruct =  KEGG_struct(modelAdapter.getPhylDistStructPath());
%Get the KEGG code for the model's organism
org_name       = params.org_name;
org_index      = find_inKEGG(org_name,phylDistStruct.names);
%build an index for genus in the phyl dist struct
%first just extract the genus (i.e. the first part of the name)
phylDistStruct.genus = lower(regexprep(phylDistStruct.names,'\s.*',''));
%create a map for the genuses
phylDistStruct.uniqueGenusList = unique(phylDistStruct.genus);
phylDistStruct.genusHashMap = containers.Map(phylDistStruct.uniqueGenusList,1:length(phylDistStruct.uniqueGenusList));
phylDistStruct.uniqueGenusIndices = cell(length(phylDistStruct.uniqueGenusList),1);

%Then for each genus create a list with indices to the names
for i = 1:length(phylDistStruct.genus)
    matchInd = cell2mat(values(phylDistStruct.genusHashMap,phylDistStruct.genus(i)));
    phylDistStruct.uniqueGenusIndices{matchInd} = [phylDistStruct.uniqueGenusIndices{matchInd};i];
end

%Allocate output
kcats = zeros(length(eccodes),1);
mM = length(eccodes);

%Create empty kcatInfo
%Legacy, no longer given as output, rather used to construct
%kcatList.wildcardLvl and kcatList.origin.
kcatInfo.info.org_s   = zeros(mM,1);
kcatInfo.info.rest_s  = zeros(mM,1);
kcatInfo.info.org_ns  = zeros(mM,1);
kcatInfo.info.rest_ns = zeros(mM,1);
kcatInfo.info.org_sa  = zeros(mM,1);
kcatInfo.info.rest_sa = zeros(mM,1);
kcatInfo.info.wcLevel = NaN(mM,1);
kcatInfo.stats.queries  = 0;
kcatInfo.stats.org_s    = 0;
kcatInfo.stats.rest_s   = 0;
kcatInfo.stats.org_ns   = 0;
kcatInfo.stats.rest_ns  = 0;
kcatInfo.stats.org_sa   = 0;
kcatInfo.stats.rest_sa  = 0;
kcatInfo.stats.wc0      = 0;
kcatInfo.stats.wc1      = 0;
kcatInfo.stats.wc2      = 0;
kcatInfo.stats.wc3      = 0;
kcatInfo.stats.wc4      = 0;
kcatInfo.stats.matrix   = zeros(6,5);

%build an EC index to speed things up a bit - many of the ECs appear
%many times - unnecessary to compare them all
%so, here, each EC string appears only once, and you get a vector with
%indices to the rows in KCATcell
[ECIndexIds,~,ic] = unique(KCATcell{1});
EcIndexIndices = cell(length(ECIndexIds),1);
for i = 1:length(EcIndexIndices)
    EcIndexIndices{i} = find(ic == i).';
end

%Apply force wildcard level
while forceWClvl > 0
    eccodes=regexprep(eccodes,'(.)*(\.\d+)(\.-)*$','$1\.-$3');
    forceWClvl = forceWClvl - 1;
end
if forceWClvl == 1
    eccodes = regexprep(eccodes,'.*','-\.-\.-\.-');
end

progressbar('Gathering kcat values by fuzzy matching to BRENDA database')
%Main loop:
for i = 1:mM
    %Match:
    EC = eccodes{i};
    if ~isempty(EC)
        EC = strsplit(EC,';');
        %Try to match direct reaction:
        if ~isempty(substrates{i})
            [kcats(i), kcatInfo.info,kcatInfo.stats] = iterativeMatch(EC,substrates{i},substrCoeffs{i},i,KCATcell,...
                kcatInfo.info,kcatInfo.stats,org_name,...
                phylDistStruct,org_index,SAcell,ECIndexIds,EcIndexIndices);
        end
    end
    progressbar(i/mM)
end

kcatList.source      = 'brenda';
kcatList.rxns        = model.ec.rxns(ecRxns);
kcatList.substrates  = substrates;
kcatList.kcats       = kcats;
kcatList.eccodes     = eccodes;
kcatList.wildcardLvl = kcatInfo.info.wcLevel;
kcatList.origin      = NaN(numel(model.ec.rxns(ecRxns)),1);
% This can be refactored, iterativeMatch and their nested functions can
% just directly report the origin number.
origin = [kcatInfo.info.org_s kcatInfo.info.rest_s kcatInfo.info.org_ns kcatInfo.info.rest_ns kcatInfo.info.org_sa kcatInfo.info.rest_sa];
for i=1:6
    kcatList.origin(find(origin(:,i))) = i;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kcat,dir,tot] =iterativeMatch(EC,subs,substrCoeff,i,KCATcell,dir,tot,...
    name,phylDist,org_index,SAcell,ECIndexIds,EcIndexIndices)
%Will iteratively try to match the EC number to some registry in BRENDA,
%using each time one additional wildcard.

kcat    = zeros(size(EC));
origin  = zeros(size(EC));
matches = zeros(size(EC));
wc_num  = ones(size(EC)).*1000;
for k = 1:length(EC)
    success  = false;
    while ~success
        %Atempt match:
        [kcat(k),origin(k),matches(k)] = mainMatch(EC{k},subs,substrCoeff,KCATcell,...
            name,phylDist,...
            org_index,SAcell,ECIndexIds,EcIndexIndices);
        %If any match found, ends. If not, introduces one extra wild card and
        %tries again:
        if origin(k) > 0
            success   = true;
            wc_num(k) = sum(EC{k}=='-');
        else
            dot_pos  = [2 strfind(EC{k},'.')];
            wild_num = sum(EC{k}=='-');
            wc_text  = '-.-.-.-';
            EC{k}    = [EC{k}(1:dot_pos(4-wild_num)) wc_text(1:2*wild_num+1)];
        end
    end
end

if sum(origin) > 0
    %For more than one EC: Choose the maximum value among the ones with the
    %less amount of wildcards and the better origin:
    best_pos   = (wc_num == min(wc_num));
    new_origin = origin(best_pos);
    best_pos   = (origin == min(new_origin(new_origin~=0)));
    max_pos    = find(kcat == max(kcat(best_pos)));
    wc_num     = wc_num(max_pos(1));
    origin     = origin(max_pos(1));
    matches    = matches(max_pos(1));
    kcat       = kcat(max_pos(1));

    %Update dir and tot:
    dir.org_s(i)   = matches*(origin == 1);
    dir.rest_s(i)  = matches*(origin == 2);
    dir.org_ns(i)  = matches*(origin == 3);
    dir.org_sa(i)  = matches*(origin == 4);
    dir.rest_ns(i) = matches*(origin == 5);
    dir.rest_sa(i) = matches*(origin == 6);
    dir.wcLevel(i) = wc_num;
    tot.org_s        = tot.org_s   + (origin == 1);
    tot.rest_s       = tot.rest_s  + (origin == 2);
    tot.org_ns       = tot.org_ns  + (origin == 3);
    tot.org_sa       = tot.org_sa  + (origin == 4);
    tot.rest_ns      = tot.rest_ns + (origin == 5);
    tot.rest_sa      = tot.rest_sa + (origin == 6);
    tot.wc0          = tot.wc0     + (wc_num == 0);
    tot.wc1          = tot.wc1     + (wc_num == 1);
    tot.wc2          = tot.wc2     + (wc_num == 2);
    tot.wc3          = tot.wc3     + (wc_num == 3);
    tot.wc4          = tot.wc4     + (wc_num == 4);
    tot.queries      = tot.queries + 1;
    tot.matrix(origin,wc_num+1) = tot.matrix(origin,wc_num+1) + 1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kcat,origin,matches] = mainMatch(EC,subs,substrCoeff,KCATcell,...
    name,phylDist,org_index,SAcell,ECIndexIds,EcIndexIndices)

%First make the string matching. This takes time, so we only want to do
%this once:
%Relaxes matching if wild cards are present:
wild     = false;
wild_pos = strfind(EC,'-');
if ~isempty(wild_pos)
    EC   = EC(1:wild_pos(1)-1);
    wild = true;
end
stringMatchesEC_cell = extract_string_matches(EC,KCATcell{1},wild,ECIndexIds,EcIndexIndices);

% Matching function prioritizing organism and substrate specificity when
% available.

origin = 0;
%First try to match organism and substrate:
[kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,name,true,false,...
    phylDist,org_index,SAcell,stringMatchesEC_cell,[]);
if matches > 0 && ~wild % If wildcard, ignore substrate match
    origin = 1;
    %If no match, try the closest organism but match the substrate:
else
    [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,'',true,false,...
        phylDist,org_index,SAcell,stringMatchesEC_cell,[]);
    if matches > 0 && ~wild % If wildcard, ignore substrate match
        origin = 2;
        %If no match, try to match organism but with any substrate:
    else
        [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,name,false,false,...
            phylDist,org_index,SAcell,stringMatchesEC_cell,[]);
        if matches > 0
            origin = 3;
            %If no match, try to match organism but for any substrate (SA*MW):
        else
            %create matching index for SA, has not been needed until now
            stringMatchesSA = extract_string_matches(EC,SAcell{1},wild,[],[]);

            [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,name,false,...
                true,phylDist,org_index,...
                SAcell,stringMatchesEC_cell,stringMatchesSA);
            if matches > 0
                origin = 4;
                %If no match, try any organism and any substrate:
            else
                [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,'',false,...
                    false,phylDist,...
                    org_index,SAcell,stringMatchesEC_cell,stringMatchesSA);
                if matches > 0
                    origin = 5;
                    %Again if no match, look for any org and SA*MW
                else
                    [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,'',...
                        false,true,phylDist,...
                        org_index,SAcell,stringMatchesEC_cell,stringMatchesSA);
                    if matches > 0
                        origin = 6;
                    end
                end

            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,organism,...
    substrate,SA,phylDist,...
    org_index,SAcell,KCATcellMatches,SAcellMatches)

%Will go through BRENDA and will record any match. Afterwards, it will
%return the average value and the number of matches attained.
kcat    = [];
matches = 0;

if SA
    %SAcell{1},wild,[],[]
    EC_indexes = extract_indexes(SAcellMatches,[],SAcell{2},subs,substrate,...
        organism,org_index,phylDist);

    kcat       = SAcell{3}(EC_indexes);
    org_cell   = SAcell{2}(EC_indexes);
    MW_BRENDA  = SAcell{4}(EC_indexes);

else
    %KCATcell{1},wild,ECIndexIds,EcIndexIndices
    EC_indexes = extract_indexes(KCATcellMatches,KCATcell{2},KCATcell{3},...
        subs,substrate,organism,org_index,...
        phylDist);
    if substrate
        for j = 1:length(EC_indexes)
            indx = EC_indexes(j);
            for k = 1:length(subs)
                if (isempty(subs{k}))
                    break;
                end
                %l = logical(strcmpi(model.metNames,subs{k}).*(model.S(:,i)~=0)); %I don't understand the .* (model.S(:,i)~=0) part, it shouldn't be needed?/JG;
                if ~isempty(subs{k}) && strcmpi(subs{k},KCATcell{2}(indx))
                    if KCATcell{4}(indx) > 0
                        coeff = min(substrCoeff);
                        kCatTmp = KCATcell{4}(indx);
                        kcat  = [kcat;kCatTmp/coeff];
                    end
                end
            end
        end
    else
        kcat = KCATcell{4}(EC_indexes);
    end
end
%Return maximum value:
if isempty(kcat)
    kcat = 0;
else
    matches        = length(kcat);
    [kcat,MaxIndx] = max(kcat);
end
%Avoid SA*Mw values over the diffusion limit rate  [Bar-Even et al. 2011]
if kcat>(1E7)
    kcat = 1E7;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the string matches of the ECs. This is heavy, so only do it once!
%
function EC_indexes = extract_string_matches(EC,EC_cell,wild,ECIndexIds,EcIndexIndices)
EC_indexes = [];
EC_indexesOld = [];
if wild
    if (~isempty(ECIndexIds)) %In some cases the EC_cell is not from KCatCell
        X = find(contains(ECIndexIds, EC));
        for j = 1:length(X)
            EC_indexes = [EC_indexes,EcIndexIndices{X(j)}];
        end
    else %Not optimized
        for j=1:length(EC_cell)
            if strfind(EC_cell{j},EC)==1
                EC_indexes = [EC_indexes,j];
            end
        end
    end
else
    if (~isempty(ECIndexIds)) %In some cases the EC_cell is not from KCatCell
        mtch = find(strcmpi(EC,ECIndexIds));
        if (~isempty(mtch))
            EC_indexes = EcIndexIndices{mtch};
        end
    else %%Not optimized
        if ~isempty(EC_cell)
            EC_indexes = transpose(find(strcmpi(EC,EC_cell)));
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract the indexes of the entries in the BRENDA data that meet the
%conditions specified by the search criteria
function EC_indexes = extract_indexes(EC_indCellStringMatches,subs_cell,orgs_cell,subs,...
    substrate,organism, org_index,...
    phylDist)

EC_indexes = EC_indCellStringMatches;%reuse so the string comparisons are not run many times

%If substrate=true then it will extract only the substrates appereances
%indexes in the EC subset from the BRENDA cell array

if substrate
    if (~isempty(EC_indexes)) %optimization
        Subs_indexes = [];
        for l = 1:length(subs)
            if (isempty(subs{l}))
                break;
            end
            Subs_indexes = horzcat(Subs_indexes,EC_indexes(strcmpi(subs(l),...
                subs_cell(EC_indexes))));
        end
        EC_indexes = Subs_indexes;
    end
end

EC_orgs = orgs_cell(EC_indexes);
%If specific organism values are requested looks for all the organism
%repetitions on the subset BRENDA cell array(EC_indexes)
if string(organism) ~= ''
    EC_indexes = EC_indexes(strcmpi(string(organism),EC_orgs));

    %If KEGG code was assigned to the organism (model) then it will look for
    %the Kcat value for the closest organism
elseif org_index~='*' %&& org_index~=''
    KEGG_indexes = [];temp = [];

    %For relating a phyl dist between the modelled organism and the organisms
    %on the BRENDA cell array it should search for a KEGG code for each of
    %these
    for j=1:length(EC_indexes)
        %Assigns a KEGG index for those found on the KEGG struct
        orgs_index = find(strcmpi(orgs_cell(EC_indexes(j)),phylDist.names),1);
        if ~isempty(orgs_index)
            KEGG_indexes = [KEGG_indexes; orgs_index];
            temp         = [temp;EC_indexes(j)];
            %For values related to organisms without KEGG code, then it will
            %look for KEGG code for the first organism with the same genus
        else
            org = orgs_cell{EC_indexes(j)};
            orgGenus = lower(regexprep(org,'\s.*',''));
            if isKey(phylDist.genusHashMap,orgGenus) %annoyingly, this seems to be needed
                matchInd = cell2mat(values(phylDist.genusHashMap,{orgGenus}));
                matches = phylDist.uniqueGenusIndices{matchInd};
                k = matches(1);
                KEGG_indexes = [KEGG_indexes;k];
                temp         = [temp;EC_indexes(j)];
            end
        end
    end
    %Update the EC_indexes cell array
    EC_indexes = temp;
    %Looks for the taxonomically closest organism and saves the index of
    %its appearences in the BRENDA cell
    if ~isempty(EC_indexes)
        distances = phylDist.distMat(org_index,KEGG_indexes);
        EC_indexes = EC_indexes(distances == min(distances));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function org_index = find_inKEGG(org_name,names)
org_index      = find(strcmpi(org_name,names));
if numel(org_index)>1
    org_index = org_index(1);
elseif isempty(org_index)
    org_name    = regexprep(org_name,'\s.*','');
    org_index   = find(strcmpi(org_name,names));
    if numel(org_index)>1
        org_index = org_index(1);
    elseif isempty(org_index)
        org_index = '*';
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phylDistStruct =  KEGG_struct(phylpath)
load(phylpath)
phylDistStruct.ids   = transpose(phylDistStruct.ids);
phylDistStruct.names = transpose(phylDistStruct.names);
phylDistStruct.names = regexprep(phylDistStruct.names,'\s*\(.*','');
end
