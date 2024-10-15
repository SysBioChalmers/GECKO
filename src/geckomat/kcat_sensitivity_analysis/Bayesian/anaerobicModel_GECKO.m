%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anaerobicModel_GECKO.m
% Converts model to anaerobic
%
% Benjamin J. Sanchez
% Feiran Li - 2019-09-24
% Feiran Li - Last update: 2019-10-02 modify the order of changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = anaerobicModel_GECKO(model)

%1th change: Refit GAM and NGAM to exp. data, change biomass composition
GAM   = 58.1988;  %Data from Nissen et al. 1997
strain = strrep(model.id,' specific model genereted from panYeast','');
if strcmp(strain,'Candida_glabrata') | strcmp(strain,'Candida_parapsilosis')
 GAM = 30;% if 
end
%P     = 0.461;  %Data from Nissen et al. 1997
NGAM  = 0;      %Refit done in Jouthen et al. 2012

model = changeGAM(model,GAM,NGAM);
%model = scaleBioMass(model,'protein',P,'carbohydrate');

%2nd change: Removes the requirement of heme a in the biomass equation
%            (not used under aerobic conditions)
%mets = {'s_3714[c]','s_1198[c]','s_1203[c]','s_1207[c]','s_1212[c]','s_0529[c]'};
mets = {'s_3714','s_1198','s_1203','s_1207','s_1212','s_0529'};
[~,met_index] = ismember(mets,model.mets);
model.S(met_index,strcmp(model.rxns,'r_4598')) = 0;

%3st change: Changes media to anaerobic (no O2 uptake and allows sterol
%            and fatty acid exchanges)
model.lb(strcmp(model.rxns,'r_1992')) = 0;        %O2
model.lb(strcmp(model.rxns,'r_1757')) = -1000;    %ergosterol
model.lb(strcmp(model.rxns,'r_1915')) = -1000;    %lanosterol
model.lb(strcmp(model.rxns,'r_1994')) = -1000;    %palmitoleate
model.lb(strcmp(model.rxns,'r_2106')) = -1000;    %zymosterol
model.lb(strcmp(model.rxns,'r_2134')) = -1000;    %14-demethyllanosterol
model.lb(strcmp(model.rxns,'r_2137')) = -1000;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol
model.lb(strcmp(model.rxns,'r_2189')) = -1000;    %oleate

% extra media set up for ura1 original but also the anaerobic growth media;
if isfield(model,'id')
    strain = strrep(model.id,' specific model genereted from panYeast','');
species_onlyura9 = {'Alloascoidea_hylecoeti;Ambrosiozyma_kashinagacola;Ambrosiozyma_monospora;Arxula_adeninivorans;Ascoidea_asiatica;Ascoidea_rubescens;Ashbya_aceri;Aspergillus_nidulans;Babjeviella_inositovora;Brettanomyces_anomalus;Candida_albicans;Candida_apicola;Candida_arabinofermentans;Candida_auris;Candida_boidinii_JCM9604;Candida_carpophila;Candida_dubliniensis;Candida_glabrata;Candida_homilentoma;Candida_infanticola;Candida_intermedia;Candida_orthopsilosis;Candida_parapsilosis;Candida_sorboxylosa;Candida_succiphila;Candida_tanzawaensis;Candida_tenuis;Candida_tropicalis;Candida_versatilis;Clavispora_lusitaniae;Cyberlindnera_fabianii_JCM3601;Cyberlindnera_jadinii;Debaryomyces_hansenii;Dekkera_bruxellensis;Eremothecium_coryli;Eremothecium_cymbalariae;Eremothecium_gossypii;Eremothecium_sinecaudum;Geotrichum_candidum;Hanseniaspora_uvarum;Hanseniaspora_valbyensis;Hanseniaspora_vinae;Hyphopichia_burtonii;Komagataella_pastoris;Kuraishia_capsulata;Lipomyces_starkeyi;Lodderomyces_elongisporus;Metschnikowia_aberdeeniae;Metschnikowia_arizonensis;Metschnikowia_bicuspidata;Metschnikowia_borealis;Metschnikowia_bowlesiae;Metschnikowia_cerradonensis;Metschnikowia_continentalis;Metschnikowia_dekortum;Metschnikowia_drakensbergensis;Metschnikowia_hamakuensis;Metschnikowia_hawaiiensis;Metschnikowia_hibisci;Metschnikowia_ipomoeae;Metschnikowia_kamakouana;Metschnikowia_kipukae;Metschnikowia_lockheadii;Metschnikowia_matae;Metschnikowia_matae_maris;Metschnikowia_mauinuiana;Metschnikowia_proteae;Metschnikowia_santaceciliae;Metschnikowia_shivogae;Metschnikowia_similis;Meyerozyma_guilliermondii;Millerozyma_acaciae;Nadsonia_fulvescens_var_elongata;Nakaseomyces_bracarensis;Nakaseomyces_castellii;Nakaseomyces_delphensis;Nakaseomyces_nivariensis;Nakazawaea_peltata;Ogataea_methanolica;Ogataea_parapolymorpha;Ogataea_polymorpha;Pachysolen_tannophilus;Pichia_membranifaciens;Priceomyces_haplophilus;Saccharomycopsis_malanga;Saprochaete_clavata;Scheffersomyces_lignosus;Scheffersomyces_stipitis;Schizosaccharomyces_pombe;Spathaspora_arborariae;Spathaspora_girioi;Spathaspora_gorwiae;Spathaspora_hagerdaliae;Spathaspora_passalidarum;Sporopachydermia_quercuum;Starmerella_bombicola_JCM9596;Sugiyamaella_lignohabitans;Tortispora_caseinolytica;Vanderwaltozyma_polyspora;Wickerhamia_fluorescens;Wickerhamiella_domercqiae;Wickerhamomyces_anomalus;Wickerhamomyces_ciferrii;Yarrowia_deformans;Yarrowia_keelungensis;Yarrowia_lipolytica;yHMPu5000026124_Ogataea_henricii;yHMPu5000026137_Ambrosiozyma_ambrosiae;yHMPu5000026142_Citeromyces_matritensis;yHMPu5000026145_Ambrosiozyma_vanderkliftii;yHMPu5000026197_Brettanomyces_custersianus;yHMPu5000026274_Komagataella_populi;yHMPu5000034594_Starmera_quercuum;yHMPu5000034597_Candida_stellimalicola;yHMPu5000034604_Sporopachydermia_lactativora;yHMPu5000034605_Spencermartinsiella_europaea;yHMPu5000034606_Priceomyces_medius;yHMPu5000034607_Saccharomycopsis_capsularis;yHMPu5000034610_Saturnispora_hagleri;yHMPu5000034611_Saturnispora_mendoncae;yHMPu5000034612_Saturnispora_saitoi;yHMPu5000034613_Saturnispora_serradocipensis;yHMPu5000034614_Saturnispora_silvae;yHMPu5000034615_Saturnispora_zaruensis;yHMPu5000034622_Pichia_occidentalis;yHMPu5000034623_Pichia_norvegensis;yHMPu5000034624_Pichia_nakasei;yHMPu5000034625_Pichia_kudriavzevii;yHMPu5000034627_Pichia_heedii;yHMPu5000034629_Pichia_exigua;yHMPu5000034631_Martiniozyma_abiesophila;yHMPu5000034632_Candida_athensensis;yHMPu5000034635_Nadsonia_fulvescens;yHMPu5000034636_Ogataea_nitratoaversa;yHMPu5000034637_Ogataea_populiabae;yHMPu5000034643_Candida_schatavii;yHMPu5000034646_Wickerhamiella_cacticola;yHMPu5000034648_Candida_restingae;yHMPu5000034654_Aciculoconidium_aculeatum;yHMPu5000034655_Botryozyma_nematodophila;yHMPu5000034660_Diddensiella_caesifluorescens;yHMPu5000034661_Dipodascus_albidus;yHMPu5000034665_Kodamaea_laetipori;yHMPu5000034667_Blastobotrys_serpentis;yHMPu5000034669_Blastobotrys_raffinofermentans;yHMPu5000034670_Blastobotrys_proliferans;yHMPu5000034671_Blastobotrys_peoriensis;yHMPu5000034673_Blastobotrys_nivea;yHMPu5000034674_Blastobotrys_muscicola;yHMPu5000034675_Blastobotrys_mokoenaii;yHMPu5000034681_Blastobotrys_americana;yHMPu5000034742_Lipomyces_suomiensis;yHMPu5000034748_Lipomyces_oligophaga;yHMPu5000034749_Lipomyces_mesembrius;yHMPu5000034754_Lipomyces_arxii;yHMPu5000034760_Lipomyces_kononenkoae;yHMPu5000034761_Lipomyces_lipofer;yHMPu5000034883_Peterozyma_xylosa;yHMPu5000034884_Peterozyma_toletana;yHMPu5000034885_Ogataea_zsoltii;yHMPu5000034886_Ogataea_trehalophila;yHMPu5000034887_Ogataea_trehaloabstinens;yHMPu5000034890_Ogataea_ramenticola;yHMPu5000034891_Ogataea_pini;yHMPu5000034892_Ogataea_pilisensis;yHMPu5000034893_Ogataea_philodendra;yHMPu5000034897_Ogataea_glucozyma;yHMPu5000034899_Ogataea_kodamae;yHMPu5000034901_Ogataea_methylivora;yHMPu5000034902_Ogataea_minuta;yHMPu5000034903_Ogataea_naganishii;yHMPu5000034904_Ogataea_nonfermentans;yHMPu5000034918_Nakazawaea_holstii;yHMPu5000034933_Kuraishia_molischiana;yHMPu5000034939_Komagataella_pseudopastoris;yHMPu5000034946_Ambrosiozyma_oregonensis;yHMPu5000034947_Ambrosiozyma_philentoma;yHMPu5000034950_Citeromyces_hawaiiensis;yHMPu5000034952_Citeromyces_siamensis;yHMPu5000034957_Hanseniaspora_osmophila;yHMPu5000034963_Hanseniaspora_clermontiae;yHMPu5000034967_Candida_freyschussii;yHMPu5000034973_Danielozyma_ontarioensis;yHMPu5000034974_Deakozyma_indianensis;yHMPu5000034978_Cyberlindnera_mrakii;yHMPu5000034979_Cyberlindnera_misumaiensis;yHMPu5000034986_Candida_oregonensis;yHMPu5000034988_Candida_fructus;yHMPu5000034990_Candida_corydali;yHMPu5000034998_Cephaloascus_albidus;yHMPu5000034999_Cephaloascus_fragrans;yHMPu5000035011_Candida_pyralidae;yHMPu5000035018_Candida_canberraensis;yHMPu5000035022_Candida_emberorum;yHMPu5000035031_Candida_kruisii;yHMPu5000035032_Candida_gatunensis;yHMPu5000035033_Candida_cretensis;yHMPu5000035037_Candida_montana;yHMPu5000035040_Ambrosiozyma_maleeae;yHMPu5000035041_Ambrosiozyma_pseudovanderkliftii;yHMPu5000035044_Barnettozyma_californica;yHMPu5000035045_Barnettozyma_hawaiiensis;yHMPu5000035046_Barnettozyma_populi;yHMPu5000035047_Barnettozyma_pratensis;yHMPu5000035048_Barnettozyma_salicaria;yHMPu5000035242_Zygoascus_ofunaensis;yHMPu5000035243_Zygoascus_meyerae;yHMPu5000035244_Candida_incommunis;yHMPu5000035252_Yamadazyma_nakazawae;yHMPu5000035261_Candida_ponderosae;yHMPu5000035268_Wickerhamomyces_hampshirensis;yHMPu5000035271_Wickerhamomyces_bovis;yHMPu5000035274_Wickerhamomyces_alni;yHMPu5000035279_Tortispora_starmeri;yHMPu5000035282_Trigonopsis_vinaria;yHMPu5000035286_Candida_azyma;yHMPu5000035296_Priceomyces_carsonii;yHMPu5000035297_Priceomyces_castillae;yHMPu5000035301_Pichia_terricola;yHMPu5000035302_Candida_fragi;yHMPu5000035318_Hyphopichia_heimii;yHMPu5000035325_Cyberlindnera_petersonii;yHMPu5000035335_Candida_blattae;yHMPu5000035629_Yueomyces_sinensis;yHMPu5000035633_Candida_hispaniensis;yHMPu5000035639_Wickerhamomyces_canadensis;yHMPu5000035640_Yamadazyma_philogaea;yHMPu5000035641_Yamadazyma_scolyti;yHMPu5000035643_Yarrowia_bubula;yHMPu5000035645_Yarrowia_divulgata;yHMPu5000035650_Trigonopsis_variabilis;yHMPu5000035654_Tortispora_ganteri;yHMPu5000035658_Starmera_amethionina;yHMPu5000035659_Saturnispora_dispora;yHMPu5000035662_Meyerozyma_caribbica;yHMPu5000035665_Middelhovenomyces_tepae;yHMPu5000035667_Kurtzmaniella_cleridarum;yHMPu5000035670_Phaffomyces_opuntiae;yHMPu5000035671_Phaffomyces_antillensis;yHMPu5000035672_Phaffomyces_thermotolerans;yHMPu5000035673_Candida_orba;yHMPu5000035674_Kregervanrija_delftensis;yHMPu5000035675_Kregervanrija_fluxuum;yHMPu5000035677_Kodamaea_ohmeri;yHMPu5000035679_Candida_rhagii;yHMPu5000035681_Candida_gotoi;yHMPu5000035684_Kloeckera_hatyaiensis;yHMPu5000035686_Cyberlindnera_saturnus;yHMPu5000035687_Cyberlindnera_suaveolens;yHMPu5000035688_Cyberlindnera_xylosilytica;yHMPu5000035689_Candida_mycetangii;yHMPu5000035690_Candida_vartiovaarae;yHMPu5000035691_Candida_salmanticensis;yHMPu5000035695_Hanseniaspora_pseudoguilliermondii;yHMPu5000035696_Hanseniaspora_singularis;yHMPu5000035699_Cyberlindnera_maclurae;yHMPu5000035703_Cyberlindnera_americana;yHMPu5000035707_Candida_heveicola;yHMPu5000041678_Debaryomyces_prosopidis;yHMPu5000041693_Debaryomyces_nepalensis;yHMPu5000041713_Debaryomyces_maramus;yHMPu5000041743_Candida_hawaiiana;yHMPu5000041818_Magnusiomyces_tetrasperma;yHMPu5000041822_Dipodascus_geniculatus;yHMPu5000041824_Debaryomyces_subglobosus;yHMPu5000041829_Debaryomyces_fabryi;yHMPu5000041833_Candida_tammaniensis;yHMPu5000041840_Candida_wancherniae;yHMPu5000041855_Candida_ascalaphidarum;yHMPu5000041862_Candida_golubevii;yHMPu5000041863_Candida_gorgasii'};
species_onlyura9 = split(species_onlyura9,';');
anaerobic = {'Sugiyamaella_lignohabitans';'Dekkera_bruxellensis';'yHMPu5000034625_Pichia_kudriavzevii';'yHMPu5000026142_Citeromyces_matritensis';'Candida_albicans';'Candida_parapsilosis';'Candida_tropicalis';...
             'Clavispora_lusitaniae';'Spathaspora_passalidarum';'Wickerhamia_fluorescens';'Wickerhamomyces_anomalus';'yHMPu5000035686_Cyberlindnera_saturnus';'Hanseniaspora_uvarum';'Hanseniaspora_valbyensis';...
             'Hanseniaspora_vinae';'yHMPu5000034957_Hanseniaspora_osmophila';'Ashbya_aceri';'Candida_glabrata';'Eremothecium_coryli';'Kluyveromyces_lactis';'Kluyveromyces_marxianus';'Lachancea_fermentati';...
             'Lachancea_kluyveri';'Lachancea_thermotolerans';'Lachancea_waltii';'Nakaseomyces_bacillisporus';'Nakaseomyces_castellii';'Nakaseomyces_delphensis';'Naumovozyma_castellii';'Naumovozyma_dairenensis';...
             'Saccharomyces_cerevisiae';'Saccharomyces_eubayanus';'Saccharomyces_paradoxus';'Saccharomyces_uvarum';'Tetrapisispora_blattae';'Tetrapisispora_phaffii';'Torulaspora_delbrueckii';'Vanderwaltozyma_polyspora';...
             'Zygosaccharomyces_bailii';'yHAB154_Kazachstania_transvaalensis';'yHMPu5000034881_Torulaspora_pretoriensis';'yHMPu5000034876_Tetrapisispora_iriomotensis';'yHMPu5000034862_Zygotorulaspora_florentina';...
             'yHMPu5000026152_Torulaspora_franciscae';'Schizosaccharomyces_pombe'};
species_intersect = intersect(species_onlyura9,anaerobic);

if ismember(strain,species_intersect)
    model = changeRxnBounds(model,'r_2090',-1000,'l');% uracil uptake
end
end
%4rd change: Blocked pathways for proper glycerol production
%Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
model.lb(startsWith(model.rxns,'r_0713')) = 0; %Mithocondria
model.lb(strcmp(model.rxns,'r_0714')) = 0; %Cytoplasm
model.lb(strcmp(model.rxns,'r_0713_rvs')) = 0; %Mithocondria
model.lb(strcmp(model.rxns,'r_0714')) = 0; %Cytoplasm
model.lb(strcmp(model.rxns,'r_0714_rvs')) = 0; %Cytoplasm
%Block glycerol dehydroginase (only acts in microaerobic conditions)
model.ub(strcmp(model.rxns,'r_0487')) = 0;
model.ub(strcmp(model.rxns,'r_0487_rvs')) = 0;
%Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
model.ub(strcmp(model.rxns,'r_0472')) = 0;
model.ub(strcmp(model.rxns,'r_0472_fwd')) = 0;


%4th change: Blocked pathways for proper glycerol production
%Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
idx = regexp(cellstr(model.rxns),'r_0713[\w*]rvs');
idx = find(cellfun(@isempty,idx)==0);
model = changeRxnBounds(model,model.rxns(idx),0,'b'); 
model.lb(strcmp(model.rxns,'r_0713')) = 0; %Mithocondria % in case this one does not have any grRule

idx = regexp(cellstr(model.rxns),'r_0714[\w*]rvs');
idx = find(cellfun(@isempty,idx)==0);
model = changeRxnBounds(model,model.rxns(idx),0,'b'); 
model.lb(strcmp(model.rxns,'r_0714')) = 0; %Cytoplasm

%Block glycerol dehydroginase (only acts in microaerobic conditions)
% model.ub(strcmp(model.rxns,'r_0487')) = 0; is blocked in blockrxns

%Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
idx = find(startsWith(model.rxns,'r_0472_'));
model = changeRxnBounds(model,model.rxns(idx),0,'b');
model.ub(strcmp(model.rxns,'r_0472')) = 0;

end

%%

function model = changeGAM(model,GAM,NGAM)

bioPos = strcmp(model.rxnNames,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,bioPos);
    isGAM = sum(strcmp({'ATP [cytoplasm]','ADP [cytoplasm]','H2O [cytoplasm]', ...
        'H+ [cytoplasm]','phosphate [cytoplasm]'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,bioPos) = sign(S_ix)*GAM;
    end
end

if nargin >1
    pos = strcmp(model.rxnNames,'non-growth associated maintenance reaction');%NGAM
    %model = setParam(model,'eq',model.rxns(pos),NGAM);% set both lb and ub to be mu
    model.lb(pos) = NGAM;
    model.ub(pos) = NGAM;
    
end

end
