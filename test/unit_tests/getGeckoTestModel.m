function testModel = getGeckoTestModel()

testModel = struct();
testModel.name = 'testModel';
testModel.id = 'testModel';
testModel.annotation.defaultLB = -Inf;
testModel.annotation.defaultUB = Inf;
testModel.date = datestr(now,29);
testModel.description = '';
testModel.version = '';
testModel.rxns = {};
testModel.S=[];
testModel.rev=[];
testModel.metNames = {'e1';'e2';'m1';'m2'};
testModel.comps = {'e';'c'};
testModel.compNames = testModel.comps;
testModel.metComps = [1;1;2;2];
testModel.mets = strcat(testModel.metNames, testModel.comps(testModel.metComps));
testModel.grRules = {};
testModel.rxnGeneMat = [];

testModel.genes = {'G1';'G2';'G3';'G4';'G5'};


testModel.ub = [];
testModel.lb = [];

rxnsToAdd = struct();
rxnsToAdd.rxns = {  'S1';'R1';'R2';'R3';'R4';'R5';'S2'};
rxnsToAdd.grRules = {'';  '';  'G1 and G2 or G3';  'G4';  ''; 'G5';  ''};
rxnsToAdd.equations = {'e1[e] <=>';...
                       'e1[e] <=> m1[c]';...
                       'm1[c] <=> m2[c]';...
                       'm1[c] => m2[c]';...
                       'm1[c] => m2[c]';...
                       'm2[c] => e2[e]';...
                       'e2[e] <=>'
                       };
testModel = addRxns(testModel,rxnsToAdd, 3);
testModel.c = [0;0;0;0;0;1;0];%optimize for output flux, if this is used
testModel.rxnNames = testModel.rxns;
testModel.b = repmat(0,length(testModel.mets),1);
testModel.eccodes = {'';'';'1.1.1.1';'1.1.2.1';'';'1.1.1.3';''};

end


