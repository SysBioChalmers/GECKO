function [pIDs,protData,fermParameters,byProds] = load_Prot_Ferm_Data(grouping)
%load_Prot_Ferm_Data
%
% Function that loads absolute proteomics data, total protein measurements and 
% fermentation data (exchange fluxes for Glucose, O2, CO2)
%
%   grouping      (vector) Number of biological replicates for each
%                 experimental condition in the dataset.
%
%   pIDs           (cell) With Uniprot IDs for the proteomics dataset
%   abs_data       (cell) Numerical values in the proteomics dataset
%   fermParameters (cell) Contains numerical values Ptot [g /gDw], dilution 
%                  rate [1/h], GUR, OUR, and CO2 production [mmol/gDw h]
%                  from the chemostat cultures.
%   byProds        (cell) metNames for the measured byProducts in the
%                  fermentations (HPLC)            
%
% Usage: [pIDs,protData,fermData,byProds] = load_Prot_Ferm_Data(grouping)
%
% Last modified.  Ivan Domenzain 2019-09-10

%Define dataset format
format = '%s %s';
for i=1:sum(grouping)
    format = [format ' %f'];
end
fileName = '../../../Databases/abs_proteomics.txt';
fID      = fopen(fileName);
protData = textscan(fID,format,'Delimiter','\t','HeaderLines',1,'TreatAsEmpty',{'NA','na','NaN'});
pIDs     = protData{1};
protData = protData(3:end);
%Load total protein content and fermentation data
fileName  = '../../../Databases/fermentationData.txt';
fID       = fopen(fileName);
formatStr = '%s';
data      = textscan(fID,formatStr,'Delimiter','\n');
fermData  = [];
for i=1:length(data{1})
    row      = data{1}(i);
    row      = strsplit(row{1},'\t');
    row      = row(1:end);
    fermData = [fermData; row]; 
end
%Extract observed byProduct names from file
[~,n] = size(fermData);
if n>6
    byProds  = fermData(1,7:end);
    byProds  = strrep(byProds,' (mmol/gDw h)','');
    byProds  = strrep(byProds,' [mmol/gDw h]','');
else
    byProds  = [];
end
fermParameters           = [];
fermParameters.conds     = fermData(2:end,1);
fermParameters.Ptot      = str2double(fermData(2:end,2));
fermParameters.Drate     = str2double(fermData(2:end,3));
fermParameters.GUR       = str2double(fermData(2:end,4));
fermParameters.CO2prod   = str2double(fermData(2:end,5));
fermParameters.OxyUptake = str2double(fermData(2:end,6));
fermParameters.byP_flux  = byProds;
if ~isempty(fermParameters.byP_flux)
    fermParameters.byP_flux = str2double(fermData(2:end,7:end));
end
end

