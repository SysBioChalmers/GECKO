%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function genesSets = getSimpleGeneSets(originalSTR)
% function that gets a cell array with all the simple geneSets in a given 
% grRule string 
%
% Last edited  Ivan Domenzain, 2018-03-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function genesSets = getSimpleGeneSets(originalSTR)
    genesSets = cell(1,1);
    % If gene rule is not empty split in all its different isoenzymes
    if ~isempty(originalSTR)
        %Remove all brackets
        originalSTR = strrep(originalSTR,'(','');
        originalSTR = strrep(originalSTR,')','');
        %Split all the different genesSets
        genesSets  = transpose(strsplit(originalSTR,' or '));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%