function MW = calculateMW(sequence)
% calculateMW  Calculate the molecular weight of a protein.
%
% Calculates the molecular weight of a protein from its amino acid sequence.
% Matching is case-insensitive. Sequences with no recognized residue
% (including an empty sequence) return NaN, rather than being mistaken for
% an 18 Da protein.
%
% Parameters
% ----------
% sequence : char
%     the amino acid sequence (can include extra things such as spaces or
%     numbers, which are ignored).
%
% Returns
% -------
% MW : double
%     the molecular weight in Da, or NaN if no residue was recognized.

% A	Alanine
% B	Aspartic acid or Asparagine
% C	Cysteine
% D	Aspartic acid
% E	Glutamic acid
% F	Phenylalanine
% G	Glycine
% H	Histidine
% I	Isoleucine
% J	Leucine or Isoleucine
% K	Lysine
% L	Leucine
% M	Methionine
% N	Asparagine
% O	Pyrrolysine
% P	Proline
% Q	Glutamine
% R	Arginine
% S	Serine
% T	Threonine
% U	Selenocysteine
% V	Valine
% W	Tryptophan
% X	any
% Y	Tyrosine
% Z	Glutamic acid or Glutamine

% The 20 standard residues, individually sourced.
aa_codes = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R', ...
            'S','T','V','W','Y'};
aa_MWs   = [71.08 103.14 115.09 129.11 147.17 57.05 137.14 113.16 128.17 ...
            113.16 131.20 114.10 97.12 128.13 156.19 87.08 101.10 99.13 ...
            186.21 163.17];

% Ambiguous / non-standard codes, derived from the standard table above
% rather than pre-rounded constants: B/Z average their two possible
% residues, X averages across all 20 standard residues, J is Leu-or-Ile
% (both 113.16, so either index gives the same value), O/U have no
% standard-residue equivalent and keep their literature values.
bMW = mean(aa_MWs(ismember(aa_codes,{'D','N'})));
zMW = mean(aa_MWs(ismember(aa_codes,{'E','Q'})));
xMW = mean(aa_MWs);
jMW = aa_MWs(strcmp(aa_codes,'L'));

allCodes = [aa_codes, {'B','J','O','U','X','Z'}];
allMWs   = [aa_MWs, bMW, jMW, 255.31, 150.04, xMW, zMW];

% Average mass of water (g/mol): the two H atoms and one OH that remain
% after peptide bond condensation.
waterMW = 18.01528;

sequence = upper(sequence);
MW = waterMW;
nMatched = 0;
for i = 1:length(allCodes)
    count = length(strfind(sequence,allCodes{i}));
    MW = MW + count*allMWs(i);
    nMatched = nMatched + count;
end
if nMatched == 0
    MW = NaN;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
