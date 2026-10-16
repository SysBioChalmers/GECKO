function n = countRxns(ecRxns)
% countRxns  Number of reactions an ecRxns argument selects (logical or indices).
if islogical(ecRxns)
    n = sum(ecRxns);
else
    n = numel(unique(ecRxns));
end
end
