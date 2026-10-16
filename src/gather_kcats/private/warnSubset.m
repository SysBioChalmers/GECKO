function warnSubset(ecRxns, nRxns, file)
% warnSubset  Say that gathering for some reactions replaces the whole file.
if ~isempty(ecRxns) && countRxns(ecRxns) < nRxns && isfile(file)
    printOrange(sprintf(['%s is replaced by kcats for %d reactions only. Give saveAs ' ...
        'to keep it.\n'], file, countRxns(ecRxns)));
end
end
