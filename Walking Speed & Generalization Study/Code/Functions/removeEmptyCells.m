function outCA = removeEmptyCells(inCA)
    outCA =  inCA(~cellfun('isempty',inCA)); 
end