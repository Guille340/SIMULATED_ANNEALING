function isTrue = isOptionsStruct(S)

isTrue = false;
if isstruct(S)
    fieldNames_valid = fieldnames(initialiseOptionsStruct);
    fieldNames = fieldnames(S);
    if all(ismember(fieldNames_valid,fieldNames)) ...
            && all(ismember(fieldNames,fieldNames_valid))
        isTrue = true;
    end
end
    