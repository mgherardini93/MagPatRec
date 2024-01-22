% -------------------------- Function Description -------------------------
% This function will extract the middle part of the window (whether it be
% contraction period or rest period). The exctraction part depends on
% variable cTp, which is set as a percentage
%
% input: 
% - data structure which contains the data with calculated distances and
% labels which indicate when is rest and when contraction
% - cTp is the percentage of the window that is exctracted
%
% output:
% - data structure which contains dataTreated and labelTreated, which has 
% reduced size with respect to dataLA
% - message 
%
function [dataStruct, message] = EffectiveContractionExtraction(dataStruct, cTp)

message = 'Error';
data = dataStruct.dataLA;
label = dataStruct.label;
nM = size(dataStruct.rawData,2); % number of movements
dataTreated = cell(1,nM);
labelTreated = cell(1,nM);

dataStruct.cTp = cTp;
for m = 1:nM
    diffLabel = diff(label{m}); % find where the contraction starts and ends
    % note: rest is always labeled as 0
    indxRep = find(diffLabel ~= 0);
    indxRep = [0; indxRep; length(label{m})];
    % for each segment (whether it is rest or contraction), find which part
    % has to be discarded
    for s = 1:(length(indxRep)-1)
        lengthSegment = indxRep(s+1) - indxRep(s);
        effectivePart = round(lengthSegment*cTp); % number of samples that is kept
        discardPart = round((lengthSegment - effectivePart)/2); % number of samples that is descarded at the beggining and at the end of the segment
        indxStart = indxRep(s) + 1 + discardPart;
        indxEnd = indxRep(s) + effectivePart + discardPart;
        dataTreated{m} = [dataTreated{m}; data{m}(indxStart:indxEnd,:)];
        labelTreated{m} = [labelTreated{m}; label{m}(indxStart:indxEnd)];
    end
end

if(isempty(dataTreated{1}) || isempty(labelTreated{1}))
    message = 'Error: Data not extracted correctly!';
    return;
end

dataStruct.dataTreated = dataTreated;
dataStruct.labelTreated = labelTreated;
message = 'Data successfully extracted.';

end