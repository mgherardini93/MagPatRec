% -------------------------- Function Description -------------------------
% This function will extract all the selected features and store them in
% the exFeatures variable within the data structure.
% If the transient features are selected, it will extract features from
% full data, which is sotred in dataLA variable.
%
% input: 
% - data structure which contains dataTreated and dataLA
% - selected features in a form of a cell row
%
% output:
% - data structure which contains exFeatures 
% - message 
%
function [dataStruct, message] = ExtractFeatures(dataStruct, selFeatures)

try
    cTp = dataStruct.cTp;
    cT = dataStruct.cT;
    nR = dataStruct.nR;
    nM = dataStruct.nM;
    data = dataStruct.dataTreated;
    % First calculate the frequency if it is unknown
    if (dataStruct.sF == 0)
        lengthData = size(data{1},1);
        timeOfContractions = nR*cT*cTp;
        dataStruct.sF = round(lengthData/timeOfContractions);
    end
    sF = dataStruct.sF;
    wSamples = round(dataStruct.tW*sF); % number of samples that correspond to a time window
    iSamples = round(dataStruct.wIncrement*sF); % number of samples for the time window increment
    
    exFeatures = cell(1,nM);
    
    for m = 1:nM
        idx = 1; % index which indicates the beggining of the time window
        idxFeatures = 1;
        lengthMovement = size(data{m},1);
        while (idx <= (lengthMovement - wSamples + 1))
            windowData = data{m}(idx:(idx+wSamples-1),:); % data that belong to one time window
            exFeatures{m}(idxFeatures,:) = GetSigFeatures(windowData,sF,[],selFeatures');
            idx = idx + iSamples;
            idxFeatures = idxFeatures + 1;
        end
        disp("Number of samples not considered: " + num2str(lengthMovement - (idx - iSamples + wSamples - 1))); % display the number of samples which were not used for extracting the time windows
    end
    dataStruct.exFeatures = exFeatures;
    dataStruct.selFeatures = selFeatures;

    if (sum(contains(selFeatures,'Tr'))) % in case transient features are selected
        % the whole data needs to be considered
        data = dataStruct.dataLA;
        label = dataStruct.label;
        nM = size(data,2);
        exFeaturesTrans = cell(1,nM);
        labelFeatures = cell(1,nM);
        for m = 1:nM
            idx = 1; % index which indicates the beggining of the time window
            idxFeatures = 1;
            lengthMovement = size(data{m},1);
            while (idx <= (lengthMovement - wSamples + 1))
                windowData = data{m}(idx:(idx+wSamples-1),:); % data that belong to one time window
                windowLabel = label{m}(idx:(idx+wSamples-1),:);
                exFeaturesTrans{m}(idxFeatures,:) = GetSigFeatures(windowData,sF,[],selFeatures');
                labelFeatures{m}(idxFeatures,:) = mode(windowLabel);
                idx = idx + iSamples;
                idxFeatures = idxFeatures + 1;
            end
            disp("Number of samples not considered: " + num2str(lengthMovement - (idx - iSamples + wSamples - 1))); % display the number of samples which were not used for extracting the time windows
        end
        dataStruct.exFeaturesTrans = exFeaturesTrans;
        dataStruct.labelFeatures = labelFeatures;
    end

    message = 'Features successfully extracted.';
    disp('Features extracted')
catch
    message = 'Error: Features not extracted!';
end

end