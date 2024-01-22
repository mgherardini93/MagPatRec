function [patRec, message] = OfflineTraining(dataStruct, patRec)

movUsed = patRec.mov;
nM = length(movUsed);
patRec.nM = nM;
selFeatures = patRec.selFeatures;
magSelFeatures = patRec.magSelFeatures;
[~, magnetPairsIdx] = ismember(patRec.magnetPairs, dataStruct.magnetPairs,'rows'); 

% Checking the selected data ----------------------------------------------
if (isempty(magnetPairsIdx))
    message = 'Error: No magnet pairs are selected!';
    return;
end
if (contains(patRec.algorithm.type,'Select'))
    message = 'Error: Classifier not chosen!';
    return;
end
if (contains(patRec.algorithm.kernel,'Select') && ~strcmp(patRec.algorithm.type,'MLR'))
    message = 'Error: Classifier kernel not chosen!';
    return;
end
if (contains(patRec.normalization.type,'Select'))
    message = 'Error: Normalization method not chosen!';
    return;
end
if (contains(patRec.topology,'Select'))
    message = 'Error: Topology not chosen!';
    return;
end
if (isempty(selFeatures))
    message = 'Error: No features are selected!';
    return;
end
if (isempty(magSelFeatures))
    message = 'Error: No magnet features are selected!';
    return;
end

% Separate selected classes -----------------------------------------------
if(patRec.algorithm.transient || strcmp(patRec.algorithm.type,'MLR')) % transient features
    if (isfield(dataStruct,'exFeaturesTrans')) % check if transient features are extracted
        if(~contains(movUsed,'Rest'))
            featuresStruct = dataStruct.exFeaturesTrans(ismember(dataStruct.mov,movUsed));
        else
            message = 'Error: Remove "Rest" class from training!';
            return;
        end
    else
        message = 'Error: No transient features extracted!';
        return;
    end
else
    featuresStruct = dataStruct.exFeatures(ismember(dataStruct.mov,movUsed));
end
if (isempty(featuresStruct))   % steady-state features
    message = 'Error: No movements selected!';
    return;
end

dataLength = 0;
for m = 1:nM
    dataLength = dataLength + length(featuresStruct{m});
end

% Based on magnet selected features and magnet pairs that are used, find
% the indexes of channels that need to be used ----------------------------
if (sum(magSelFeatures) > sum(dataStruct.magSelFeatures))
    message = 'Error: Magnet feature selected is not extracted!';
    return;
end
if((sum(magSelFeatures)==sum(dataStruct.magSelFeatures))&&(sum(magSelFeatures)==1)) % in case only one magnet feature is extracted and selected
    chIdx = magnetPairsIdx;
end
if (sum(dataStruct.magSelFeatures)==2) % in case both magnet features are extracted
    if (sum(magSelFeatures)==2)        % in case both magnet features are selected
        chIdx = [magnetPairsIdx*2; magnetPairsIdx*2-1];
        chIdx = sort(chIdx);
    else                               % in case only one magnet feature is selected
        if (magSelFeatures(1))         % linear dist
            chIdx = magnetPairsIdx*2-1;
        elseif(magSelFeatures(2))      % angular dist
            chIdx = magnetPairsIdx*2;
        else
            message = 'Error: No magnet features selected!';
            return;
        end
    end
end

% Stack selected features in a matrix -------------------------------------
nCh = length(chIdx);
patRec.nCh = nCh;
patRec.vCh = chIdx;
dataSet = zeros(dataLength,nCh*length(selFeatures));
dataOut = zeros(dataLength,nM);
dSIdx = 0; % index that indicated at which point to enter data in dataSet
for m = 1:nM
    movementLength = length(featuresStruct{m});
    for s = 1:movementLength
        for f = 1:length(selFeatures)
            dataSet(dSIdx+s,(1:nCh) + (f-1)*nCh) = featuresStruct{m}(s).(selFeatures{f})(chIdx); 
        end
    end
    if (strcmp(patRec.algorithm.type,'MLR'))
        dataOut(:,m) = -1*ones(dataLength,1);
        dataOut((dSIdx+1):(dSIdx + movementLength),m) = dataStruct.labelFeatures{m};
    else
        dataOut((dSIdx+1):(dSIdx + movementLength),m) = 1;
    end
    dSIdx = dSIdx + movementLength;
end
patRec.nR = dataStruct.nR;
patRec.minData = dataStruct.restDataLA;

% Randomize data set if needed --------------------------------------------
if (patRec.randomize)
    if (patRec.algorithm.transient)
        message = 'Error: Transent and Randomize data are both selected!';
        return;
    end
    indD = 0;
    for m = 1:nM
        movementLength = sum(dataOut(:,m));
        randInd = randperm(movementLength);
        dataSet((1:movementLength) + indD,:) = dataSet(randInd + indD,:);
        indD = indD + movementLength;
    end
end

% Train and test data depending on the chosen topology --------------------
switch patRec.topology
    case 'Single Classifier'
        [patRec, message] = SingleClassifier(patRec,dataSet,dataOut);
    case 'Cross-validation'
        [patRec, message] = CrossValidationClassifier(patRec,dataSet,dataOut);
    case 'Single Classifier - test full'
        [patRec, message] = SingleClassifierTestFull(patRec,dataSet,dataOut,dataStruct);
    case 'Cross-validation - test full'
        [patRec, message] = CrossValidationClassifierTestFull(patRec,dataSet,dataOut,dataStruct);
    case 'Bottleneck test'
        [patRec, message] = BottleneckTest(patRec,dataSet,dataOut);
    case 'Train model'
        [patRec, message] = SingleClassifier(patRec,dataSet,dataOut);
        if (contains(message,'Error'))
            return;
        end
        [patRec, message] = TrainModel(patRec,dataSet,dataOut);
    case 'Test model'
        [patRec, message] = TestModel(patRec,dataSet,dataOut);
    otherwise
        message = 'Error: Selected topology not found';
        return;
end

if (contains(message,'Error'))
    return;
end

patRec = rmfield(patRec,'plotConfMat');
patRec = rmfield(patRec,'plotFigures');
patRec = rmfield(patRec,'randomize');

% message = 'Classifer trained and tested successfully';

end