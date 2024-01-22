function [patRec, message] = SingleClassifierTestFull(patRec, dataSet, dataOut, dataStruct)

movUsed = patRec.mov;
if(sum(strcmp(movUsed,'Rest')))
    movUsed(strcmp(movUsed,'Rest')) = [];
end
featuresStruct = dataStruct.exFeaturesTrans(ismember(dataStruct.mov,movUsed));
tP = patRec.trainingPercentage;
dataLength = size(dataSet,1);
nCh = size(dataSet,2);      % number of channels and features
Ntr = round(tP*dataLength); % number of samples for training
nM = patRec.nM;
nMt = length(movUsed);

trSet = zeros(Ntr,nCh); % training set
trOut = zeros(Ntr,nM);  % training out
idxD = 0;  % index to follow the dataSet
idxDT = 0; % index to follow the testDataSet
idxTr = 0; % index to follow trSet
idxT = 0;  % index to follow tSet

% Stack selected features in a matrix for test ----------------------------
nPair = patRec.nCh;
chIdx = patRec.vCh;
testDataLength = 0;
for m = 1:nMt
    testDataLength = testDataLength + length(featuresStruct{m});
end
testDataSet = zeros(testDataLength,nCh);
testDataOut = zeros(testDataLength,nM);
dSIdx = 0; % index that indicated at which point to enter data in dataSet
for m = 1:nMt
    movementLength = length(featuresStruct{m});
    for s = 1:movementLength
        for f = 1:length(patRec.selFeatures)
            testDataSet(dSIdx+s,(1:nPair) + (f-1)*nPair) = featuresStruct{m}(s).(patRec.selFeatures{f})(chIdx); 
        end
    end
    label = dataStruct.label{m};
    labelLength = length(label);
    if (labelLength > movementLength)
        ratio = floor(labelLength/movementLength);
        labelDownsampled = downsample(label,ratio);
        labelLength = length(labelDownsampled);
        if (labelLength > movementLength)
            samplesToRemove = labelLength - movementLength;
            samplesToRemoveIdx = randi(labelLength,1,samplesToRemove);
            labelDownsampled(samplesToRemoveIdx) = [];
        end
        contractionIdx = find(labelDownsampled);
        restIdx = find(labelDownsampled == 0);
        testDataOut(contractionIdx + dSIdx,m) = 1;
        testDataOut(restIdx + dSIdx, nM) = 1;
    end
%     testDataOut((dSIdx+1):(dSIdx + movementLength),m) = 1;
    dSIdx = dSIdx + movementLength;
end
Nt = round((1-tP)*testDataLength); % number of samples for testing
tSet = zeros(Nt,nCh);              % test set
tOut = zeros(Nt,nM);               % test out

% Separate training form test data ----------------------------------------
for m = 1:nM
    movementLength = sum(dataOut(:,m) == 1); % find which data correspond to the movement
    trLength = round(movementLength*tP);     % length of the data that goes into training 
    trSet((1:trLength)+idxTr,:) = dataSet((1:trLength)+idxD,:);
    trOut((1:trLength)+idxTr,m) = 1;
    % increase indeces to reach the next movement
    idxD = idxD + movementLength; 
    idxTr = idxTr + trLength;     
end
for m = 1:nMt
    movementLength = length(featuresStruct{m}); % find which data correspond to the movement
    tLength = round(movementLength*(1 - tP));     % length of the data that goes into test
    tSet((1:tLength)+idxT,:) = testDataSet((1:tLength)+idxDT+(movementLength - tLength),:);
    tOut((1:tLength)+idxT,:) = testDataOut((1:tLength)+idxDT+(movementLength - tLength),:);  
    % increase indeces to reach the next movement
    idxDT = idxDT + movementLength;   
    idxT = idxT + tLength; 
end
Ntr = idxTr;
Nt = idxT;
if (Nt < size(tOut,1))
    tOut = tOut(1:Nt,:);
    tSet = tSet(1:Nt,:);
end
if (Ntr < size(trOut,1))
    trOut = trOut(1:Ntr,:);
    trSet = trSet(1:Ntr,:);
end


% Train and test the classifier -------------------------------------------
[normParams, model, performance, message] = ClassifierTrainAndTest(patRec,trSet,trOut,tSet,tOut);

if (contains(message,'Error'))
    return;
end

if (patRec.plotFigures)
    figure('Name','Continuous performance of the classifier');
    plot(-performance.prediction,'b-')
    hold on
    [label,~] = find(tOut');
    plot(-label,'r-')
    legend('Classifier predicitons','Ground truth');
    hold off
    xlabel('Sampeles')
    ylabel('Classes')
    yticks(-patRec.nM:-1)
    yticklabels(patRec.mov(end:-1:1));
end

patRec.normalization.param1 = normParams.param1;
patRec.normalization.param2 = normParams.param2;
patRec.algorithm.model = model;
patRec.performance = performance;
disp("Training and testing finished")

end