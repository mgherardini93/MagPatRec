function [patRec, message] = CrossValidationClassifierTestFull(patRec, dataSet, dataOut, dataStruct)

movUsed = patRec.mov;
if(sum(strcmp(movUsed,'Rest')))
    movUsed(strcmp(movUsed,'Rest')) = [];
end
featuresStruct = dataStruct.exFeaturesTrans(ismember(dataStruct.mov,movUsed)); % struct which contains full data
% dataLength = size(dataSet,1);
nM = patRec.nM;
nMt = length(movUsed); % test data doesn't have a separate collumn for rest class
nR = patRec.nR;
nRTest = patRec.testingRep;
nFolds = fix(nR/nRTest);
% nCh = size(dataSet,2);      % number of channels and features
% Nt = round(dataLength/nFolds); % number of samples for training
% Ntr = dataLength - Nt;      % number of samples for testing
if (mod(nR, nRTest))
    disp('Warning: The number of folds is not a round number!')
    disp("The number of repetitions that is left out is: " + num2str(nR - nFolds*nRTest));
end

% Stack selected features in a matrix for test ----------------------------
nCh = patRec.nCh;
chIdx = patRec.vCh;
testDataLength = 0;
for m = 1:nMt
    testDataLength = testDataLength + length(featuresStruct{m});
end
testDataSet = zeros(testDataLength,nCh*length(patRec.selFeatures));
testDataOut = zeros(testDataLength,nM);
dSIdx = 0; % index that indicated at which point to enter data in dataSet
for m = 1:nMt
    movementLength = length(featuresStruct{m});
    % Getting test data set
    for s = 1:movementLength
        for f = 1:length(patRec.selFeatures)
            testDataSet(dSIdx+s,(1:nCh) + (f-1)*nCh) = featuresStruct{m}(s).(patRec.selFeatures{f})(chIdx); 
        end
    end
    % Getting test labels
    label = dataStruct.label{m};
    labelLength = length(label);
    if (labelLength > movementLength) % labels need to be downsamples to correspond to the windowed data
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
        testDataOut(contractionIdx + dSIdx,m) = 1; % label the contractions
        testDataOut(restIdx + dSIdx, nM) = 1;      % label the rest part (rest class is always the last collumn)
    end
%     testDataOut((dSIdx+1):(dSIdx + movementLength),m) = 1;
    dSIdx = dSIdx + movementLength;
end
% Nt = round((1-tP)*testDataLength); % number of samples for testing
% tSet = zeros(Nt,nCh);              % test set
% tOut = zeros(Nt,nM);               % test out

% Going fold by fold ------------------------------------------------------
accuracy = 0;
classAccuracy = zeros(nM,1);
precision = zeros(nM,1);
f1 = zeros(nM,1);
confMat = zeros(nM,nM);
label = [];
prediction = [];
for f = 1:nFolds
    disp("Fold " + num2str(f));
    idxD = 0;  % index to follow the dataSet
    idxDT = 0;  % index to follow the testDataSet
    idxTr = 0; % index to follow trSet
    idxT = 0;  % index to follow tSet
    trSet = []; trOut =[];
    tSet =[]; tOut = [];
    % Separate data for training ------------------------------------------
    for m = 1:nM
        movementLength = sum(dataOut(:,m) == 1); % find which data correspond to the movement
        foldLength = fix(movementLength/nFolds);
        tIndex = (1:foldLength) + (f-1)*foldLength;  % index of the data that goes into test
        trIndex = setdiff(1:movementLength,tIndex);  % index of the data that goes into training    
        trSet((1:length(trIndex))+idxTr,:) = dataSet(trIndex+idxD,:);
        trOut((1:length(trIndex))+idxTr,m) = 1; 
        % increase indeces to reach the next movement
        idxD = idxD + movementLength; 
        idxTr = idxTr + length(trIndex);     
    end
    % Separate data for testing -------------------------------------------
    for m = 1:nMt
        movementLength = length(featuresStruct{m}); % find which data correspond to the movement
        foldLength = fix(movementLength/nFolds);
        tIndex = (1:foldLength) + (f-1)*foldLength;     % length of the data that goes into test
        tSet((1:length(tIndex))+idxT,:) = testDataSet(tIndex+idxDT,:);
        tOut((1:length(tIndex))+idxT,:) = testDataOut(tIndex+idxDT,:);  
        % increase indeces to reach the next movement
        idxDT = idxDT + movementLength;   
        idxT = idxT + length(tIndex); 
    end

    % Onset detector ------------------------------------------------------
%     if (patRec.algorithm.transient)
%         if (strcmp(patRec.mov,'Rest'))
%             message = 'Error: Remove the "Rest" class!';
%             return;
%         end
%         [trSet, trOut, threshold, normParams, dNoiseMAV, message] = OnsetDetectorTraining(patRec, trSet, trOut,nR - nRTest);
%         if (contains(message,'Error'))
%             return;
%         end
%         [tSet, tOut] = OnsetDetectorTest(patRec, tSet, tOut, threshold, normParams, dNoiseMAV, nRTest);
%         patRec.normalization.param1 = normParams.param1;
%         patRec.normalization.param2 = normParams.param2;
%     end

    % Train and test the classifier ---------------------------------------
    [normParams, model, performanceFold, message] = ClassifierTrainAndTest(patRec,trSet,trOut,tSet,tOut);
    
    if (contains(message,'Error'))
        return;
    end

%     if (patRec.plotFigures)
%         figure('Name',"Continuous performance of the classifier for fold " + num2str(f));
%         plot(-performanceFold.prediction,'b-')
%         hold on
%         [label,~] = find(tOut');
%         plot(-label,'r-')
%         legend('Classifier predictions','Ground truth');
%         hold off
%         xlabel('Samples')
%         ylabel('Classes')
%         yticklabels(patRec.mov(end:-1:1));
%     end
    prediction = [prediction; performanceFold.prediction];
    [temp,~] = find(tOut');
    label = [label; temp];

    performance.fold{f} = performanceFold;
    accuracy = accuracy + performanceFold.accuracy;
    classAccuracy = classAccuracy + performanceFold.classAccuracy;
    precision = precision + performanceFold.precision;
    f1 = f1 + performanceFold.f1;
    confMat = confMat + performanceFold.confMat;
end

if (patRec.plotFigures)
    figure('Name',"Continuous performance of the classifier");
    plot(-prediction,'b-')
    hold on
    plot(-label,'r-')
    legend('Classifier predictions','Ground truth');
    hold off
    xlabel('Samples')
    ylabel('Classes')
    yticks(-patRec.nM:-1)
    yticklabels(patRec.mov(end:-1:1));
end

if (patRec.plotConfMat)
    figure
    cm = confusionchart(confMat,patRec.mov);
    sortClasses(cm,patRec.mov)
end

performance.accuracy = accuracy/nFolds;
performance.classAccuracy = classAccuracy/nFolds;
performance.precision = precision/nFolds;
performance.f1 = f1/nFolds;
performance.confMat = confMat;
patRec.performance = performance;
disp("Training and testing finished")

end