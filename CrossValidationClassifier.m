function [patRec, message] = CrossValidationClassifier(patRec, dataSet, dataOut)

dataLength = size(dataSet,1);
nM = patRec.nM;
nR = patRec.nR;
nRTest = patRec.testingRep;
nFolds = fix(nR/nRTest);
nCh = size(dataSet,2);      % number of channels and features
Nt = round(dataLength/nFolds); % number of samples for training
Ntr = dataLength - Nt;      % number of samples for testing
if (mod(nR, nRTest))
    disp('Warning: The number of folds is not a round number!')
    disp("The number of repetitions that is left out is: " + num2str(nR - nFolds*nRTest));
end

% restLength = sum(dataOut(:,end) == 1);
% ind = randperm(restLength);
% restSet = dataSet(dataOut(:,end) == 1,:);
% restSet = restSet(ind,:);
% dataSet(dataOut(:,end) == 1,:) = restSet;

% Going fold by fold ------------------------------------------------------
accuracy = 0;
classAccuracy = zeros(nM,1);
precision = zeros(nM,1);
f1 = zeros(nM,1);
confMat = zeros(nM,nM);
for f = 1:nFolds
    disp("Fold " + num2str(f));
    idxD = 0;  % index to follow the dataSet
    idxTr = 0; % index to follow trSet
    idxT = 0;  % index to follow tSet
    trSet = []; trOut =[];
    tSet =[]; tOut = [];
    % Separate data for testing and training ------------------------------
    for m = 1:nM
        movementLength = sum(dataOut(:,m) == 1); % find which data correspond to the movement
        foldLength = fix(movementLength/nFolds);
        tIndex = (1:foldLength) + (f-1)*foldLength;  % index of the data that goes into test
        trIndex = setdiff(1:movementLength,tIndex);  % index of the data that goes into training    
        trSet((1:length(trIndex))+idxTr,:) = dataSet(trIndex+idxD,:);
        tSet((1:length(tIndex))+idxT,:) = dataSet(tIndex+idxD,:);
        trOut((1:length(trIndex))+idxTr,m) = 1;
        tOut((1:length(tIndex))+idxT,m) = 1;  
        % increase indeces to reach the next movement
        idxD = idxD + movementLength; 
        idxTr = idxTr + length(trIndex);     
        idxT = idxT + length(tIndex); 
    end

    % Onset detector ------------------------------------------------------
    if (patRec.algorithm.transient)
        if (strcmp(patRec.mov,'Rest'))
            message = 'Error: Remove the "Rest" class!';
            return;
        end
        [trSet, trOut, threshold, normParams, dNoiseMAV, message] = OnsetDetectorTraining(patRec, trSet, trOut,nR - nRTest);
        if (contains(message,'Error'))
            return;
        end
        [tSet, tOut] = OnsetDetectorTest(patRec, tSet, tOut, threshold, normParams, dNoiseMAV, nRTest);
        patRec.normalization.param1 = normParams.param1;
        patRec.normalization.param2 = normParams.param2;
    end

    % Train and test the classifier ---------------------------------------
    [normParams, model, performanceFold, message] = ClassifierTrainAndTest(patRec,trSet,trOut,tSet,tOut);
    
    if (contains(message,'Error'))
        return;
    end

    performance.fold{f} = performanceFold;
    accuracy = accuracy + performanceFold.accuracy;
    classAccuracy = classAccuracy + performanceFold.classAccuracy;
    precision = precision + performanceFold.precision;
    f1 = f1 + performanceFold.f1;
    confMat = confMat + performanceFold.confMat;
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