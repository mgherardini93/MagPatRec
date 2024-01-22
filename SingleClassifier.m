function [patRec, message] = SingleClassifier(patRec, dataSet, dataOut)

tP = patRec.trainingPercentage;
dataLength = size(dataSet,1);
nCh = size(dataSet,2);      % number of channels and features
Ntr = round(tP*dataLength); % number of samples for training
Nt = dataLength - Ntr;      % number of samples for testing
nM = patRec.nM;

trSet = zeros(Ntr,nCh); % training set
tSet = zeros(Nt,nCh);   % test set
trOut = zeros(Ntr,nM);  % training out
tOut = zeros(Nt,nM);    % test out
idxD = 0;  % index to follow the dataSet
idxTr = 0; % index to follow trSet
idxT = 0;  % index to follow tSet
notAClassLabel = 0;
if strcmp(patRec.algorithm.type,'MLR')
    notAClassLabel = -1;
end

% Separate training form test data ----------------------------------------
for m = 1:nM
    movementLength = sum(dataOut(:,m) ~= notAClassLabel); % find which data correspond to the movement
    trLength = round(movementLength*tP);     % length of the data that goes into training 
    tLength = movementLength - trLength;     % length of the data that goes into test
    trSet((1:trLength)+idxTr,:) = dataSet((1:trLength)+idxD,:);
    tSet((1:tLength)+idxT,:) = dataSet((1:tLength)+idxD+trLength,:);
    trOut((1:trLength)+idxTr,m) = dataOut((1:trLength)+idxD,m);
    tOut((1:tLength)+idxT,m) = dataOut((1:tLength)+idxD+trLength,m); 
    % increase indeces to reach the next movement
    idxD = idxD + movementLength; 
    idxTr = idxTr + trLength;     
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

% Onset detector ----------------------------------------------------------
if (patRec.algorithm.transient)
    if (strcmp(patRec.mov,'Rest'))
        message = 'Error: Remove the "Rest" class!';
        return;
    end
    nR = patRec.nR*patRec.trainingPercentage;
    if mod(nR,1)
        warning('The traning percentage does not give full number of repetitions!');
        nR = round(nR);
    end
    [trSet, trOut, threshold, normParams, dNoiseMAV, message] = OnsetDetectorTraining(patRec, trSet, trOut,nR);
    if (contains(message,'Error'))
        return;
    end
    nR = patRec.nR - nR;
    [tSet, tOut] = OnsetDetectorTest(patRec, tSet, tOut, threshold, normParams, dNoiseMAV, nR);
    patRec.normalization.param1 = normParams.param1;
    patRec.normalization.param2 = normParams.param2;
end

% Train and test the classifier -------------------------------------------
[normParams, model, performance, message] = ClassifierTrainAndTest(patRec,trSet,trOut,tSet,tOut);

if (contains(message,'Error'))
    return;
end

patRec.normalization.param1 = normParams.param1;
patRec.normalization.param2 = normParams.param2;
patRec.algorithm.model = model;
if isfield(performance,'coeff')
    patRec.featSelection.coeff = performance.coeff;
    performance = rmfield(performance,'coeff');
end
patRec.performance = performance;
disp("Training and testing finished")

end