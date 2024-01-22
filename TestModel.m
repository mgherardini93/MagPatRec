function [patRec, message] = TestModel(patRec, dataSet, dataOut)
tSet = dataSet;
tOut = dataOut;

nM = patRec.nM;
Nt = size(tSet,1);
performance = [];

% Normalize the test data -------------------------------------------------
switch patRec.normalization.type
    case 'Midrange0Range2'
        tSet = (tSet - patRec.normalization.param1)./(patRec.normalization.param2/2);
    case 'Mean0Std1'
        tSet = (tSet - patRec.normalization.param1)./patRec.normalization.param2;
end

% Feature selection -------------------------------------------------------
switch patRec.featSelection.algorithm  
    case 'PCA' %-----------------------------------------------------------
        tSet = tSet*patRec.featSelection.coeff;
    otherwise
end

if ~isfield(patRec.algorithm,'model')
    message = 'Error: classifier not trained!';
    return;
end

% Test the classifier -----------------------------------------------------
model = patRec.algorithm.model;
score = zeros(Nt,nM);
for m = 1:nM
    [~, f] = predict(model{m},tSet);
    score(:,m) = f(:,2);
end
[~,prediction] = max(score,[],2);
[label,~] = find(tOut');

confMat = confusionmat(label,prediction);
if (patRec.plotConfMat)
    figure
    cm = confusionchart(confMat,patRec.mov);
    sortClasses(cm,patRec.mov)
end
performance.confMat = confMat;
performance.accuracy = trace(confMat)/sum(sum(confMat)).*100;
performance.classAccuracy = diag(confMat)./sum(confMat,2);
performance.precision = diag(confMat)./sum(confMat,1)';
performance.f1 = (2.*performance.precision.*performance.classAccuracy)./(performance.precision + performance.classAccuracy);
performance.classAccuracy = performance.classAccuracy.*100;
performance.precision = performance.precision.*100;
performance.prediction = prediction;

patRec.performance = performance;

message = 'Correctly tested.';
disp("Testing is finished")

end