function [patRec, message] = TrainModel(patRec, dataSet, dataOut)

trSet = dataSet;
trOut = dataOut;
message = 'The model is trained!';

% Onset detector ----------------------------------------------------------
if (patRec.algorithm.transient)
    if (strcmp(patRec.mov,'Rest'))
        message = 'Error: Remove the "Rest" class!';
        return;
    end
    [trSet, trOut, threshold, normParams, ~, message] = OnsetDetectorTraining(patRec, trSet, trOut,patRec.nR);
    if (contains(message,'Error'))
        return;
    end
    patRec.normalization.param1 = normParams.param1;
    patRec.normalization.param2 = normParams.param2;
    patRec.algoritm.threshold = threshold;
end

nM = patRec.nM;
Ntr = size(trSet,1);

if (patRec.algorithm.transient) % In case of transient, the data is normalizaed within the ODA
    normParams.param1 = patRec.normalization.param1;
    normParams.param2 = patRec.normalization.param2;
else                            % In case of steady state, the data is normalized here
% Normalize the training data ---------------------------------------------
    [trSet, normParams] = NormalizeData(trSet, patRec.normalization.type);
    if (isempty(normParams.param1) || isempty(normParams.param2))
        message = 'Error: within normalization';
        return;
    end  
end

% Feature selection -------------------------------------------------------
switch patRec.featSelection.algorithm
    case 'MRMR' %----------------------------------------------------------
        [labelTr,~] = find(trOut');
        [idxF,scores] = fscmrmr(trSet,labelTr);
        idx = 0;
        sumPrc = 0;
        sumScore = sum(scores);
        while sumPrc < patRec.featSelection.scorePercentage
            idx = idx + 1;
            sumPrc = sumPrc + scores(idxF(idx))/sumScore*100;
        end
        trSet = trSet(:,idxF(1:idx-1));
        disp(['Number of features selected: ', num2str(idx-1)]);

    case 'Chi-square' %----------------------------------------------------
        [labelTr,~] = find(trOut');
        [idxF,scores] = fscchi2(trSet,labelTr);
        idx = 0;
        sumPrc = 0;
        idxInf = isinf(scores);
        if isempty(idxInf)
            sumScore = sum(scores);
            while sumPrc < patRec.featSelection.scorePercentage
                idx = idx + 1;
                sumPrc = sumPrc + scores(idxF(idx))/sumScore*100;
            end
            trSet = trSet(:,idxF(1:idx-1));
            disp(['Number of features selected: ', num2str(idx-1)]);
        else
            trSet = trSet(:,idxInf);
            disp(['Scores were infinite. Number of features selected: ', num2str(sum(idxInf))]);
        end
        
    case 'PCA' %-----------------------------------------------------------
        [coeff,~,~,~,scorePrct] = pca(trSet);
        score = 0;
        idx = 0;
        while score < patRec.featSelection.scorePercentage
            idx = idx + 1;
            score = score + scorePrct(idx);
        end
        trSet = trSet*coeff(:,1:(idx-1));
        disp(['First ' num2str(idx-1) ' components of PCA are taken']);
        PCAcoeff = coeff(:,1:(idx-1));

    case 'Correlation'
        [trSet, ~, message] = FeatureSelection_Correlation(patRec,trSet,trOut,trSet);
        if contains(message,'Error')
            return;
        end
    otherwise
end

% Train the classifier ----------------------------------------------------
%%%%%%%%%%%% TO DO: add validation part %%%%%%%%%%%
model = cell(1,nM);
switch patRec.algorithm.type
    case 'SVM'
        for m = 1:nM
            trGroup = zeros(Ntr,1);
            trGroup(trOut(:,m) == 1) = 1;
            switch patRec.algorithm.kernel
                case 'Linear'
                    model{m} = fitcsvm(trSet,trGroup,'KernelFunction','linear','BoxConstraint',1,'Standardize',false);
                case 'Polynomial'
                    model{m} = fitcsvm(trSet,trGroup,'KernelFunction','polynomial','BoxConstraint',1,'PolynomialOrder',2,'Standardize',false);
                case 'Gaussian'
                    model{m} = fitcsvm(trSet,trGroup,'KernelFunction','rbf','BoxConstraint',5,'KernelScale',1,'Standardize',false);
            end
        end
    otherwise
end

if isempty(model)
    message = 'Error: classifier not trained!';
    return;
end

if (contains(message,'Error'))
    return;
end

patRec.normalization.param1 = normParams.param1;
patRec.normalization.param2 = normParams.param2;
patRec.algorithm.model = model;
if (exist('PCAcoeff')==1)
    patRec.featSelection.coeff = PCAcoeff;
end
disp("Training is finished")

end