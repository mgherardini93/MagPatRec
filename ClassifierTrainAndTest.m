function [normParams, model, performance, message] = ClassifierTrainAndTest(patRec,trSet,trOut,tSet,tOut)

nM = patRec.nM;
Ntr = size(trSet,1);
Nt = size(tSet,1);
model = [];
performance = [];

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
    
% Normalize the test data -------------------------------------------------
    switch patRec.normalization.type
        case 'Midrange0Range2'
            tSet = (tSet - normParams.param1)./(normParams.param2/2);
        case 'Mean0Std1'
            tSet = (tSet - normParams.param1)./normParams.param2;
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
        tSet = tSet(:,idxF(1:idx-1));
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
            tSet = tSet(:,idxF(1:idx-1));
            disp(['Number of features selected: ', num2str(idx-1)]);
        else
            trSet = trSet(:,idxInf);
            tSet = tSet(:,idxInf);
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
        tSet = tSet*coeff(:,1:(idx-1));
        disp(['First ' num2str(idx-1) ' components of PCA are taken']);
        performance.coeff = coeff(:,1:(idx-1));

    case 'Correlation' %---------------------------------------------------
        [trSet, tSet, message] = FeatureSelection_Correlation(patRec,trSet,trOut,tSet);
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
                    model{m} = fitcsvm(trSet,trGroup,'KernelFunction','polynomial','BoxConstraint',1,'PolynomialOrder',3,'Standardize',false);
                case 'Gaussian'
                    model{m} = fitcsvm(trSet,trGroup,'KernelFunction','rbf','BoxConstraint',5,'KernelScale',1,'Standardize',false);
            end
        end
    case 'MLR'
        [model{1}.coeff, ~, model{1}.residuals, model{1}.CovB] = mvregress([ones(size(trSet,1),1), trSet],trOut);
    otherwise
end

if isempty(model)
    message = 'Error: classifier not trained!';
    return;
end

% Test the classifier -----------------------------------------------------
if (~strcmp(patRec.algorithm.type,'MLR'))
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
else
% Test the regressor ------------------------------------------------------
    if (patRec.plotFigures)
        figure()
        for dof = 1:nM
            subplot(nM,1,dof)
            hold on
            plot(tOut(:,dof),'b');
            plot([ones(size(tSet,1),1), tSet]*model{1}.coeff(:,dof),'r')
            title("DOF" + num2str(dof))
            hold off
        end
    end
    performance.accuracy = 0;
    performance.classAccuracy = 0;
    performance.precision = 0;
    performance.f1 = 0;
    performance.classAccuracy = 0;
    performance.precision = 0;
    performance.prediction = 0;
end

message = 'Correctly trained and tested.';

end