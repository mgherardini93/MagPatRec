function [trSetFinal, tSetFinal, message] = FeatureSelection_Correlation(patRec,trSet,trOut, tSet)

trSetFinal = [];
tSetFinal = [];
message = 'Features Selected!';
[labels,~] = find(trOut');
rhoL = corr(trSet,labels); % calculating the correlation with the class
rhoL = abs(rhoL);
nPairs = size(patRec.magnetPairs,1);
nCol = nPairs*sum(patRec.magSelFeatures);
nFeat = length(patRec.selFeatures)*sum(patRec.magSelFeatures);
rhoL_feat = zeros(nFeat,nPairs);

% Averaging the absolute correlation across magnet pairs
for i = 1:nPairs
    if (nCol ~= nPairs)
        rhoL_feat(1:2:end,i) = rhoL((1:nCol:end)+2*(i-1));
        rhoL_feat(2:2:end,i) = rhoL((2:nCol:end)+2*(i-1));
    else
        rhoL_feat(:,i) = rhoL((1:nCol:end)+(i-1));
    end
end
rhoL_avg = mean(rhoL_feat,2); % contains the averaged absolute correlation coef. for each feature

featPassed = rhoL_avg > patRec.featSelection.scoreLabel; % checking which features passed the threshold

if (nCol ~= nPairs) % in case both linear and angular distances are used --
    l_featPassed = featPassed(1:2:end);
    a_featPassed = featPassed(2:2:end);

    if (patRec.plotFigures) % plotting the correlation with the class
        figure('Position',[600 300 250 500])
        yValues = cell(1,nFeat);
        yValues(1:2:end) = strcat(patRec.selFeatures,'_l');
        yValues(2:2:end) = strcat(patRec.selFeatures,'_a');
        heatmap({'label'},yValues,rhoL_avg)
        title('Correlation of the features with label')
    end
%     disp('Linear features that are correlated with label: '); disp(patRec.selFeatures(l_featPassed)');
%     disp('Angular features that are correlated with label: '); disp(patRec.selFeatures(a_featPassed)');

    linFeat_col = 1:2:length(rhoL); angFeat_col = 2:2:length(rhoL);
    featPassed_col = [];
    % extracting the column numbers of linear distance features which passed the threshold
    for i = 1:length(l_featPassed)
        if(l_featPassed(i))
            featPassed_col = [featPassed_col, linFeat_col((1:nPairs)+(nPairs)*(i-1))];
        end
    end
    trSet_lin = trSet(:,featPassed_col);
    tSet_lin = tSet(:,featPassed_col);

    % extracting the column numbers of angular distance features which passed the threshold
    featPassed_col= [];
    for i = 1:length(a_featPassed)
        if(a_featPassed(i))
            featPassed_col = [featPassed_col, angFeat_col((1:nPairs)+(nPairs)*(i-1))];
        end
    end
    trSet_ang = trSet(:,featPassed_col);
    tSet_ang = tSet(:,featPassed_col);

    nFeat_new = sum(l_featPassed)+sum(a_featPassed);
    if (~nFeat_new) % no features passed the threshold
        message = 'Error: Label correlation coefficients are too high!';
        return;
    end

    % calculating the feature correlation matrix for each magnet pair
    corrMat = zeros(nFeat_new,nFeat_new,nPairs);
    for i = 1:nPairs
        tempSet = [trSet_lin(:,(1:nPairs:end)+i-1), trSet_ang(:,(1:nPairs:end)+i-1)];
        corrMat(:,:,i) = abs(corr(tempSet));
    end
    corrMat = mean(corrMat,3); % average the correlation matrix across magnet pairs
    corrMat = tril(corrMat,-1); % make the matrix triangular 

    if (patRec.plotFigures) % plot the correlation matrix
        figure()
        tempCorrMat = corrMat; 
        tempCorrMat(tempCorrMat == 0) = NaN;
        heatmap([strcat(patRec.selFeatures(l_featPassed),'_l'),...
            strcat(patRec.selFeatures(a_featPassed),'_a')],...
            [strcat(patRec.selFeatures(l_featPassed),'_l'),...
            strcat(patRec.selFeatures(a_featPassed),'_a')], tempCorrMat);
        title('Correlation between features')
    end

    scoresFeat = sum(corrMat > patRec.featSelection.scoreFeat,2); % check how many features have the correlation above the threshold
    trSetFinal = [trSet_lin, trSet_ang];
    tSetFinal = [tSet_lin, tSet_ang];
    % features which are not correlated with others will have 0 in the scoresFeat variable
    featPassed_col = (1:nPairs) + (find(scoresFeat==0)-1)*nPairs; % select features which are not correlated with other features
    if isempty(featPassed_col)
        message = 'Error: Feature correlation coefficient are too low!';
        return;
    end
    featPassed_col = reshape(featPassed_col,1,[]);
    trSetFinal = trSetFinal(:,featPassed_col);
    tSetFinal = tSetFinal(:,featPassed_col);

    disp('Final features: ');
    disp('Linear: '); temp = scoresFeat(1:sum(l_featPassed)); disp(patRec.selFeatures(temp==0)); 
    patRec.featSelection.linFeatures = patRec.selFeatures(temp==0);
    disp('Angular: '); temp = scoresFeat(sum(l_featPassed)+1:end); disp(patRec.selFeatures(temp==0));  
    patRec.featSelection.angFeatures = patRec.selFeatures(temp==0);

else % in case just linear or angular distance is used --------------------
    if (patRec.plotFigures) % plotting the correlation with the class
        figure('Position',[600 300 250 500])
        yValues = patRec.selFeatures;
        heatmap({'label'},yValues,rhoL_avg)
        title('Correlation of the features with label')
    end
    
    featPassed_col = [];
    cols = 1:length(rhoL);
    % extracting the column numbers of features that passed the threshold
    for i = 1:length(featPassed)
        if(featPassed(i))
            featPassed_col = [featPassed_col, cols((1:nPairs)+(nPairs)*(i-1))];
        end
    end
    trSetFinal = trSet(:,featPassed_col);
    tSetFinal = tSet(:,featPassed_col);

    nFeat_new = sum(featPassed);
    if (~nFeat_new)
        message = 'Error: Label correlation coefficients are to high!';
        return;
    end

    % calculating the feature matrix for each magnet pair
    corrMat = zeros(nFeat_new,nFeat_new,nPairs);
    for i = 1:nPairs
        tempSet = trSetFinal(:,(1:nPairs:end)+i-1);
        corrMat(:,:,i) = abs(corr(tempSet));
    end
    corrMat = mean(corrMat,3); % average the correlation matrix across magnet pairs
    corrMat = tril(corrMat,-1); % make the matrix triangular

    if (patRec.plotFigures) % plot the correlation matrix
        figure()
        tempCorrMat = corrMat; 
        tempCorrMat(tempCorrMat == 0) = NaN;
        heatmap(patRec.selFeatures(featPassed), patRec.selFeatures(featPassed), tempCorrMat);
        title('Correlation between features')
    end

    scoresFeat = sum(corrMat > patRec.featSelection.scoreFeat,2);
    % features which are not correlated with others will have 0 in the scoresFeat variable
    featPassed_col = (1:nPairs) + (find(scoresFeat==0)-1)*nPairs;
    featPassed_col = reshape(featPassed_col',1,[]);
    trSetFinal = trSetFinal(:,featPassed_col);
    tSetFinal = tSetFinal(:,featPassed_col);
    disp('Final features: '); disp(patRec.selFeatures(scoresFeat==0));
    patRec.featSelection.selFeatures = patRec.selFeatures(scoresFeat==0);
end

end

