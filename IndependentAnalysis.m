clear all
clc
%% Variables
nMag = 6;  % Number of magnets
folderName = "C:\Users\UTENTE\OneDrive - Scuola Superiore Sant'Anna\Desktop\Work\MYKI\Experiment data\";
day = "03"; month = "05";
folderDate = month + day + "\"; % date of the data aquisition
dataDate = day + month;
dataName = {"data_77Mov07" + dataDate + ".txt", "data_78Mov08" + dataDate + ".txt"...
    "data_79Mov10" + dataDate + ".txt", "data_80Mov09" + dataDate + ".txt"...
    "data_81Mov01" + dataDate + ".txt", "data_82Mov05" + dataDate + ".txt", "data_83Mov06" + dataDate + ".txt"};
% ClassNames = {'Cylindrical', 'Lateral','Pinch','Open hand','Wrist flexion','Wrist extension','Pronation','Supination','Radial deviation','Ulnar Deviation','Rest'};
ClassNames = {'Io', 'Middle Finger','Flex all fingers','Open hand','Wrist flexion','Wrist extension','Pronation','Supination','Radial deviation','Ulnar Deviation','Thumb flexion','Tu','Rest'};
% cT = 5; % contraction time
% rT = 5; % rest time
% nR = 15; % number of repetitions per movement
% orderOfMov = [1 2 3 4 5 6 7 8 9 10]; % movement files that need to be read
orderOfMov = [7,8,12,1,3,2,6];
nMov = length(dataName);  % Number of classes

%% Reading .txt files and saving .mat
% restData = load(['Rest', dataDate, '.txt']);
% restData = readmatrix(['Rest', dataDate, '.txt'], 'DecimalSeparator','.');
% restData = restData(:,1:(6*nMag));
% restData = table2array(restData);
data =[]; labels = [];
for m = 1:nMov
    temp = readmatrix(folderName + folderDate + dataName{m},'DecimalSeparator','.');
%     temp = readmatrix(['Mov0', num2str(orderOfMov(m)), dataDate, '.txt'], 'DecimalSeparator',',');
%     temp = table2array(temp);
    labels = [labels; ones(size(temp,1),1).*m, temp(:,37)];
    data = [data; temp(:,1:(6*nMag))];
end

%% Abostolute magnet distance
magPairs = nchoosek(1:nMag,2);
idxCoo = (1:3)' + (0:(nMag-1))*6;
idxCoo = reshape(idxCoo,1,[]);
absDist = cell(1,nMov);
for classToCheck = 1:nMov
    xyz = data(labels(:,1) == classToCheck,idxCoo);
    figure('Name',['Class', num2str(classToCheck)])
    for i = 1:nMag
        iMP = (1:3) + 3*(i-1); % get the magnet pair
        mP1 = xyz(:,iMP);
        dist = mP1(1,1:3) - mP1(:,1:3);
%         distN = (d_rest(i) - vecnorm(dist,2,2)); % calculate linear distance
        distN = vecnorm(dist,2,2); % calculate linear distance
        absDist{classToCheck} = [absDist{classToCheck}, distN];
        subplot(3,2,i)
        plot(distN*1000)
%         hold on
%         plot(repmat((0:cTsamples:l)',1,2)',repmat(ylim,length(0:cTsamples:l),1)','r--')
%         hold off
%         disp(max(distN)*1000)
%         temp_max = [temp_max, max(distN)*1000];
        title("Magnet " + num2str(i));
        ylabel('Distance (m)')
    end
%     displ_max = [displ_max; temp_max];
%     temp_max = [];
end

%%
nR = 3;
offsetR = 0;
% dataStruct = struct([]);
dataStruct.("data" + dataDate).rest = cell(1,nMov);
dataStruct.("data" + dataDate).contr = cell(1,nMov);
dataStruct.("data" + dataDate).classes = ClassNames(orderOfMov);
for m = 1:nMov
    labelMov = labels(labels(:,1) == m,2);
    dL = diff(labelMov);
    restStop = find(dL>0);
    contrStop = find(dL<0);
    restData = [];
    contrData = [];
    if (restStop(1)<contrStop(1))
        contrStop = [0; contrStop];
        for r = (1:nR) + offsetR
            temp = absDist{m}(contrStop(r)+1:restStop(r),:)*1000;
            l = length(temp);
            temp1 = temp(round(l/4):round(l/4)+round(l/2),:);
            restData = [restData; temp1];
            temp = absDist{m}(restStop(r)+1:contrStop(r+1),:)*1000;
            l = length(temp);
            temp1 = temp(round(l/4):round(l/4)+round(l/2),:);
            contrData = [contrData; temp1];
        end
    else
        restStop = [0; restStop];
        for r = (1:nR) + offsetR
            temp = absDist{m}(contrStop(r)+1:restStop(r+1),:)*1000;
            l = length(temp);
            temp1 = temp(round(l/4):round(l/4)+round(l/2),:);
            restData = [restData; temp1];
            temp = absDist{m}(restStop(r)+1:contrStop(r),:)*1000;
            l = length(temp);
            temp1 = temp(round(l/4):round(l/4)+round(l/2),:);
            contrData = [contrData; temp1];
        end
    end
    dataStruct.("data" + dataDate).rest{m} = restData;
    dataStruct.("data" + dataDate).contr{m} = contrData;
end

%% Check idividual movement signal (linear distance is plotted)
% magPairs = nchoosek(1:nMag,2);
% idxCoo = (1:3)' + (0:(nMag-1))*6;
% idxCoo = reshape(idxCoo,1,[]);
% % cTsamples = cT*sF;
% % rTsamples = rT*sF;
% if (size(magPairs,1) > 7)
%     subPlotcollums = 2;
%     subPlotrows = fix(size(magPairs,1)/2) + 1;
% else
%     subPlotcollums = 1;
%     subPlotrows = size(magPairs,1);
% end
% % displ_max = [];
% % temp_max = [];
% for classToCheck = 1:nMov
%     xyz = data(labels(:,1) == classToCheck,idxCoo);
%     figure('Name',['Class', num2str(classToCheck)])
%     for i = 1:size(magPairs,1)
%         iMP = [(1:3) + 3*(magPairs(i,2)-1), (1:3) + 3*(magPairs(i,1)-1)]; % get the magnet pair
%         mP1 = xyz(:,iMP);
%         dist = mP1(:,1:3) - mP1(:,4:6);
% %         distN = (d_rest(i) - vecnorm(dist,2,2)); % calculate linear distance
%         distN = vecnorm(dist,2,2); % calculate linear distance
%         subplot(subPlotrows,subPlotcollums,i)
%         plot(distN*1000)
% %         hold on
% %         plot(repmat((0:cTsamples:l)',1,2)',repmat(ylim,length(0:cTsamples:l),1)','r--')
% %         hold off
% %         disp(max(distN)*1000)
% %         temp_max = [temp_max, max(distN)*1000];
%         title("Magnet pair " + num2str(i));
%         ylabel('Distance (m)')
%     end
% %     displ_max = [displ_max; temp_max];
% %     temp_max = [];
% end

%% Analyse the data struct
dataFieldNames = fieldnames(dataStruct);

for d = 1:length(dataFieldNames)
    meanData1{d} = abs(dataStruct.(dataFieldNames{d}).contr{1} - mean(dataStruct.(dataFieldNames{d}).rest{1}));
    meanData2{d} = abs(dataStruct.(dataFieldNames{d}).contr{2} - mean(dataStruct.(dataFieldNames{d}).rest{2}));
    mD1(d,:) = mean(meanData1{d});
    sD1(d,:) = std(meanData1{d});
    mD2(d,:) = mean(meanData2{d});
    sD2(d,:) = std(meanData2{d});
    classes1{d} = dataStruct.(dataFieldNames{d}).classes{1};
    classes2{d} = dataStruct.(dataFieldNames{d}).classes{2};
end

%%
lineColor = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880];
titles = {'ED_p','ED_d','FPL_p','FPL_d','FCU_p','FCU_d'};
xValues = 1:length(dataFieldNames);
figure
for m = 1:nMag
    subplot(3,2,m)
    plot(xValues,mD1(:,m),'Color',lineColor(m,:))
    hold on
    patch([xValues fliplr(xValues)], [mD1(:,m)'-sD1(:,m)' fliplr(mD1(:,m)'+sD1(:,m)')], lineColor(m,:),'FaceAlpha',0.3,'EdgeColor','none')
    hold off
    xticks(xValues);
    xticklabels(classes1);
    title(titles{m});
    ylabel('Displacement (mm)')
    ylim([0 8])
end

figure
for m = 1:nMag
    subplot(3,2,m)
    plot(xValues,mD2(:,m),'Color',lineColor(m,:))
    hold on
    patch([xValues fliplr(xValues)], [mD2(:,m)'-sD2(:,m)' fliplr(mD2(:,m)'+sD2(:,m)')],lineColor(m,:),'FaceAlpha',0.3,'EdgeColor','none')
    hold off
    xticks(xValues);
    xticklabels(classes2);
    title(titles{m});
    ylabel('Displacement (mm)')
    ylim([0 8])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse the extended data struct
nMag = 6;
dataFieldNames = fieldnames(dataStruct);
% dataFieldNames(5) = [];
mD = []; sD =[];
for d = 1:length(dataFieldNames)
    for m = 1:length(dataStruct.(dataFieldNames{d}).classes)
        meanData{d,m} = abs(dataStruct.(dataFieldNames{d}).contr{m} - mean(dataStruct.(dataFieldNames{d}).rest{m}));
        mD(m,:,d) = mean(meanData{d,m});
        sD(m,:,d) = std(meanData{d,m});
    end
end

% for d = 1:length(dataFieldNames)
%     for m = 1:length(dataStruct.(dataFieldNames{d}).classes)
%         classR = dataStruct.(dataFieldNames{d}).classes{m};
%         if (strcmp(classR,'Supination')) || (strcmp(classR,'Pronation'))
%             mD(m,:,d) = 0;
%             sD(m,:,d) = 0;
%         end
%     end
% end

[maxValue, idxMax] = max(mD);
stdMax =  [];
classMax = {};
for mm = 1:nMag
    for d = 1:length(dataFieldNames)
    stdMax(mm,d) = sD(idxMax(:,mm,d),mm,d);
    classMax{mm,d} = dataStruct.(dataFieldNames{d}).classes{idxMax(:,mm,d)};
    end
end
maxValue = squeeze(maxValue);
%%
classMaxAbr = cell(6,12);
for m = 1:nMag
    classMaxAbr(m,strcmp(classMax(m,:),'Wrist extension')) = {'WE'};
    classMaxAbr(m,strcmp(classMax(m,:),'Wrist flexion')) = {'WF'};
    classMaxAbr(m,strcmp(classMax(m,:),'Supination')) = {'SU'};
    classMaxAbr(m,strcmp(classMax(m,:),'Pronation')) = {'PR'};
    classMaxAbr(m,strcmp(classMax(m,:),'Io')) = {'RD'};
    classMaxAbr(m,strcmp(classMax(m,:),'Tu')) = {'UD'};
    classMaxAbr(m,strcmp(classMax(m,:),'Middle Finger')) = {'ME'};
end
%%
lineColor = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880];
titles = {'ED_p','ED_d','FPL_p','FPL_d','FCU_p','FCU_d'};
% xValues = 1:length(dataFieldNames);
xValues = [13,14,19,26,28,33,34,35,39,40,41];
% xValues = [13,14,19,26,33,34,35,39,40,41];
% maxValue(:,[5,6]) = [];
% stdMax(:,[5,6]) = [];
% classMaxAbr(:,[5,6]) = [];
% xValues = [14,26,28,33,34,35,39,40,41];
figure
for m = 1:nMag
    subplot(3,2,m)
    plot(xValues,maxValue(m,:),'Color',lineColor(m,:))
    hold on
    patch([xValues fliplr(xValues)], [maxValue(m,:)-stdMax(m,:) fliplr(maxValue(m,:)+stdMax(m,:))], lineColor(m,:),'FaceAlpha',0.3,'EdgeColor','none')
    text(xValues,maxValue(m,:)+0.5,classMaxAbr(m,:),'Color',lineColor(m,:))
    hold off
%     xticks(xValues);
    xticks(min(xValues):max(xValues));
%     xticklabels(erase(dataFieldNames,"data"));
    title(titles{m});
    ylabel('Displacement (mm)')
    xlabel('Days from the surgery')
    axis tight
    ylim([0 8])
end
%%
lineColor = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880];
titles = {'ED_p','ED_d','FPL_p','FPL_d','FCU_p','FCU_d'};
classesAvailable = {'Io','Tu','Wrist extension','Wrist flexion','Supination','Pronation',...
    'Thumb flexion','Middle Finger','Flex all fingers'};
% xValues = 1:length(dataFieldNames);
%xValues = [13,14,19,26,27,28,33,34,35,39,40,41];
xValues = [2:6];
Xvalues = [xValues;xValues;xValues;xValues;xValues;xValues;xValues];
cMap = colormap(turbo(9));
% xValues = [13,14,19,26,33,34,35,39,40,41];
% maxValue(:,[5,6]) = [];
% stdMax(:,[5,6]) = [];
% classMaxAbr(:,[5,6]) = [];
% xValues = [14,26,28,33,34,35,39,40,41];
jitter = (rand(length(xValues),length(classesAvailable)) - 0.5).*0.3;
figure()
for m = 1:nMag
    subplot(3,2,m)
%     plot(xValues,squeeze(mD(:,m,:)),'Color',lineColor(m,:))
    hold on
    for c = 1:9
        classData = [];
        stdClassData = [];
        for d = 1:length(dataFieldNames)
            clIdx = strcmp(dataStruct.(dataFieldNames{d}).classes,classesAvailable{c});
            if(sum(clIdx))
                classData = [classData, mD(clIdx,m,d)];
                stdClassData = [stdClassData, sD(clIdx,m,d)];
            else
                classData = [classData, NaN];
                stdClassData = [stdClassData, NaN];
            end
        end
        cData = [nanmean(classData(1:2)), nanmean(classData(3)), nanmean(classData(4:6)),...
            nanmean(classData(7:9)), nanmean(classData(10:12))];
        sData = [std(classData(1:2)), std(classData(3)), std(classData(4:6)),...
            std(classData(7:9)), std(classData(10:12))];
%         sData = [nanmean(stdClassData(1:2)), nanmean(stdClassData(3)), nanmean(stdClassData(4:6)),...
%             nanmean(stdClassData(7:9)), nanmean(stdClassData(10:12))];
        scatter(xValues + jitter(:,c)',cData,120,cMap(c,:),'filled','MarkerFaceAlpha',0.6)
%         errorbar(xValues + jitter(:,c)', cData, sData./2,"LineStyle","none",'Color',cMap(c,:))
    end
    hold off
%     bMatrix = squeeze(mD(:,m,:)); bMatrix(bMatrix == 0) = NaN;
%     boxplot(bMatrix,xValues)
%     scatter(Xvalues,bMatrix)
%     xticks(xValues);
    xticks(min(xValues):max(xValues));
%     xticklabels(erase(dataFieldNames,"data"));
    title(titles{m});
    ylabel('Displacement (mm)')
    xlabel('# weeks after implant')
    axis tight
    ylim([0 8])
    xlim([min(xValues)-1 max(xValues)+1])
    if (m == 1)
        legend(classesAvailable)
    end
end
%%
classMaxCat = categorical(classMax);
figure
for m = 1:nMag
    subplot(3,2,m)
    histogram(classMaxCat(m,:));
    ylim([0 5])
end

%%
figure
for i = 1:7
    hold on
 plot(squeeze(mD(i,4,:)))
end