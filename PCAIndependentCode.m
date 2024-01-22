dd{1} = load('DD1.mat');
dd{2} = load('DD11.mat');
dd{3} = load('DD2.mat');
dd{4} = load('DD3.mat');
% variableName = {'dd1','dd11','dd2','dd3'};
l = 0;
dataLength = 0;
dataSet = [];
for i = 1:4
    dataSet = [dataSet; dd{i}.dataSet];
    l = l + length(dd{i}.dataSet);
    dataLength = [dataLength, l];
end
%%
coeff = pca(dataSet);
dataSet3DAll = dataSet*coeff(:,1:3);

for i = 1:4
    dd{i}.dataSet3D = dataSet3DAll(dataLength(i) + 1:dataLength(i+1),:);
end

cMap = [];
for i = 1:3
    value = ones(1,3); 
    value = value - 0.2; value(i) = 1;
    cMap = [cMap; value];
    value = ones(1,3); 
    value = value - 0.4; value(i) = 1;
    cMap = [cMap; value];
    value = ones(1,3); 
    value = value - 0.6; value(i) = 1;
    cMap = [cMap; value];
    value = ones(1,3); 
    value = value - 0.8; value(i) = 1;
    cMap = [cMap; value];
end
cMap(7,2) = 0.9; cMap(8,2) = 0.8; %cMap(5:8,[1 3]) = cMap(5:8,[1 3]) - 0.12;

figure
for m = 1:3
    hold on
    for i = 1:4
        dataSet3D = dd{i}.dataSet3D;
        dataOut = dd{i}.dataOut;
        idxDD = find(dataOut == m);
        scatter3(dataSet3D(idxDD,1), dataSet3D(idxDD,2), dataSet3D(idxDD,3),36,cMap(i+(m-1)*4,:),'filled');
    end
%                 scatter3(dataSet3D(dataOut == m,1), dataSet3D(dataOut == m,2), dataSet3D(dataOut == m,3),36,C(m,:),'filled');
%                 scatter(dataSet3D(dataOut == m,1), dataSet3D(dataOut == m,2),36,C(m,:),'filled');
end
hold off
grid on
xlabel('1st Component') 
ylabel('2nd Component')
zlabel('3rd Component')
%             legend(app.magData.mov)
legend('Pronation_1','Pronation_1_1','Pronation_2','Pronation_3','Supination_1','Supination_1_1','Supination_2','Supination_3','Rest_1','Rest_1_1','Rest_2','Rest_3')

%%
dd{1} = load('magData0905.mat');
dd{2} = load('magData1005.mat');
dd{3} = load('magData1105.mat');
%%
dataSet = [];
dataOut = [];
dataSetAll = [];
for i = 1:size(dd,2)
    selFeatures = 'tmd';
    magData = dd{i}.magData;
    featuresStruct = magData.exFeatures;
    nM = magData.nM;
    nCh = magData.nCh*sum(magData.magSelFeatures);
    dataLength = 0;
    for m = 1:nM
        dataLength = dataLength + length(featuresStruct{m});
    end    
    dataSet{i} = zeros(dataLength,nCh);
    dataOut{i} = zeros(dataLength,1);
    dSIdx = 0; % index that indicated at which point to enter data in dataSet
    for m = 1:nM
        movementLength = length(featuresStruct{m});
        for s = 1:movementLength
            dataSet{i}(dSIdx+s,:) = featuresStruct{m}(s).(selFeatures); 
        end
        dataOut{i}((dSIdx+1):(dSIdx + movementLength)) = m;
        dSIdx = dSIdx + movementLength;
    end
    dataSetAll = [dataSetAll; dataSet{i}];
end
l = 0;
dataLength = 0;
for i = 1:size(dd,2)
    l = l + length(dataSet{i});
    dataLength = [dataLength, l];
end
%%
coeff = pca(dataSetAll);
dataSet3DAll = dataSetAll*coeff(:,1:3);

for i = 1:size(dd,2)
    dd{i}.dataSet3D = dataSet3DAll(dataLength(i) + 1:dataLength(i+1),:);
end

cMap = [];
for i = 1:3
    value = ones(1,3); 
    value = value - 0.2; value(i) = 1;
    cMap = [cMap; value];
    value = ones(1,3); 
    value = value - 0.4; value(i) = 1;
    cMap = [cMap; value];
    value = ones(1,3); 
    value = value - 0.6; value(i) = 1;
    cMap = [cMap; value];
    value = ones(1,3); 
    value = value - 0.8; value(i) = 1;
    cMap = [cMap; value];
end
cMap(7,2) = 0.9; cMap(8,2) = 0.8; %cMap(5:8,[1 3]) = cMap(5:8,[1 3]) - 0.12;

figure
for m = 1:3
    hold on
    for i = 1:size(dd,2)
        dataSet3D = dd{i}.dataSet3D;
        dataOutMov = dataOut{i};
        idxDD = find(dataOutMov == m);
        scatter3(dataSet3D(idxDD,1), dataSet3D(idxDD,2), dataSet3D(idxDD,3),36,cMap(i+(m-1)*4,:),'filled');
    end
%                 scatter3(dataSet3D(dataOut == m,1), dataSet3D(dataOut == m,2), dataSet3D(dataOut == m,3),36,C(m,:),'filled');
%                 scatter(dataSet3D(dataOut == m,1), dataSet3D(dataOut == m,2),36,C(m,:),'filled');
end
hold off
grid on
xlabel('1st Component') 
ylabel('2nd Component')
zlabel('3rd Component')
%             legend(app.magData.mov)
legend('Io_0_9','Io_1_0','Io_1_1','Sup_0_9','Sup_1_0','Sup_1_1','Rest_0_9','Rest_1_0','Rest_1_1')
% legend('Io_0_9','Io_1_0','Sup_0_9','Sup_1_0','Rest_0_9','Rest_1_0')
