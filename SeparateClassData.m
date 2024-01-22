% -------------------------- Function Description -------------------------
% This function will store only the contraction part of every movement.
% Incase the flagAddRest is on, the rest part between the contractions will
% be kept as well as a separate class.
%
% input: 
% - data structure which contains dataTreated
% - flag which decides if rest class needs to be created from the data
% inbetween the movement contractions
%
% output:
% - data structure which contains dataTreated with separated classes
% - message 
%
function [dataStruct, message] = SeparateClassData(dataStruct, flagAddRest)
try
    data = dataStruct.dataTreated;
    label = dataStruct.labelTreated;
    nM = size(dataStruct.rawData,2);
    nMov = nM;
    if (flagAddRest) % in case the rest need to be added, number of movements needs to be increased
        nMov = nMov + 1;
        restData = [];
    end
    dataTreated = cell(1,nMov);
    
    % note: rest is always labeled as 0
%     if (flagAddRest)
%         restData = [restData; data{2}(label{2} == 0,:)];
%     end
    for m = 1:nM
        dataTreated{m} = data{m}(label{m} ~= 0,:);
        if (flagAddRest)
            restData = [restData; data{m}(label{m} == 0,:)];
        end
    end
    if (flagAddRest)
        dataTreated{nMov} = restData;
        dataStruct.nM = nMov; % update the nM variable
    end
    dataStruct.dataTreated = dataTreated;
    message = 'Classes separated.';
catch
    message = 'Error: In separating classes!';
end

end