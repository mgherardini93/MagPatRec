% -------------------------- Function Description -------------------------
% This function will calculate the linear and angular distances of all
% pairs of available magnets. 
% It will store the data as:
% (magnet pair 1 linear dist, magnet pair 1 angular dist, magnet pair 2
% linear dist...)
%
% input: 
% - data structure which contains the raw data
% - selFeatures contains the info on which distances should be caclulated
% (first value is for linear, second for angular)
% - flag which allows the calculatation of the displacement of magnet pairs
% from the rest position instead of distances
%
% output:
% - data structure which contains dataLA in which the distances are stored
% - message 
%
function [magData, message] = ExtractLinearAngularDistance(magData,selFeatures,flagDisplacement)

nM = size(magData.rawData,2); % number of movements
nMag = magData.nMag; % number of available magnets
mDim = 6; % data dimension corresponding to one magnet
mPairs = nchoosek(1:nMag,2); % all possible magnet pairs 
nPairs = size(mPairs,1);
laData = cell(1,nM); % temporary variable to store distances
message = 'Error:';
if (nMag ~= 0)
    magnetPoses = magData.rawData;
    if isfield(magData,'restData')
        restData = magData.restData;
    else
        labelRest = magData.label{1};
        restData = magData.rawData{1}(labelRest==0,:);
        restData = restData(1:200,:);
        magData.restData = restData;
    end
    restDataLA = [];
    for i = 1:nPairs
        distR = restData(:,(1:3)+mDim*(mPairs(i,1)-1)) - restData(:,(1:3)+mDim*(mPairs(i,2)-1)); 
        distRN = vecnorm(distR,2,2); % rest linear distance
        angR  = atan2(vecnorm(cross(restData(:,(4:6)+mDim*(mPairs(i,1)-1)),restData(:,(4:6)+mDim*(mPairs(i,2)-1)),2),2,2),dot(restData(:,(4:6)+mDim*(mPairs(i,1)-1)),restData(:,(4:6)+mDim*(mPairs(i,2)-1)),2)); % rest angular distance
        meanRestDist = median(distR,1);
        normMeanRestDist = vecnorm(meanRestDist); % median rest linear distance
        normMeanRestAng = vecnorm(median(angR,1));
        % calculating linear and angular distance for each movement
        for m = 1:nM
            dist = magnetPoses{m}(:,(1:3)+mDim*(mPairs(i,1)-1)) - magnetPoses{m}(:,(1:3)+mDim*(mPairs(i,2)-1));
            temp1 = magnetPoses{m}(:,(4:6)+mDim*(mPairs(i,1)-1));
            temp2 = magnetPoses{m}(:,(4:6)+mDim*(mPairs(i,2)-1));
            ang = atan2(vecnorm(cross(temp1,temp2,2),2,2),dot(temp1,temp2,2));
            distN = vecnorm(dist,2,2);
            if(selFeatures(1)) % if linear distance is chosen
                if (flagDisplacement) % store the displacement
                    laData{m} = [laData{m}, normMeanRestDist - distN];
                    restDataLA = [restDataLA, normMeanRestDist - distRN];
                else                  % store the distance
                    laData{m} = [laData{m}, distN];
                    restDataLA = [restDataLA, distRN];
                end
            end
            if(selFeatures(2)) % if angular distance is chosen
                if (flagDisplacement)
                    laData{m} = [laData{m}, normMeanRestAng - ang];
                    restDataLA = [restDataLA, normMeanRestAng - angR];
                else
                    laData{m} = [laData{m}, ang];
                    restDataLA = [restDataLA, angR];
                end
            end
        end
    end
    magData.magnetPairs = mPairs;
    magData.nCh = nPairs;
    if isempty(laData{1})
        message = 'Error: distances not calculated!';
        return;
    end
    if isempty(restDataLA)
        message = 'Error: rest distances not calculated!';
        return;
    end
    magData.magSelFeatures = selFeatures;
    magData.dataLA = laData;
    magData.restDataLA = restDataLA;
    message = 'Distances successfully calculated.';
else
    message = 'Error: There are no magnets selected!';
end

end