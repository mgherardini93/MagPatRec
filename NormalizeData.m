function [normData, normParams] = NormalizeData(data, type)

switch type
    case 'Midrange0Range2'
        minData = min(data);
        maxData = max(data);
        range = maxData - minData;
        midRange = (maxData + minData)/2;
        normData = (data - midRange)./(range/2);
        normParams.param1 = midRange;
        normParams.param2 = range;
    case 'Mean0Std1'
        meanData = mean(data);
        stdData = std(data);
        normData = (data - meanData)./stdData;
        normParams.param1 = meanData;
        normParams.param2 = stdData;
    otherwise
        normParams.param1 = [];
        normParams.param2 = [];
end

end