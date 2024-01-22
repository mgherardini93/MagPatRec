function [onsetSet, onsetOut] = OnsetDetectorTest(patRec, dataSet, dataOut, threshold, normParams, dNoiseMAV, nR)

waitNextPeak = 0.9;
nM = patRec.nM;

% Normalize the test data -------------------------------------------------
switch patRec.normalization.type
    case 'Midrange0Range2'
        dataSet = (dataSet - normParams.param1)./(normParams.param2/2);
    case 'Mean0Std1'
        dataSet = (dataSet - normParams.param1)./normParams.param2;
end

% Finding indeces of MAV features -----------------------------------------
mavIdx = find(strcmp(patRec.selFeatures,'tmabsTr1'));

nCh = patRec.nCh;
mavChIdx = (1:nCh) + (mavIdx - 1)*nCh;
if (patRec.magSelFeatures(2)) % if angular distance is chosen remove it 
    mavChIdx = mavChIdx(1):2:mavChIdx(end);
end

onsetSig = mean(dataSet(:,mavChIdx),2); % mean value of MAV feature of all linear distance channels
onsetSigDer = high_pass(onsetSig);      % derivative of the onset signal

% Finding the onsets by checking the threshold ----------------------------
len = sum(dataOut(:,1));
minDist = waitNextPeak*(len/nR);
[~, onsets] = findpeaks(double(onsetSigDer > threshold), 'MinPeakDistance', minDist);
nPeaks = numel(onsets); disp(nPeaks);
if numel(onsets) ~= nR * nM
    warning('Onset count does not match number of movements!')
end

% Plotting the signal and onsets ------------------------------------------
if patRec.plotFigures
    figure; plot(onsetSig,'k-'); hold on; plot(repmat(onsets,1,2)',repmat(ylim,numel(onsets),1)','r--'); hold off;
    figure; plotSignal = plot(onsetSigDer,'k-'); hold on;...
        plotOnsets = plot(repmat(onsets,1,2)',repmat(ylim,numel(onsets),1)','r--'); ...
        plotOnsetThreshold = line([0 length(onsetSigDer)],[threshold threshold],'color','g');...
        plotNoise = line([0 length(onsetSigDer)],[dNoiseMAV dNoiseMAV],'color','y');
    ylim([0 max(onsetSigDer)]);
    hold off;
end

% Taking just the samples after the found onsets --------------------------
onsetSet = dataSet(onsets,:);
onsetOut = dataOut(onsets,:);

end