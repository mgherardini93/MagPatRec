function [onsetSet, onsetOut, threshold, normParams, dNoiseMAV, message] = OnsetDetectorTraining(patRec, dataSet, dataOut, nR)

waitNextPeak = 0.9;
origDataSet = dataSet;
nM = patRec.nM;

% Normalize the training data ---------------------------------------------
[dataSet, normParams] = NormalizeData(dataSet, patRec.normalization.type);
if (isempty(normParams.param1) || isempty(normParams.param2))
    message = 'Error: within normalization';
    return;
end

% Finding indeces of MAV features -----------------------------------------
mavIdx = find(strcmp(patRec.selFeatures,'tmabsTr1'));
if (isempty(mavIdx))
    onsetSet = []; onsetOut = []; threshold = 0; dNoiseMAV = 0;
    message = 'Error: tmabsTr1 feature is not chosen!';
    return;
end
if (~patRec.magSelFeatures(1))
    onsetSet = []; onsetOut = []; threshold = 0; dNoiseMAV = 0;
    message = 'Error: linear distance is not chosen!';
    return;
end

nCh = patRec.nCh;
mavChIdx = (1:nCh) + (mavIdx - 1)*nCh;
if (patRec.magSelFeatures(2)) % if angular distance is chosen remove it 
    mavChIdx = mavChIdx(1):2:mavChIdx(end);
end

onsetSig = mean(dataSet(:,mavChIdx),2); % mean value of MAV feature of all linear distance channels
onsetSigDer = high_pass(onsetSig);      % derivative of the onset signal

if (~sum(sum(patRec.minData)))
    warning('Check if you loaded rest file!');
    disp('Rest values are all 0.');
end

% Calculate rest level ----------------------------------------------------
minRampMAV = abs(patRec.minData);
a = 1;
b = 1/50.*ones(50,1)';
minRampMAV = filter(b,a,minRampMAV);
% minRampMAV = downsample(minRampMAV,25);
minRampMAV = minRampMAV(:,mavChIdx);
minFeatSet = min(origDataSet(:,mavChIdx));
maxFeatSet = max(origDataSet(:,mavChIdx));
scaledMinRamp = 2*(minRampMAV - minFeatSet)./(maxFeatSet - minFeatSet) - 1;
dminRampMAV = high_pass(scaledMinRamp);
dNoiseMAV = 6*std(mean(dminRampMAV(10:end-10,:),2))+mean(mean(dminRampMAV(10:end-10,:),2)); % derivative noise level at rest
disp(dNoiseMAV);

% Find the threshold for the onset detector --------------------------------
thrVals = [];
for j = 1 : nM
    len = sum(dataOut(:,j));
    sig = onsetSigDer(logical(dataOut(:,j)));
    sig = smooth(sig,3); % smooth MAV j movement
    sig(1:5) = 0; % AM 20190117
    minDist = waitNextPeak*(len/nR);
    [pks,~] = findpeaks( sig, 'MinPeakDistance', minDist, 'SortStr', 'descend', 'NPeaks', nR);
    pks(pks<dNoiseMAV) = pks(pks<dNoiseMAV)+100;
    thrVals(1,j) = median(pks);
    
    nIt = 800;
    res = zeros(nIt,2);
    res(:,1) = linspace(dNoiseMAV, thrVals(1,j), nIt)'; %From "dNoise" instead of min(sig) becuase it's the derivative. Creating 200 thresholds from -1 to min peak of j movement

    for i = 1:nIt
        if sig(1) < res(i,1)                                    % the testing threshold has to be higher than -1
            [~,pos] = findpeaks( double(sig > res(i,1)));       % find peaks for signal higher than threshold i
%             [~,pos] = findpeaks(double(sig > res(i,1)),'MinPeakDistance', minDist); %Better results for the derivative
            res(i,2) = numel(pos);                              % count found peaks
        end
    end

    try
        thrVals(2,j) = median(res(res(:,2) == nR, 1)); % finding the min threshold that is able to localize the target onsets 
        nValR = nR;
        while (isnan(thrVals(2,j)) && (nValR > 0))
            nValR = nValR - 1;
            thrVals(2,j) = median(res(res(:,2) == nValR, 1));
        end
    catch
        thrVals(2,j) = NaN;
    end
end
threshold = min(thrVals(2,:)); % treshold for the onset detector

[~, onsets] = findpeaks(double(onsetSigDer > threshold), 'MinPeakDistance', minDist); % the found onsets
nPeaks = numel(onsets); disp(nPeaks);
if numel(onsets) ~= nR * nM
    warning('Onset count does not match number of movements!')
end

% Plot the signal and the positions of the found onsets -------------------
if patRec.plotFigures
    figure; plot(onsetSig,'k-'); hold on; plot(repmat(onsets,1,2)',repmat(ylim,numel(onsets),1)','r--'); hold off;
    figure; plotSignal = plot(onsetSigDer,'k-'); hold on;...
        plotOnsets = plot(repmat(onsets,1,2)',repmat(ylim,numel(onsets),1)','r--'); ...
        plotOnsetThreshold = line([0 length(onsetSig)],[threshold threshold],'color','g');...
        plotNoise = line([0 length(onsetSig)],[dNoiseMAV dNoiseMAV],'color','y');
    ylim([0 max(onsetSigDer)]);
    for mov = 1:nM
        plotMedian = line([round(length(onsetSig)/nM)*(mov-1) round(length(onsetSig)/nM)*mov],[thrVals(1,mov) thrVals(1,mov)],'color','b');
        plotSelectedThreshold = line([round(length(onsetSig)/nM)*(mov-1) round(length(onsetSig)/nM)*mov],[thrVals(2,mov) thrVals(2,mov)],'color','c');
        line([round(length(onsetSig)/nM)*(mov-1) round(length(onsetSig)/nM)*mov],[thrVals(1,mov) thrVals(1,mov)],'color','b');
        line([round(length(onsetSig)/nM)*(mov-1) round(length(onsetSig)/nM)*mov],[thrVals(2,mov) thrVals(2,mov)],'color','c');
    end
    leg = legend([plotSignal plotOnsets(1) plotOnsetThreshold plotNoise plotMedian],...
        'dEMG envelope','Detected Onsets', 'Onset Threshold', 'Rest Noise', 'Median Peaks');
    hold off;
end

% Take only the data rigth after onsets -----------------------------------
onsetSet = dataSet(onsets,:);
onsetOut = dataOut(onsets,:);
message = 'Threshold for the ODA calculated.';

end