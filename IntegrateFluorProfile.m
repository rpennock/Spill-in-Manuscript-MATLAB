edit IntegrateFluorProfile

clear

%script calculates cumulative fluorescence for line scan images from two
%pathways on same dendrite
%compares cumulative fluorescence using KS test

CF=struct2array(load('020221_Cell1_CF_dGoverR.mat'));
PF=struct2array(load('020221_Cell1_PF_dGoverR.mat'));

fxFilterWidth=1; %filter Width in microns
fxFilterAlpha=0.5; %scales sharpness of filter
pixelSize=0.0349; %pixel size in microns, can find in Prairie metadata

peakPosition=[5.76510624431158 5.32082013066293]; %[CF PF], position
sigma=[0.791296960969077 1.79788651644356]; %[CF PF], sigma

peakPosCorrection=1; %in pixels
                      %the cropping of the images is not always consistent 
                      %between the array and when Gaussian fits were
                      %calculated. This correction will adjust the peak
                      %position that the peaks are in the center of the
                      %fluorescence profile
timePointUsed=180; % in ms
samplingRate=0.5; %in kHz

cutoffSigmaLeft=1.5; %distance lateral to peak to include in analysis
cutoffSigmaRight=2.5; %multiple of the sigma of the Gaussian fit
%----------------------------------------
%-Extract average dG/R array

avgProfileIndex=size(CF,3)-1; %Average of successful sweeps on second to last page
CFAvgProfile=CF(:,peakPosCorrection:end,avgProfileIndex); 
PFAvgProfile=PF(:,peakPosCorrection:end,avgProfileIndex); 

gaussWinFx=gausswin(floor(fxFilterWidth/pixelSize),fxFilterAlpha)...
            /(fxFilterWidth/pixelSize); %Gaussian kernel
        
filteredCFAvgProfile=filter(gaussWinFx, 1, CFAvgProfile, [], 2); %apply Gaussian filter
filteredPFAvgProfile=filter(gaussWinFx, 1, PFAvgProfile, [], 2);

filterOffset=floor(grpdelay(gaussWinFx)); %calculate group delay from filter and offset profiles
filterOffsetCFAvgProfile=filteredCFAvgProfile(:, filterOffset:end);
filterOffsetPFAvgProfile=filteredPFAvgProfile(:, filterOffset:end);


timePointUsedRow=timePointUsed*samplingRate; %calculate the row of array to use

%determine which profile has a peak further to the left on the profile
%calculate the cutoff to the left to be used for KS test
[minPosition sigmaMinIndex]=min(peakPosition);
minSigma=sigma(sigmaMinIndex);
leftCutoff=floor((minPosition-(minSigma*cutoffSigmaLeft))/pixelSize);

%same but to the right
[maxPosition sigmaMaxIndex]=max(peakPosition);
maxSigma=sigma(sigmaMaxIndex);
rightCutoff=floor((maxPosition+(maxSigma*cutoffSigmaRight))/pixelSize);

xAxisValues=leftCutoff*pixelSize:pixelSize:rightCutoff*pixelSize;

%plot the cropped fluorescence profiles for both pathways
plot(xAxisValues, filterOffsetCFAvgProfile(timePointUsedRow,leftCutoff:rightCutoff))
set(gca, 'TickDir', 'out')
hold on
plot(xAxisValues, filterOffsetPFAvgProfile(timePointUsedRow,leftCutoff:rightCutoff))

%calculat the cumulative sum for each profile in the range indicated
intCF=cumsum(filterOffsetCFAvgProfile(timePointUsedRow,leftCutoff:rightCutoff));
intPF=cumsum(filterOffsetPFAvgProfile(timePointUsedRow,leftCutoff:rightCutoff));

%normalize the cum sum
intNormCF=intCF/max(intCF);
intNormPF=intPF/max(intPF);

%plot normalized cum sum for both pathways
figure
plot(xAxisValues, intNormCF)
set(gca, 'TickDir', 'out')
hold on
plot(xAxisValues, intNormPF)

%two sample KS test, default alpha=0.05
[boolian p]=kstest2(intNormCF, intNormPF)

rowNames={'Left Peak', 'Left Sigma', 'Left Sigma Cutoff',...
          'Right Position', 'Right Sigma' 'Right Sigma Cutoff', 'p value'};
analysisValues={minPosition, minSigma, cutoffSigmaLeft,...
                maxPosition, maxSigma, cutoffSigmaRight, p};
analysisValuesTable(:,1)=rowNames';
analysisValuesTable(:,2)=analysisValues';

writecell(analysisValuesTable, 'KS_test.txt', 'Delimiter', 'tab')

