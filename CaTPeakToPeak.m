edit CaTPeakToPeak
clear
%measure the distance between different sites of Ca2+ influx in cells where
%more than one site was contained within a single scan

%use 'AlignImages.m' to align, then average the images
redChannel=dir('040323_Cell2_Alexa594_Baseline_avg.tif');
greenChannel=dir('040323_Cell2_Fluo5F_Baseline_avg.tif');

pathway='CF_Baseline';

stimOnset=170; %stim onset in ms
samplingRate=0.5; %sampling rate in kHz

pixelSize=0.03086; %pixel size in microns

filterWidth=20; %width of filter used to smooth red channel prior to alignment

fxFilterWidth=1; %filter Width in microns
fxFilterAlpha=0.5; %scales sharpness of filter

timePoints=[1 2 3 4 5]; %scans post stim used to fit with Gaussians

%--------------------------------------------------------------------------
%-Calculate dG/R for average image-----------------------------------------

stimOnset=stimOnset*samplingRate; %convert onset in ms to # of scans

dimensions=imfinfo(redChannel(1).name); %import width/height of images

redImage=double(imread(redChannel(1).name)); %import images
greenImage=double(imread(greenChannel(1).name));

%calculate deltaG
greenBaseline=sum(greenImage(1:stimOnset, :), 1)./stimOnset;
greenBaseline=repmat(greenBaseline, dimensions.Height,1);
deltaG=greenImage-greenBaseline;
imshow(deltaG, [0 20]);
title("deltaG")
figure
imshow(redImage, [0 210])
title("redImage")
%smooth red image to increase signal to noise
redImage=movmean(redImage, filterWidth, 1);
%calculate dG/R
deltaGoverR=deltaG./redImage;
%show image of dG/R as grayscale
figure
imshow(deltaGoverR, [0 0.2]);
title("deltaGoverR")

disp('***Hit Enter to Proceed***')
pause

%--------------------------------------------------------------------------
%-Select scan and fit with Gaussians---------------------------------------

timePoints=timePoints+stimOnset; %assign time points to be plotted

fxdGoverRsweeps=deltaGoverR(timePoints, :); %pull those points out of dG/R

gaussWinFx=gausswin(floor(fxFilterWidth/pixelSize),fxFilterAlpha)...
            /(fxFilterWidth/pixelSize); %construct Gaussian filter

filteredScans=filter(gaussWinFx, 1, fxdGoverRsweeps, [],2); %apply filter
filterOffset=floor(grpdelay(gaussWinFx)); %offset filter delay

filteredScans=filteredScans(:, filterOffset:end); %crop plot to account for
                                                  %filter delay

filteredImageLength=size(filteredScans, 2)*pixelSize;
lengthAxis=0:pixelSize:filteredImageLength-pixelSize; %create x axis (um)

figure
plot(lengthAxis, filteredScans') %plot scans to pick which one to be used

selectScan=inputdlg({'Which scan do you want to fit?'}, 'Select Scan',...
                     [1 15]);
selectScan=str2double(selectScan); %select scan, convert from string

selectedScan=filteredScans(selectScan,:);
plot(lengthAxis, selectedScan', 'Color', 'black',...
                                'LineWidth', 1) %visualize selected scan

%when looking at selected scan estimate number of sites and choose the 
%appropriate number of Gaussians to use in fit
gaussianTypes={'gauss1', 'gauss2', 'gauss3', 'gauss4', 'gauss5', 'gauss6'};
selectGaussian=listdlg('ListString', gaussianTypes);
chosenGaussian=string(gaussianTypes(selectGaussian));

gaussChar=chosenGaussian{1}; %convert string to character array
peakQuant=str2double(gaussChar(6)); %pick out the last char, # of peaks
%create arrays to fill with:
peakConstraints=zeros(peakQuant, 2); %estimated location of each peak 
coeffLowerBounds=zeros(1, peakQuant*3); %lower bounds for coefficients
coeffUpperBounds=repmat([Inf 0 Inf], 1, peakQuant); %upper bounds

for ii=1:peakQuant
    %use cursor to bracket each peak, extract those bounds and insert into
    %arrays to set constraints
    peakLoc=ginput(2);
    coeffLowerBounds(1, (3*ii-1))=peakLoc(1,1);
    coeffUpperBounds(1, (3*ii-1))=peakLoc(2,1);
    
end

%using Least Absolute Residuals to fit, this sets bounds for fit
fitParameters=fitoptions(chosenGaussian,...
                         'Robust', 'LAR',...
                         'Lower', coeffLowerBounds,...
                         'Upper', coeffUpperBounds);
%fit data with Gaussians     
gaussFit=fit(lengthAxis', selectedScan', chosenGaussian, fitParameters);


%--------------------------------------------------------------------------
%-Use coefficient values to calculate peak to peak distance----------------

%get coefficients used in fit (a, b, c; a=amplitude, b=x location, c=width)
gaussCoefficients=coeffvalues(gaussFit);
coeffMatrix=zeros(4, peakQuant);
%move coefficients into an easier to read array
for ii=1:peakQuant
    
    lastColumn=ii*3;
    coeffMatrix(1:3,ii)=gaussCoefficients(1,(lastColumn-2):lastColumn)';
    
end
%convert c to microns
coeffMatrix(4,:)=sqrt(coeffMatrix(3,:).^2./2);
%get rid of repeated values of matrix that has values for b
peakCoord=coeffMatrix(2,1:peakQuant);
peakCoordMatrix=repmat(peakCoord, peakQuant, 1);
peakCoordDiff=nonzeros(abs(peakCoordMatrix-peakCoord'));
peakCoordDiffUnique=unique(peakCoordDiff);

%--------------------------------------------------------------------------
%-Plot Gaussians and Peak to Peak Difference-------------------------------

gaussianPlots=zeros(length(lengthAxis), peakQuant);
%Gauss fit shows the sum of all the Gaussians used in fit
%this uses the coefficients to plot each Gaussian individually
for ii=1:peakQuant
    
    a=coeffMatrix(1, ii);
    b=coeffMatrix(2, ii);
    c=coeffMatrix(3,ii);
    gaussianPlots(:,ii)=a*exp(-((lengthAxis-b)./c).^2);
    
end
%plot scan used in fit
plot(lengthAxis, selectedScan', 'Color', [0.5020    0.5020    0.502],...
                                'LineWidth', 2)
hold on
%plot each Gaussian
for ii=1:peakQuant
    
    plot(lengthAxis', gaussianPlots(:,ii), 'r--')
    hold on
 
end
%plot the summed Gaussian
f=plot(gaussFit);
set(f, 'LineWidth', 2)
%choose peaks to draw lines between and display distance from peaks
examplePeaks=inputdlg({'Peaks to draw arrows from 1:', '2:'}, 'Draw arrows',...
                      [1 10]);
examplePeaks=str2double(examplePeaks);
%select peak with higher amplitude
if coeffMatrix(1,examplePeaks(2))>coeffMatrix(1,examplePeaks(1))
    
    linePlotPeaks(:,1)=coeffMatrix(:,examplePeaks(2));
    linePlotPeaks(:,2)=coeffMatrix(:,examplePeaks(1));
    
else
    
    linePlotPeaks(:,1)=coeffMatrix(:,examplePeaks(1));
    linePlotPeaks(:,2)=coeffMatrix(:,examplePeaks(2));
end
%draw lines at and between the peaks
line([linePlotPeaks(2,1) linePlotPeaks(2,1)],...
     [linePlotPeaks(1,1)-0.005 linePlotPeaks(1,1)+0.005],...
     'Color', 'black', 'LineWidth', 2.5)
line([linePlotPeaks(2,2) linePlotPeaks(2,2)],...
     [linePlotPeaks(1,2)-0.005 linePlotPeaks(1,1)+0.005],...
     'Color', 'black', 'LineWidth', 2.5) 
line([linePlotPeaks(2,1) linePlotPeaks(2,2)],...
     [linePlotPeaks(1,1) linePlotPeaks(1,1)],...
     'Color', 'black', 'LineStyle', '--', 'LineWidth', 2)
%add text showing the distance between the peaks
distance=num2str(abs(linePlotPeaks(2,1)-linePlotPeaks(2,2)));
centerPtoP=linePlotPeaks(2,1)+(linePlotPeaks(2,2)-linePlotPeaks(2,1))/2;
text(centerPtoP, linePlotPeaks(1,1)+0.002,...
     [distance,'\mum'],...
     'Color', 'black', 'FontSize', 12,...
     'HorizontalAlignment', 'center')
     
%
%-Write Files--------------------------------------------------------------
%add time point used to coefficients, some variability from experiment to
%experiment depending on when the peaks are most clear
timeUsed=timePoints(selectScan)/samplingRate;
coeffMatrix(5,:)=timeUsed;

%create a table with labeled rows before saving as txt file
rowNames={'Amplitude (dG/R)', 'Position (um)', 'c', 'sigma (um)',...
          'Time point used (ms)'};

coeffTable=array2table(coeffMatrix, 'RowNames', rowNames);
%use imported file name and string find to name files
expNameRaw=redChannel(1).name;
underscoreLocs=strfind(expNameRaw,'_');
dateAndCell=expNameRaw(1:underscoreLocs(2));
timePointUsed=int2str(timeUsed);
tableName=[dateAndCell pathway '_' timePointUsed '_Gaussian_Fit_Coefficients.txt'];
%write file
writetable(coeffTable, tableName,...
           'WriteRowName', true,...
           'Delimiter', '\t')

%
%-Write red, deltaG, deltaG/R as images------------------------------------
%convert matrices from double to uint16 and save as .tif files

%red image shows red image filtered with movemean function
redImageFileName=[dateAndCell pathway '_' 'RedImage.tif'];
redImageTiff=uint16(redImage);
imwrite(redImageTiff, redImageFileName);

deltaGFileName=[dateAndCell pathway '_' 'deltaGImage.tif'];
deltaGImageTiff=uint16(deltaG);
imwrite(deltaGImageTiff, deltaGFileName);

%added scalar to deltaGoverR so that the raw values were high enough to
%convert to uint16, all 0s without scalar due to low raw values
deltaGoverRFileName=[dateAndCell pathway '_' 'deltaGoverRImage.tif'];
deltaGoverRImageTiff=uint16(deltaGoverR*1000);
imwrite(deltaGoverRImageTiff, deltaGoverRFileName);













