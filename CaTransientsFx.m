edit CaTransientsFx

clear

%Filter and analyze CaTs in spatial dimension
%Select time points to be analyzed and use curve fitting tool to fit to
%Gaussian curves

%relevant variables were saved and will be loaded into current workspace

%%%---Only for deltaG/R at the moment-----%%%


load('120319_Cell2_30_power_dGoverR.mat');%load parameters from CalculateTransients
%load('031020_Cell1_Control_dGoverG.mat');
load('ExperimentNameDeltaGoverR.mat')
%load('ExperimentNameDeltaGoverG.mat')
load('samplingRate.mat');
load('stimOnset.mat');
load('sweepsWithCaT.mat');
%--------------------------------------------------------------------------
timePostStim=150; %ms post stimulus to analyze fx data
samplePlotTimePoints=[10 11 12 13 14 15 16]; %sample rows to show after filter

fxFilterWidth=1; %filter Width in microns
fxFilterAlpha=0.5; %scales sharpness of filter

pixelSize=0.0291; %pixel size in microns, can find in Prairie metadata
%--------------------------------------------------------------------------
%construct and apply Gaussianfilter to data--------------------------------
gaussWinFx=gausswin(floor(fxFilterWidth/pixelSize),fxFilterAlpha)...
            /(fxFilterWidth/pixelSize);

timePoints=(stimOnset+1):stimOnset+timePostStim*samplingRate;
fxdGoverRsweeps=deltaGoverR(timePoints, :, :);
%fxdGoverGsweeps=deltaGoverG(timePoints, :, :);

filterdGoverR=filter(gaussWinFx,1,fxdGoverRsweeps, [],2);
%filterdGoverG=filter(gaussWinFx,1,fxdGoverGsweeps, [],2);

filterOffset=floor(grpdelay(gaussWinFx));
filtdGoverRoffsetCorr=filterdGoverR(:, filterOffset:end,:);
%filtdGoverGoffsetCorr=filterdGoverG(:, filterOffset:end,:);

%--------------------------------------------------------------------------
%fit data with Gaussian----------------------------------------------------

filteredImageL=size(filtdGoverRoffsetCorr, 2)*pixelSize;
dAlongScan=0:pixelSize:filteredImageL-pixelSize;
sweepsUsed=samplePlotTimePoints;

averageSuccessColumn=size(filtdGoverRoffsetCorr,3)-1;

for ii=sweepsUsed
        
    plot(fxdGoverRsweeps(ii,:,averageSuccessColumn)',...
         'Color', [0.7 0.7 0.7])
    hold on
end

for jj=sweepsUsed
    
    plot(filtdGoverRoffsetCorr(jj,:,averageSuccessColumn)',...
        'LineWidth', 2);
 
    hold on
    
end
hold off

disp('****Hit Enter to Proceed****')
pause

contToPlot=questdlg('Do you want to continue to curve fitting?');

switch contToPlot
    
    case 'Yes'
        
gaussFitCoeff=zeros(4,timePostStim*samplingRate);
subplotDim=ceil(sqrt(length(timePoints)));

figure
for ii=1:timePostStim*samplingRate
    
    gaussFit=fit(dAlongScan',...
                 filtdGoverRoffsetCorr(ii,:,averageSuccessColumn)',...
                 'gauss1');
    gaussFitCoeff(2:4,ii)=coeffvalues(gaussFit)';
    
    subplot(subplotDim, subplotDim, ii)
    plot(dAlongScan',filtdGoverRoffsetCorr(ii,:,...
                                          size(filtdGoverRoffsetCorr,3)-1)')
    hold on
    plot(gaussFit)
    legend off
    hold on
end
hold off

discardSweeps=inputdlg({'Discard fit from sweeps 1 through'},...
                       'Discard poor fit',[1 50]);
discardSweeps=str2double(discardSweeps);
gaussFitCoeff(:,1:discardSweeps)=NaN;    
gaussFitCoeff(4,:)=sqrt(gaussFitCoeff(4,:).^2./2);
gaussFitCoeff(1,:)=timePoints./samplingRate;
coeffNames={'Time (ms)'; 'Amplitude'; 'Peak_Location'; 'Sigma'};
gaussFitSave=[coeffNames, num2cell(gaussFitCoeff)];

figure
plot(gaussFitCoeff(1,:)', gaussFitCoeff(4,:)', '-o')
xlabel('time (ms)', 'FontSize', 14)
ylabel('\sigma (\mum)', 'FontSize', 14)

dGoverRfxtxtfile=filtdGoverRoffsetCorr(:,:,averageSuccessColumn);

txtdGRFx=['dGoverR/' dGoverRexperimentName '_fx_filtered' '.txt'];
writematrix(dGoverRfxtxtfile, txtdGRFx, 'Delimiter', 'tab');
txtdGRFxCoeff=['dGoverR/' dGoverRexperimentName...
               '_fx_filtered_coefficients' '.txt'];
writecell(gaussFitSave, txtdGRFxCoeff, 'Delimiter', 'tab');

    case 'No'
end

