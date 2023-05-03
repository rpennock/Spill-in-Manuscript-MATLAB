edit CalculateTransients

clear
%Using aligned images find amplitude of Ca2+ transients
%Check images in ImageJ before moving on to Ca2+ transient analysis

%Use images saved to the folder for correct experimental condition
imageInfoRed=dir('041923_Cell1_Alexa594_TBOA_Corrected*'); 
imageInfoGreen=dir('041923_Cell1_Fluo5F_TBOA_Corrected*'); 
                                                         
savetxtFiles='y'; %save .txt files containing CaT(t), yes (y) or no

stimOnset=85; %place of first stimulus in pixels (time dimension)
samplingRate=0.5; %sampling rate in kHz
sweepsWithCaT=1:8; %Enter sweeps that successfully evoked CaT

redSmoothWidth=20; %filter width for red channel

ftWidth=10; %use odd number, pixels taken to calc transient as f(t)

filterWidthFt=10; %width of Gaussian filter used for CaT(t)
filterAlphaFt=4.2; %scalar of Gaussian filter

dGoverRScaleBar=0.05; %__dG/R
dGoverGScaleBar=1; %__dG/G
timeScaleBar=50; %in pixels

%--------------------------------------------------------------------------
%Import red and green images-----------------------------------------------

dimensions=imfinfo(imageInfoRed(1).name);
sweepQuantity=length(imageInfoRed);
sweepQuantitySuccess=length(sweepsWithCaT);

redChannel=zeros(dimensions.Height, dimensions.Width,...
                 sweepQuantity, 'uint16');
greenChannel=zeros(dimensions.Height, dimensions.Width,...
                 sweepQuantity, 'uint16');
              
for ii=1:sweepQuantity
    
    redChannel(:,:,ii)=imread(imageInfoRed(ii).name);
    greenChannel(:,:,ii)=imread(imageInfoGreen(ii).name);
    
end

%--------------------------------------------------------------------------
%Calculate dG/R and dG/G for entire images---------------------------------
redChannelSmooth=zeros(dimensions.Height, dimensions.Width,...
                       sweepQuantity);
                   
backgroundGreen=zeros(dimensions.Height, dimensions.Width,...
                      sweepQuantity);
deltaG=zeros(dimensions.Height, dimensions.Width,...
             sweepQuantity);
         
deltaGoverR=zeros(dimensions.Height, dimensions.Width,...
                  sweepQuantity);
deltaGoverG=zeros(dimensions.Height, dimensions.Width,...
                  sweepQuantity);

for ii=1:sweepQuantity

    greenTemp=double(greenChannel(:,:,ii));
    greenTempAvg=sum(greenTemp(1:stimOnset,:), 1)/stimOnset;
    
    backgroundGreen(:,:,ii)=repmat(greenTempAvg,dimensions.Height,1);
    
    deltaG(:,:,ii)=double(greenTemp)-backgroundGreen(:,:,ii);
    deltaGoverG(:,:,ii)=deltaG(:,:,ii)./backgroundGreen(:,:,ii);
    
    redTemp=double(redChannel(:,:,ii));
    redChannelSmooth(:,:,ii)=movmean(redTemp, redSmoothWidth, 1);
    deltaGoverR(:,:,ii)=deltaG(:,:,ii)./redChannelSmooth(:,:,ii);
    
end

%--------------------------------------------------------------------------
%separate sweeps with failures---------------------------------------------

dGoverRsuccess=zeros(dimensions.Height, dimensions.Width,...
                     sweepQuantitySuccess);
dGoverGsuccess=zeros(dimensions.Height, dimensions.Width,...
                     sweepQuantitySuccess);
              
dGoverRfailures=zeros(dimensions.Height, dimensions.Width,...
                      sweepQuantity-sweepQuantitySuccess);
dGoverGfailures=zeros(dimensions.Height, dimensions.Width,...
                      sweepQuantity-sweepQuantitySuccess);

for ii=1:sweepQuantity
    if ii<=sweepQuantitySuccess
        dGoverRsuccess(:,:,ii)=deltaGoverR(:,:,ii);
        dGoverGsuccess(:,:,ii)=deltaGoverG(:,:,ii);
        
        elseif ii>sweepQuantitySuccess
            index=ii-sweepQuantitySuccess;
            dGoverRfailures(:,:,index)=deltaGoverR(:,:,ii);
            dGoverGfailures(:,:,index)=deltaGoverG(:,:,ii);
    end
end

%--------------------------------------------------------------------------
%Average separated traces--------------------------------------------------

dGoverRsuccessAvg=sum(dGoverRsuccess,3)/sweepQuantitySuccess;
dGoverGsuccessAvg=sum(dGoverGsuccess,3)/sweepQuantitySuccess;

dGoverRfailuresAvg=sum(dGoverRfailures,3)/(sweepQuantity...
                                           -sweepQuantitySuccess);
dGoverGfailuresAvg=sum(dGoverGfailures,3)/(sweepQuantity...
                                           -sweepQuantitySuccess);

%--------------------------------------------------------------------------
%plots average dG/R, ask how many influx sites are present-----------------

imshow(dGoverRsuccessAvg, [0 0.1])
axis on
hold on
line ([1, dimensions.Width], [stimOnset-2, stimOnset-2],...
      'color', 'g', 'LineStyle', '--', 'LineWidth', 2);
for ii=25:25:(floor(dimensions.Width/25)-1)*50
    line([ii, ii], [1, dimensions.Height],...
      'color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
end
hold off

disp('****Hit Enter to Proceed****')
pause

transientPeakNum=inputdlg({'How many sites are present?'},...
                          'Number of Sites', [1 50]);
transientPeakNum=str2double(transientPeakNum);
transientPeakLoc=zeros(1, transientPeakNum);

%--------------------------------------------------------------------------
%Calculate CaT(t) for dG/R and dG/G) and apply Gaussian filter-------------

deltaGoverRft=zeros(dimensions.Height, sweepQuantity,...
                    transientPeakNum);
deltaGoverRftFiltered=zeros(dimensions.Height, sweepQuantity,...
                    transientPeakNum);
dGoverRftAvg=zeros(dimensions.Height, 2, transientPeakNum);
dGoverRftAvgFiltered=zeros(dimensions.Height, 2, transientPeakNum);

deltaGoverGft=zeros(dimensions.Height, sweepQuantity,...
                    transientPeakNum);
deltaGoverGftFiltered=zeros(dimensions.Height, sweepQuantity,...
                    transientPeakNum);
dGoverGftAvg=zeros(dimensions.Height, 2, transientPeakNum);
dGoverGftAvgFiltered=zeros(dimensions.Height, 2, transientPeakNum);

gaussianWindowft=gausswin(filterWidthFt,filterAlphaFt);
filterDelay=mean(ceil(grpdelay(gaussianWindowft)));

dGoverRftAvgFiltDelayCorr=zeros(dimensions.Height-filterDelay,...
                                2, transientPeakNum);
dGoverRftFiltDelayCorr=zeros(dimensions.Height-filterDelay,...
                            sweepQuantity,transientPeakNum);...
dGoverGftAvgFiltDelayCorr=zeros(dimensions.Height-filterDelay,...
                                2, transientPeakNum);
dGoverGftFiltDelayCorr=zeros(dimensions.Height-filterDelay,...
                            sweepQuantity,transientPeakNum);...

for ii=1:transientPeakNum
    
    dialog=sprintf('Center of site %d on x axis (estimate)', ii);
    location=inputdlg({dialog}, 'Center of signal', [1 50]);
    transientPeakLoc(1,ii)=str2double(location);
    
    ftColumnsLeft=transientPeakLoc(1,ii)-(ftWidth-1)/2;
    ftColumnsRight=transientPeakLoc(1,ii)+(ftWidth-1)/2;
    
    for jj=1:sweepQuantity
        %dG/R, calculate CaT(t), apply filter and correct for delay
        tempdGoverR=deltaGoverR(:,:,jj);
        transientCalc=sum(tempdGoverR(:,...
                          ftColumnsLeft:ftColumnsRight),2)...
                          /ftWidth;
        deltaGoverRft(:,jj,ii)=transientCalc;
        deltaGoverRftFiltered(:,jj,ii)=filter(gaussianWindowft,...
                                              filterWidthFt/3.8,...
                                              deltaGoverRft(:,jj,ii),[],1);
       
        dGoverRftFiltDelayCorr(:,jj,ii)=deltaGoverRftFiltered...
                                        (filterDelay:dimensions.Height-1,...
                                         jj,ii);  
        %dG/G                             
        tempdGoverG=deltaGoverG(:,:,jj);
        transientCalcGoverG=sum(tempdGoverG(:,...
                            ftColumnsLeft:ftColumnsRight),2)...
                            /ftWidth;
        deltaGoverGft(:,jj,ii)=transientCalcGoverG;
        deltaGoverGftFiltered(:,jj,ii)=filter(gaussianWindowft,...
                                              filterWidthFt/3.8,...
                                              deltaGoverGft(:,jj,ii),[],1);
       
        dGoverGftFiltDelayCorr(:,jj,ii)=deltaGoverGftFiltered...
                                        (filterDelay:dimensions.Height-1,...
                                         jj,ii);
    end
    
    %dG/R, calculate mean CaT(t), separate sweeps with failures
    tempdGoverRii=deltaGoverRft(:,:,ii);
    dGoverRftAvg(:,1,ii)=sum(tempdGoverRii(:,1:sweepQuantitySuccess),2)...
                             /sweepQuantitySuccess;
    dGoverRftAvg(:,2,ii)=sum(tempdGoverRii(:,...
                             sweepQuantitySuccess+1:...
                             sweepQuantity),2)...
                             /(sweepQuantity-sweepQuantitySuccess);  
                         
     dGoverRftAvgFiltered(:,1,ii)=filter(gaussianWindowft,...
                                         filterWidthFt/3.8,...
                                         dGoverRftAvg(:,1,ii),[],1);
     dGoverRftAvgFiltered(:,2,ii)=filter(gaussianWindowft,...
                                         filterWidthFt/3.8,...
                                         dGoverRftAvg(:,2,ii),[],1);                     
     
     dGoverRftAvgFiltDelayCorr(:,1,ii)=dGoverRftAvgFiltered...
                                        (filterDelay:dimensions.Height-1,...
                                         1,ii);
     dGoverRftAvgFiltDelayCorr(:,2,ii)=dGoverRftAvgFiltered...
                                        (filterDelay:dimensions.Height-1,...
                                         2,ii);
    tempdGoverGii=deltaGoverGft(:,:,ii);
    dGoverGftAvg(:,1,ii)=sum(tempdGoverGii(:,1:sweepQuantitySuccess),2)...
                             /sweepQuantitySuccess;
    dGoverGftAvg(:,2,ii)=sum(tempdGoverGii(:,...
                             sweepQuantitySuccess+1:...
                             sweepQuantity),2)...
                             /(sweepQuantity-sweepQuantitySuccess);  
                         
     dGoverGftAvgFiltered(:,1,ii)=filter(gaussianWindowft,...
                                         filterWidthFt/3.8,...
                                         dGoverGftAvg(:,1,ii),[],1);
     dGoverGftAvgFiltered(:,2,ii)=filter(gaussianWindowft,...
                                         filterWidthFt/3.8,...
                                         dGoverGftAvg(:,2,ii),[],1);                     
     
     dGoverGftAvgFiltDelayCorr(:,1,ii)=dGoverGftAvgFiltered...
                                        (filterDelay:dimensions.Height-1,...
                                         1,ii);
     dGoverGftAvgFiltDelayCorr(:,2,ii)=dGoverGftAvgFiltered...
                                        (filterDelay:dimensions.Height-1,...
                                         2,ii);
end

%-------------------------------------------------------------------------
%plot CaT(t) from each sweep with average overlaid------------------------

if transientPeakNum>1 %choose site to plot if more than 1 present
    dataSet=inputdlg({'Plot transients from site:'}, 'Plot CaT(t)', [1 30]);
    dataSet=str2double(dataSet);
    
   else
    dataSet=1;
end
tAxis=0:samplingRate:((dimensions.Height-(filterDelay+1))*samplingRate);

%creates dialog to choose whether to plot dG/R or dG/G
transientUsed=questdlg('Do you want to plot dG/R or dG/G?',...
                       'Choose data to plot',...
                       'dG/R', 'dG/G','dG/R');
switch transientUsed
    case 'dG/R' %Assigns dG/R data to plots
        sweepsToPlot=dGoverRftFiltDelayCorr(:,:,dataSet);
        sweepstoPlotAvgS=dGoverRftAvgFiltDelayCorr(:,1,dataSet);
        sweepstoPlotAvgF=dGoverRftAvgFiltDelayCorr(:,2,dataSet);
        yScaleBar=dGoverRScaleBar;
        dGoverX='\DeltaG/R';
        
    case 'dG/G' %Assigns dG/G data to plots
        sweepsToPlot=dGoverGftFiltDelayCorr(:,:,dataSet);
        sweepstoPlotAvgS=dGoverGftAvgFiltDelayCorr(:,1,dataSet);
        sweepstoPlotAvgF=dGoverGftAvgFiltDelayCorr(:,2,dataSet);
        yScaleBar=dGoverGScaleBar;
        dGoverX='\DeltaG/G';
end

figure
plot(tAxis, sweepsToPlot, 'Color', [0.6 0.6 0.6],... %plots all sweeps
                          'LineWidth', 0.3)
hold on

f=plot(tAxis, sweepstoPlotAvgF, ... %average of failures
                                'Color', [0.4 0.4 0.4],...
                                'LineWidth', 2,...
                                'DisplayName', 'Failures');
hold on
s=plot(tAxis, sweepstoPlotAvgS,... %average of successes
                               'Color', [0.75 0 0],...
                               'LineWidth', 2,...
                               'DisplayName', 'Successes');
hold off

tAxisEnd=max(tAxis); %calculate origin of time scale bar
dFromtAxisEnd=tAxisEnd*0.02;
scaleBarOriginTime=(tAxisEnd-dFromtAxisEnd)-timeScaleBar;

yAxisEnd=max(max(sweepsToPlot)); %dGoverR scale bar
dFromyAxisEnd=yAxisEnd*0.02;
scaleBarOriginCa=(yAxisEnd-dFromyAxisEnd)-yScaleBar;

%draw scale bars
line([scaleBarOriginTime scaleBarOriginTime],... %draw dGoverR scale bar
     [scaleBarOriginCa scaleBarOriginCa+yScaleBar],...
     'Color', 'black',...
     'LineWidth', 2);
line([scaleBarOriginTime scaleBarOriginTime+timeScaleBar*samplingRate],... %time
     [scaleBarOriginCa scaleBarOriginCa],...
     'Color', 'black',...
     'LineWidth', 2);
%label scale bars
percdGoverX=yScaleBar*100; 

dGoverXscaleLabelt=scaleBarOriginTime-(0.15*(timeScaleBar*samplingRate));
dGoverXscaleLabely=scaleBarOriginCa+(0.9*yScaleBar);

text(dGoverXscaleLabelt,dGoverXscaleLabely,... 
                                  [num2str(percdGoverX),'%', dGoverX],... 
                                  'FontSize', 12,...
                                  'HorizontalAlignment', 'right',...
                                  'Rotation', 90);
                              
timeScaleLabelt=scaleBarOriginTime+(0.75*...
                                   (timeScaleBar*samplingRate)); %placement of scale
timeScaleLabely=scaleBarOriginCa-(0.15*yScaleBar);               %bars
                              
text(timeScaleLabelt, timeScaleLabely,... 
                                  [num2str(timeScaleBar/samplingRate),' ms'],...
                                  'FontSize', 12,...
                                  'HorizontalAlignment', 'right');

plotTitle=sprintf('Site %d', dataSet); %site plotted                           
title(plotTitle)
axes=gca; %turn off axes
axes.XColor='none';
axes.YColor='none';
legend([f s], 'Location', 'southeast') %turn on legend
legend('boxoff')

%--------------------------------------------------------------------------
%Save files as .mat and .txt files as needed-------------------------------

if savetxtFiles=='Y'||'y'

mkdir('dGoverR');
mkdir('dGoverG');
addpath('dGoverR','dGoverG');

for ii=1:transientPeakNum 
        
    txtNameRaw=[imageInfoRed(1).name];
    underScoreLoc=strfind(txtNameRaw,'_');
    underScoreQuant=length(underScoreLoc);
    periodLoc=strfind(txtNameRaw,'.');
    dateAndCell=txtNameRaw(1:underScoreLoc(2)-1);
    expCond=txtNameRaw(underScoreLoc(3):...
                           underScoreLoc(underScoreQuant-1)-1);
                       
    siteNum=sprintf('_Site %d', ii);  
    
    txtNamedGR=['dGoverR/' dateAndCell expCond '_dGoverR_ft'...
                siteNum '.txt'];
    txtNamedGRfilt=['dGoverR/' dateAndCell expCond '_dGoverR_ft'...
                    '_filtered' siteNum '.txt'];
    txtNamedGG=['dGoverG/' dateAndCell expCond '_dGoverG_ft'...
                siteNum '.txt'];
    txtNamedGGfilt=['dGoverG/' dateAndCell expCond '_dGoverG_ft'...
                    '_filtered' siteNum '.txt'];
    
    dGoverRcombinedTraces(:,1:sweepQuantity)=deltaGoverRft(:,:,ii);
    dGoverRcombinedTraces(:,sweepQuantity+1:sweepQuantity+2)=dGoverRftAvg...
                                                               (:,:,ii);
    dGoverRcombinedTracesFilt(:,1:sweepQuantity)=dGoverRftFiltDelayCorr...
                                                               (:,:,ii);
    dGoverRcombinedTracesFilt(:,sweepQuantity+1:sweepQuantity+2)=...
                                               dGoverRftAvgFiltDelayCorr...
                                                               (:,:,ii);
                                                            
    dGoverGcombinedTraces(:,1:sweepQuantity)=deltaGoverGft(:,:,ii);
    dGoverGcombinedTraces(:,sweepQuantity+1:sweepQuantity+2)=dGoverGftAvg...
                                                               (:,:,ii);
    dGoverGcombinedTracesFilt(:,1:sweepQuantity)=dGoverGftFiltDelayCorr...
                                                               (:,:,ii);
    dGoverGcombinedTracesFilt(:,sweepQuantity+1:sweepQuantity+2)=...
                                               dGoverGftAvgFiltDelayCorr...
                                                               (:,:,ii);                                                          
        
    writematrix(dGoverRcombinedTraces, txtNamedGR, 'Delimiter', 'tab');
    writematrix(dGoverRcombinedTracesFilt, txtNamedGRfilt,...
                 'Delimiter', 'tab');
             
    writematrix(dGoverGcombinedTraces, txtNamedGG, 'Delimiter', 'tab');
    writematrix(dGoverGcombinedTracesFilt, txtNamedGGfilt,...
                 'Delimiter', 'tab');
end

workSpaceSaveGR=[dateAndCell expCond '_dGoverR' '.mat'];
workSpaceSaveGG=[dateAndCell expCond '_dGoverG' '.mat'];
dGoverRexperimentName=[dateAndCell expCond '_dGoverR'];
dGoverGexperimentName=[dateAndCell expCond '_dGoverG'];

deltaGoverR(:,:,sweepQuantity+1)=dGoverRsuccessAvg;
deltaGoverR(:,:,sweepQuantity+2)=dGoverRfailuresAvg;
deltaGoverG(:,:,sweepQuantity+1)=dGoverGsuccessAvg;
deltaGoverG(:,:,sweepQuantity+2)=dGoverGfailuresAvg;

save('stimOnset.mat', 'stimOnset')
save('samplingRate.mat', 'samplingRate')
save('sweepsWithCaT.mat', 'sweepsWithCaT')
save(workSpaceSaveGR, 'deltaGoverR')
save(workSpaceSaveGG, 'deltaGoverG')
save('ExperimentNameDeltaGoverR.mat', 'dGoverRexperimentName')
save('ExperimentNameDeltaGoverG.mat', 'dGoverGexperimentName')

else
end


