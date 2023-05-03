edit AlignImages
clear
%For aligning and cropping line scans from a given experiment
%Uses average fluorescence intensity of Alexa XXX (or whatever you use to
%fill the cells)over a bout of line scans to align images

%User will be prompted input the peak they want to use to align the images
%as well as the number of pixels to keep on either side of the peak
%Code is simply eliminating drift by aligning the cropped images around the
%selected peak

%Cell name should be formated:
%'Date_CellX_Dye_ExperimentalCondition_pClampFile_sweep.tif'
%e.g. '010120_Cell1_Alexa594_NBQX_0001_025.tif'

alexa=dir('041923_Cell1_Alexa594_TBOA_0001*'); %Import images using filenames
caIndicator=dir('041923_Cell1_Fluo5F_TBOA_0001*');

%Works best if limited to line scans taken close in time (e.g. only
%baseline sweeps or sweeps after a drug has taken effect)

stackNameAlexa='041923_Cell1_Alexa594'; %Name stack of images
stackNameCa='041923_Cell1_Fluo5F';      %Keep consistent with original
                                        %file names
                                        
expCondition='TBOA'; %e.g. Control, NBQX, CPP, etc.

rowsUsed=25:250; %Rows averaged to create fluorescence profiles, cut off
                 %end of image
minPeakProm=0.15;
movMeanWidth=50; %width of filter used for smoothing each scan

%No need to change code below here unless you are altering plots
%------------------------------------------------------------------------
%------------------------------------------------------------------------

dimensions=imfinfo(alexa(1).name); %use first image to get dimensions
alexaImages=zeros(dimensions.Height, dimensions.Width, length(alexa),...
                  'uint16'); %preallocate 3D array using sample image dim
caImages=zeros(dimensions.Height, dimensions.Width, length(alexa),...
               'uint16'); %specify 16 bit or images will not be written
                          %correctly later in code
%import images in arrays
for ii=1:length(alexa)
    alexaImages(:,:,ii)=imread(alexa(ii).name);
    caImages(:,:,ii)=imread(caIndicator(ii).name);
end

%smooth images along x axis
alexaImagesSmoothed=movmean(alexaImages, movMeanWidth, 1);

%use sample fluorescence profile to get dimensions for next loop
sampleProfile=sum(alexaImagesSmoothed(rowsUsed,:,1)/length(rowsUsed));
pks=findpeaks(sampleProfile,'MinPeakProminence',minPeakProm);

fluorProfiles=zeros(dimensions.Width,length(alexa));
amp_xPos=zeros(floor((length(pks)/length(alexa)))+3,2,length(alexa));
%create fluor profiles for all images and find peaks
for ii=1:length(alexa)
    fluorProfiles(:,ii)=sum(alexaImagesSmoothed(rowsUsed,:,ii)...
                      /length(rowsUsed));              
    fluorProfiles(:,ii)=fluorProfiles(:,ii)/max(fluorProfiles(:,ii));
    
    [pk, locs]=findpeaks(fluorProfiles(:,ii),...
                         'MinPeakProminence', minPeakProm);
    amp_xPos(1:length(pk),1,ii)=pk;
    amp_xPos(1:length(locs),2,ii)=locs;
    
    findpeaks(fluorProfiles(:,ii),'MinPeakProminence' ,minPeakProm);
    hold on
end
%creates text pointing to peaks
text(amp_xPos(1,2,1)+5, amp_xPos(1,1,1),...
    '\bf \leftarrow Peak1','FontSize', 16);
text(amp_xPos(2,2,1)+5, amp_xPos(2,1,1),...
    '\bf \leftarrow Peak2','FontSize', 16);
text(amp_xPos(3,2,1)+5, amp_xPos(3,1,1),...
    '\bf \leftarrow Peak3','FontSize',16);

hold off

manualInput=questdlg('Do any peaks need to be manually adjusted?');
switch manualInput
    case 'Yes'
       numberOfPeaks=inputdlg({'How many peaks need to be adjusted?'},...,
                               '# of peaks to adjust', [1 20]);
       numberOfPeaks=str2double(numberOfPeaks);
        for ii=1:numberOfPeaks
            peakToAdjust=inputdlg({'Trace to adjust','Peak to Adjust',...
                                   'Pixel value of peak'}, ...
                                   'Manual Peak Input',[1 20; 1 20; 1 20]);
            peakToAdjust=str2double(peakToAdjust);
            amp_xPos(peakToAdjust(2),2,peakToAdjust(1))=peakToAdjust(3);
        
           
        
       
        end
    case 'No'
end
                                

peaks=inputdlg({'Select peak to be used for correction (0 for none)'},...
               'Peak Selection',[1 30]); %dialog box asking for peak
peaks=str2double(peaks); %converts out of dialog to usable integer

if peaks>=1

peakCoordinates=zeros(1, length(alexa));
%import x axis coordinate of peaks into array
for ii=1:length(alexa)
    peakCoordinates(1,ii)=amp_xPos(peaks,2,ii);
end

cropWidth=inputdlg({'Width of imaged left of peak',...
                    'Width of Image right of peak'},...
                   'Crop Image', [1 50; 1 50]); %dialog box asking for             
cropWidth=str2double(cropWidth);                %image width
imageWidth=sum(cropWidth)+1;

croppedAlexaImages=zeros(dimensions.Height,imageWidth,length(alexa),...
                         'uint16');
croppedCaImages=zeros(dimensions.Height,imageWidth,length(alexa),...
                      'uint16');

%crop images according to given parameters, place in 3D array
%save as tiffs

mkdir(expCondition)
addpath(expCondition)

for ii=1:length(alexa)
    
    imageStart=peakCoordinates(ii)-cropWidth(1);
    imageEnd=peakCoordinates(ii)+cropWidth(2);
    
    croppedAlexaImages(:,:,ii)=alexaImages(:,imageStart:imageEnd,ii);
    croppedCaImages(:,:,ii)=caImages(:,imageStart:imageEnd,ii);
    
    if ii<10 %asigns sequential file names to saved images, always 3 digits
        fileNum=strcat('00', num2str(ii));
        
        elseif ii>=10
            fileNum=strcat('0', num2str(ii));
        
            else ii>=100
                fileNum=num2str(ii);
    end

    stackNameAlexaFinal=[stackNameAlexa '_' expCondition '_'...
                         'Corrected' '_' fileNum '.tif'];
    alexaStackDestination=fullfile(expCondition, stackNameAlexaFinal);
    
    imwrite(croppedAlexaImages(:,:,ii),...
            alexaStackDestination,'tif');
  
end

%save Ca2+ indicator images as tiffs, didn't work nesting this in
%other loop
for jj=1:length(alexa)

    if jj<10
        fileNum=strcat('00', num2str(jj));
        
        elseif jj>=10
            fileNum=strcat('0', num2str(jj));
        
            else jj>=100
                fileNum=num2str(jj);
    end

    stackNameCaFinal=[stackNameCa '_' expCondition '_'...
                         'Corrected' '_' fileNum '.tif'];
    caStackDestination=fullfile(expCondition, stackNameCaFinal);
    
    imwrite(croppedCaImages(:,:,jj), caStackDestination, 'tif')
end

else
end

