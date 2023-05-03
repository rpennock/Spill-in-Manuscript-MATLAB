edit CaT_Localization

%plot location of soma relative to the purkinje cell layer and the site
%where CF CaT was localized on the dendritic arbor

inputLocations=importdata('CF_CaT_Localization.txt');

%Column 1 in txt is the length of a straight line from the top of the PCL
%to the middle of the MLI soma
dSomaFromPCL=inputLocations(:,1);
%Column 2 is a straight line from the middle of the MLI soma to the CaT
%site on the dendrite
%multiple nearby sites (within several microns) are counted as one site
dInputFromSoma=inputLocations(:,2);

%angle of straight line from soma to dendrite, theta=0 for the line
%connecting soma to PCL
angleOfProjection=inputLocations(:,3);
absAngleOfProjection=deg2rad(abs(angleOfProjection));

%use sin and cos to calculate the x and y coord of CaT site
xCoordCaT=sin(absAngleOfProjection).*dInputFromSoma;
%during data entry sites projecting into quadrants 2 and 3 were given
%negative value for angles (e.g. -5 degrees is projecting down and slightly
%into quadrant 3), these were converted to positive values to calculate x
%Coord, this converts them back so they are plotted correctly
xCoordScalar=angleOfProjection./abs(angleOfProjection);
xCoordCaT=xCoordCaT.*xCoordScalar;

yCoordCaT=dSomaFromPCL-dInputFromSoma.*cos(absAngleOfProjection);

%place coordinates in matrix, column 1 and 3 are soma, 3 and 4 are dendrite
finalCoords=zeros(size(inputLocations, 1), 5);
finalCoords(:,2)=dSomaFromPCL;
finalCoords(:,3)=xCoordCaT;
finalCoords(:,4)=yCoordCaT;
finalCoords(:,5)=yCoordCaT-dSomaFromPCL;

%%-Calculate average cell-%%-----------------------------------------------

%avg soma location
avgSomaYCoord=mean(finalCoords(:,2));
sdSomaYCoord=std(finalCoords(:,2));
seSomaYCoord=sdSomaYCoord/sqrt(size(finalCoords, 1));

%avg coordinates of CaT
avgDendYCoord=mean(finalCoords(:,4));
sdDendYCoord=std(finalCoords(:,4));
seDendYCoord=sdDendYCoord/sqrt(size(finalCoords, 1));

avgDendXCoord=mean(finalCoords(:,3));
sdDendXCoord=std(finalCoords(:,3));
seDendXCoord=sdDendXCoord/sqrt(size(finalCoords, 1));

%%-Create Plot-%%----------------------------------------------------------

%distance markers
line([-75 75], [25 25], 'Color', 'white', 'LineStyle', ':',...
                                        'LineWidth', 0.5)
line([-75 75], [50 50], 'Color', 'white', 'LineStyle', ':',...
                                        'LineWidth', 1)
line([-75 75], [75 75], 'Color', 'white', 'LineStyle', ':',...
                                        'LineWidth', 0.5)
line([-75 75], [100 100], 'Color', 'white', 'LineStyle', ':',...
                                        'LineWidth', 1)
%plot lines connecting soma and dendrites
for ii=1:size(finalCoords, 1)

    x1=finalCoords(ii,1); x2=finalCoords(ii,3);
    y1=finalCoords(ii,2); y2=finalCoords(ii,4);
    line([x1 x2], [y1 y2], 'Color', 'red', 'LineStyle', '-',...
                           'LineWidth', 1.5)
    
end

%plot soma
%hold on
%plot(finalCoords(:,1), finalCoords(:,2),'o',...
%                                        'MarkerSize', 12,...
%                                        'MarkerEdgeColor', 'red',...
%                                        'MarkerFaceColor', 'red')

%hold on
%plot(finalCoords(:,1), finalCoords(:,2),'o',...
%                                       'MarkerSize', 10,...
%                                        'MarkerEdgeColor', 'red',...
%                                       'MarkerFaceColor', 'white')
%plot dendrites
%hold on
%plot(finalCoords(:,3), finalCoords(:,4), 'o',...
%                                         'LineWidth', 3,...
%                                         'MarkerSize', 5,...
%                                         'MarkerEdgeColor', 'red',...
%                                         'MarkerFaceColor', 'red')

%mark PCL and add text to distance markers
line([-75 75], [0 0], 'Color', 'white', 'LineStyle', '--',...
                                        'LineWidth',1.5)
line([-75 75], [-20 -20], 'Color', 'white', 'LineStyle', '--',...
                                        'LineWidth',1.5)
text(-49, -10, 'PCL', 'Color', 'white', 'FontSize', 20,...
                    'HorizontalAlignment', 'left')
%text(-49, 30, '25 \mum', 'Color', 'white', 'FontSize', 16,...
%                    'HorizontalAlignment', 'left')
text(-49, 55, '50 \mum', 'Color', 'white', 'FontSize', 16,...
                    'HorizontalAlignment', 'left')
%text(-49, 80, '75 \mum', 'Color', 'white', 'FontSize', 16,...
%                    'HorizontalAlignment', 'left')
text(-49, 105, '100 \mum', 'Color', 'white', 'FontSize', 16,...
                    'HorizontalAlignment', 'left')
n=ii; n=num2str(n);

text(0, 12.5, ['n=', n], 'Color', 'red', 'FontSize', 16,...
                      'HorizontalAlignment', 'center')

%plot average soma and CaT location
%line([0 avgDendXCoord], [avgSomaYCoord avgDendYCoord], 'Color', 'green',...
%     'LineStyle', ':', 'LineWidth', 4)
%line([0 0], [avgSomaYCoord-(sdSomaYCoord/2) avgSomaYCoord+(sdSomaYCoord/2)],...
%     'Color', 'green', 'LineWidth', 3)
%line([avgDendXCoord avgDendXCoord],...
%     [avgDendYCoord-(sdDendYCoord/2) avgDendYCoord+(sdDendYCoord/2)],...
%     'Color', 'green', 'LineWidth', 3)
%line([avgDendXCoord-(sdDendXCoord/2) avgDendXCoord+(sdDendXCoord/2)],...
%     [avgDendYCoord avgDendYCoord],...
%     'Color', 'green', 'LineWidth', 3) 
 
%hold on                
%plot(0, avgSomaYCoord,'o',...
%                      'MarkerSize', 12,...
%                      'MarkerEdgeColor', 'green',...
%                      'MarkerFaceColor', 'green')
%hold on                  
%plot(0, avgSomaYCoord,'o',...
%                      'MarkerSize', 8,...
%                      'MarkerEdgeColor', 'green',...
%                      'MarkerFaceColor', 'white')
%hold on
%plot(avgDendXCoord, avgDendYCoord,'o',...
%                      'MarkerSize', 12,...
%                      'MarkerEdgeColor', 'green',...
%                      'MarkerFaceColor', 'green')
                  
axes=gca;
axes.XLim=[-75 75]; 
axes.YLim=[-30 125];
axes.Color=[0 0 0];
axes.XTick=[];
axes.YTick=[];

%%-Save_txt_file---------------------------------------------------------%%

writematrix(finalCoords, 'CF_CaT_Localization_x_y_xdiff.txt', 'Delimiter',...
            'tab');
type 'CF_CaT_Localization_x_y_xdiff.txt';