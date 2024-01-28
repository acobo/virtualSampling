% get 1D sequence from our 2D experiments and also from LEIZA CSV files
% 3dic2023  read both kind of 2D maps (CSV and FIG) and export in CSV format (1D sequences)
% 16dic2023 first try to get virtualSampling(R) working
% "dynamic" version to try dynamic shape change
% 22dic2023 
% 25dic2023 dev04: tiene un bug con lastshape que no se actualiza bien cuando cambia de segmento, fixed
%           dev05: tracking of relative angle with the border
%           dev06: I am having problems with rotate and polyshape (ill-defined, NaNs, duplicate points...),trying with custom code and rotation matrix
% 12ene2024 styles of plots can be changed for publish-ready figures.
% 24ene2024 new chart witht he min/max ratio values along the transverse shape in each spatial point, to have an overview of "possibilities"
% 28ene2024 draw point and shape every N points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%febrero
% Specify the path to your CSV file (2D) or FIG (2D) or leave as '' and a dialog window will ask for the file
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\17082023-130945Nueva medida 2D LAN177\figure24ratio.fig';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\29082023-085942Nueva medida 2D LAN176\figure24ratio.fig';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20012023-094702LAN-178-2D\figure24ratio.fig';
%input_filename= 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\18082023-122613Nueva medida 2D LAN182\figure24ratio.fig'
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\29082023-131210Nueva medida 2D LAN183\figure24ratio.fig';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\24012023-184427LAN-184-2D\figure24ratio.fig';

%agosto
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\17012023-111955_LAN254_seccion-izquierda_lado-anterior_medida_2D\figure24ratio.fig'
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\25012023-150202LAN-255-2D\figure24ratio.fig';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\26012023-091827LAN-256-2D\figure24ratio.fig';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\30012023-181006LAN-258-2D\figure24ratio.fig';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\13012023-105703_LAN259_seccion-izquierda_lado-anterior_medida_2D\figure24ratio.fig'
%input_filename = 'C:\Users\Adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\10012023-165455_LAN261_seccion-izquierda_lado-anterior2D\figure24ratio.fig';

%junio
%input_filename = 'C:\Users\Adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\24082023-113152Nueva medida 2D LAN229\figure24ratio.fig';
%input_filename = 'C:\Users\Adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\25082023-102547Nueva medida 2D LAN230\figure24ratio.fig';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\02022023-193824LAN-231-2D\figure24ratio.fig';
%input_filename = 'C:\Users\Adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\28082023-114458Nueva medida 2D LAN235\figure24ratio.fig';
%input_filename = 'C:\Users\Adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\31012023-162620LAN-236-2D\figure24ratio.fig';
%input_filename = 'C:\Users\Adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\01022023-105931LAN-237-2D\figure24ratio.fig';

%noviembre
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\21022023-094620LAN96-2D\figure24ratio.fig';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20022023-170249LAN98-2D\figure24ratio.fig';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\22082023-130128Nueva medida 2D LAN100\figure24ratio.fig'
%input_filename = 'C:\Users\Adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\22022023-122312LAN105-2D\figure24ratio.fig';
%input_filename = 'C:\Users\Adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\23082023-111426Nueva medida 2D LAN109\figure24ratio.fig';
%LAN122 is not measured

%%%%%%% medidas en Alemania - LEIZA
%febrero
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Lanegre2012_repulidas_MedidasNiklas\febrero\LAN176\LAN176L_20120223_30_ANT_analyzed_data.csv';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Lanegre2012_repulidas_MedidasNiklas\febrero\LAN177\LAN177L_20120223_30_ANT_analyzed_data.csv';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Lanegre2012_repulidas_MedidasNiklas\febrero\LAN184\LAN184L_20120223_30_POST_analyzed_data.csv';
%input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Lanegre2012_repulidas_MedidasNiklas\febrero\LAN184\LAN184L_20120223_30_ANT2_analyzed_data.csv';

%%%%% LAN178
input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20012023-094702LAN-178-2D\figure24ratio.fig';
pathFileToRead = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20012023-094702LAN-178-2D\roiPath_LAN178.mat'; % C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\13012023-105703_LAN259_seccion-izquierda_lado-anterior_medida_2D\ratio_1D_25122023_231724_roiPath_VS.mat'; % "C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\02022023-193824LAN-231-2D\ratio_1D_25122023_230202_roiPath_VS.mat" % 'C:\Users\adolf\Downloads\ratio_1D_25122023_224117_roiPath_VS.mat'; % C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20012023-094702LAN-178-2D\roiPath_LAN178.mat';  %set pathFileToRead=''  to draw, a file to load that file with the linePath
borderFileToRead = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20012023-094702LAN-178-2D\roiBorder_LAN178.mat'; % C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\13012023-105703_LAN259_seccion-izquierda_lado-anterior_medida_2D\ratio_1D_25122023_231724_roiBorder_VS.mat'; % "C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\02022023-193824LAN-231-2D\ratio_1D_25122023_230202_roiBorder_VS.mat" % 'C:\Users\adolf\Downloads\ratio_1D_25122023_224117_roiBorder_VS.mat'; % 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20012023-094702LAN-178-2D\roiBorder_LAN178.mat'; % set to a ROI with the border (isochronous) shape

%%%%% LAN183 paper VirtualSampling
input_filename = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\28022023-115257_LAN183_2D\figure24ratio.fig';
pathFileToRead = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\28022023-115257_LAN183_2D\ratio_1D_again_28012024_184251_roiPath_VS.mat'; %C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\28022023-115257_LAN183_2D\roiPath_LAN178.mat'; % C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\13012023-105703_LAN259_seccion-izquierda_lado-anterior_medida_2D\ratio_1D_25122023_231724_roiPath_VS.mat'; % "C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\02022023-193824LAN-231-2D\ratio_1D_25122023_230202_roiPath_VS.mat" % 'C:\Users\adolf\Downloads\ratio_1D_25122023_224117_roiPath_VS.mat'; % C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20012023-094702LAN-178-2D\roiPath_LAN178.mat';  %set pathFileToRead=''  to draw, a file to load that file with the linePath
borderFileToRead = 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\28022023-115257_LAN183_2D\ratio_1D_again_28012024_184251_roiBorder_VS.mat'; %C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20012023-094702LAN-178-2D\roiBorder_LAN178.mat'; % C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\13012023-105703_LAN259_seccion-izquierda_lado-anterior_medida_2D\ratio_1D_25122023_231724_roiBorder_VS.mat'; % "C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\02022023-193824LAN-231-2D\ratio_1D_25122023_230202_roiBorder_VS.mat" % 'C:\Users\adolf\Downloads\ratio_1D_25122023_224117_roiBorder_VS.mat'; % 'C:\Users\adolf\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\20012023-094702LAN-178-2D\roiBorder_LAN178.mat'; % set to a ROI with the border (isochronous) shape


%%%%% ASK FOR IT
%input_filename = ''; % ask for it

%configuration of virtualSamping procedure
seq_deltaR = 0.05; % spatial resolution of the virtual sampling along the measurement path (deltaR) every 30 micras to match LEIZA resolution
iso_deltaTBorder = 0.05; % spatial resolution of the virtual sampling in the transverse (border shape) curve (deltaT)
seq_firstAbsR = 0; %  first value of the deltaR sequence, only used here for the plot, can be set to other values to match other measurements
actualNpulses = 5; % the average number of LASER pulses in each spatial point is needed to get the number of degree of freedom (pooled STD https://www.isixsigma.com/dictionary/pooled-standard-deviation/m)
maxShapeAngle = 2; % maximum degrees that can be rotated the new shape to search for the minimun RSD
maxShapeXYShift = 0*0.01; % FACTOR that multiplies the distance to the sequence point, to get the maximum shift in x,y for each point in the transverse path (deltaT)
maxAngleExcursion = 15; % this is the maximum difference of each shape angle with respect the starting one at the border, in some limpets, minimum rsd shape could be 90deg and it is wrong, for no effect, put 180
minRsdVariation = 0.02; % this is the minimum value of "rsd variation" to actually changes the shape.  "rsd variation" is the rsd of the rsd vector for all parameters trials in each point, so it is an estimation of how defined are the growth lines around that point, for "shapeless noisy" areas it is better not to change anything. set to 0 to remove this control.
drawEveryNPoints = 5; % set to one to draw all of them
fixedYaxis = [0 0]; % set to min and max values to fix de vertical scale, or [0 0] to do nothing

DEBUG=0;

%% styling
title2D = 'Mg/Ca 2D LIBS map'; % title and x,y labels for 2D charts
xCoord2DLabel = 'X position (mm)'; 
yCoord2DLabel = 'Y position (mm)';
title1D = ''; % title and x,y labels for 1D charts
xCoord1Dlabel = 'Distance to the shell edge (mm)'; %
%yCoord1DLabel = 'Mg/Ca ratio (a.u.)'; % y label is automatically selected as (a.u.) or (mmol/mol) depending on the origin of the measurements
textSize = 15;
lineWidth = 2;

set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', textSize);
set(groot, 'defaultTextFontSize', textSize);
set(groot, 'defaultLineLineWidth', lineWidth);
set(groot, 'defaultAxesLineWidth', lineWidth);
set(groot, 'DefaultaxesFontWeight', 'normal') 
set(groot,'defaultAxesXGrid','off'); 
set(groot,'defaultAxesYGrid','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

haveInterX = dir('interX.m');
haveFillMissing2 = dir('fillmissing2.m');
if isempty(haveInterX) || isempty(haveFillMissing2)
    error('Este script necesita las funciones InterX.m y fillmissing2.m, deberían poder copiarse de la carpeta macros');
end

if isequal(input_filename,'')
    %ask for file
    [file,path] = uigetfile('*.csv;*.fig');
    input_filename = [path file];
end
[fPath, fName, fExt] = fileparts(input_filename);
if isequal(fExt,'.fig') % our measurement
    h_opened_fig_2d=openfig(strcat(input_filename,''));
    % Assume 'figHandle' is your figure handle
    axesHandle = get(h_opened_fig_2d, 'CurrentAxes');  % or simply use gca if it's the current figure

    % Get axes limits
    xLimits = get(axesHandle, 'XLim');
    yLimits = get(axesHandle, 'YLim');

    % Get the children of the axes, which are the line objects
    lines = get(axesHandle, 'Children');

    % Assuming the first child is the one you want
    ejeParam1 = get(lines(1), 'XData');
    ejeParam2 = get(lines(1), 'YData');
    ratio_2D = get(lines(1), 'CData');
    ratio_2D= ratio_2D'; %XXXXXXX CData is trasposed in the figure

    % std errors: they are not in figure, but with lastest script "REPEAT_cflibs_v34_lapas_referencia.m" a .MAT file is generated with the std and relative std values
    % check if a std_error file is in the same folder
    error_filename_wildcard = [ fPath '\procesado_*.mat'];
    error_filenames = dir(error_filename_wildcard);
    if size(error_filenames,1) > 0
        disp(['using std error data from file ' error_filenames(end).name])
        load([ fPath '\' error_filenames(end).name]); % variable is ratio_AB_error_std and ratio_AB_relative_error_std
        errors_abs_2D = ratio_2D .* ratio_AB_relative_error_std; % the absolute STD error is derived from the relative of raw ratios. 
        errors_rel_2D = ratio_AB_relative_error_std; 
    else
        disp('Warning!! no data for std errors found, using zeros');
        errors_abs_2D = zeros(size(ejeParam1,2),size(ejeParam2,2)); % para las barras de error
        errors_rel_2D = zeros(size(ejeParam1,2),size(ejeParam2,2)); % para las barras de error
    end


else % CSV LEIZA format, or exported in CSV format from our experiments
    % Read the CSV file
    disp('Reading CSV file, could take a while...');
    data = readtable(input_filename);
    % Convert table to array
    dataArray = table2array(data);
    % Assign columns to variables
    allx = dataArray(:,1);  % x coordinates
    ally = dataArray(:,2);  % y coordinates
    ratio = dataArray(:,4); % ratio values
    error_abs = dataArray(:,5); % <<<<<<<< este es el número de columna donde está el error, si no es 5, puede ser 6 o 7 ¿?¿?
    error_rel = dataArray(:,6); 
    ejeParam1 = unique(allx)';
    ejeParam2 = unique(ally)';
    ratio_2D = zeros(size(ejeParam1,1),size(ejeParam2,1));
    errors_abs_2D = zeros(size(ejeParam1,2),size(ejeParam2,2)); % para las barras de error
    errors_rel_2D = zeros(size(ejeParam1,2),size(ejeParam2,2)); % para las barras de error

    %ratio_2D = reshape(ratio,[size(ejeParam1,1) size(ejeParam2,1)]) it is not so easy

    for i=1:size(ratio,1)
        [ d, ix ] = min( abs( ejeParam1-allx(i) ) );
        [ d, iy ] = min( abs( ejeParam2-ally(i) ) );
        ratio_2D(ix,iy) = ratio(i);
        errors_abs_2D(ix,iy) = error_abs(i);
        errors_rel_2D(ix,iy) = error_rel(i);
    end

end


% Plot the data as a scatter plot, where the color of each point is determined by the ratio value
h_fig_2d = figure(25);
%h_linePath=scatter(allx, ally, 20, ratio, 'filled')
alpha_channel = ratio_2D>0;
h_data_2d = imagesc(ejeParam1, ejeParam2, ratio_2D','AlphaData', alpha_channel'); % plot is trasposed, but keep data "normal"
set(gca,'YDir','normal');
pbaspect([size(ejeParam1,2) size(ejeParam2,2) 1]); % this is the real aspect ratio
%change caxis
values = nonzeros(sort(reshape(ratio_2D,1,[]))); % sorted values to get axis size without outliers
values = values(~isnan(values));
if ~isempty(values) && floor(size(values,1)*0.98)> 0
    caxis([values(floor(size(values,1)*0.05)) values(floor(size(values,1)*0.95))]); % color axis from x to x % percentile
end

colormap(viridis); % he descargado la función viridis.m de aquí:  https://es.mathworks.com/matlabcentral/answers/1654790-using-a-viridis-colormap
colorbar  % Add a colorbar to the plot
title('2D Scatter Plot of Ratio Values')  % Add a title to the plot
xlabel('X Coordinate')  % Add an x-axis label
ylabel('Y Coordinate')  % Add a y-axis label
set(gcf,'Color','w');


%%%% capture linepath from a 2d image
% v01  8jan2023 variables are "i2s" (image to sequence) to isolate and save for further use in other 2d images
% v02 10jan2023 se puede leer un ROI
%i2s_firstAbsR = 0; % the absolute value of the first point
%figureNumber=24; % <<<<<<<<<<<<<<<<<<< elegir la figura que tiene el mapa 2D
%h_linePath= figure(figureNumber);

ratio_filtered = imgaussfilt(ratio_2D,1);
%ratio_filtered = ratio_2D
ratio_filtered = fillmissing2(ratio_filtered,'nearest'); % NaN at the edges give NaN using interp(), this way NaN are removed

% read or draw the measurement path
if isempty(pathFileToRead)
    ann=annotation('textbox', [0 0.1 0.8 0.1], 'String', 'Draw the (red) measurement path from border clicking on points, double click to end...', 'Color', [1 0.5 0], 'FontWeight', 'bold', 'EdgeColor', 'none');
    roiPath = drawpolyline('Color','red');
    delete(ann);
else
    disp(['Loadig line path already stored...' pathFileToRead]);
    load(strcat(pathFileToRead,''));
    roiPath=drawpolyline('Position',roiPath.Position,'Color','red'); % place the polyline again in the figure
end

% read or draw the BORDER shape (isochronous curve)
if isempty(borderFileToRead)
    ann=annotation('textbox', [0 0.1 0.8 0.1], 'String', 'Draw the (blue) border shape clicking on points, double click to end...', 'Color', [1 0.5 0], 'FontWeight', 'bold', 'EdgeColor', 'none');
    roiBorder = drawpolyline('Color','blue');
    delete(ann);
else
    disp(['Loadig border curve already stored...' borderFileToRead]);
    load(strcat(borderFileToRead,''));
    roiBorder=drawpolyline('Position',roiBorder.Position,'Color','blue'); % place the polyline again in the figure
end

drawnow;
if ~ishandle(h_fig_2d)
    return % if window is closed, let's end the script
end

disp(['Calculating the measurement path using deltaR = ' num2str(seq_deltaR) ' ... ']);
seq_realPos = zeros(0,2); % the same vector as recorded_path but in real position, so we can draw point on the stitched image
seq_realR = zeros(0,1); % values for the R axis
seq_ratios = zeros(0,1); % ratios calculated at each sampling point from the average of isochronous shape
seq_errors_abs = zeros(0,1); % the same for errors
seq_errors_rel = zeros(0,1); % the same for errors
seqVS_ratios = zeros(0,1); % ratios but with the virtual sampling values 
seqVS_errors_abs = zeros(0,1); % the same for errors
seqVS_errors_rel = zeros(0,1); % the same for errors
seqVS_rsdEachPoint = zeros(0,1); % the rsd of each transverse (deltaT) sequence, to compare different VS configurations
seqVS_minEachPoint = zeros(0,1); % minimum value of ratio of each transverse sequence
seqVS_maxEachPoint = zeros(0,1); % maximum value of ratio of each transverse sequence


%%%% crosspoint between the border and the measurement path, if here is no intersection, the path is extrapolated 10% in X
crossPoint = InterX(roiPath.Position',roiBorder.Position'); % this is the intersection of measurement path and border
if isempty(crossPoint)
    % new dev02 extrapolate 10% left and right the path
    minX = min(roiPath.Position(:,1));
    maxX = max(roiPath.Position(:,1));
    extrapMinX = minX - abs((maxX-minX)/10);
    extrapMaxX = maxX + abs((maxX-minX)/10);
    extraMinY =  interp1(roiPath.Position(:,1),roiPath.Position(:,2),extrapMinX,'linear','extrap');
    extraMaxY =  interp1(roiPath.Position(:,1),roiPath.Position(:,2),extrapMaxX,'linear','extrap');
    % thing is, where to append both?
    if roiPath.Position(1,1) < roiPath.Position(end,1) % if the path grows to +x
      extrapRoiPath = [ [extrapMinX , extraMinY] ; roiPath.Position ; [extrapMaxX , extraMaxY]];
    else
      extrapRoiPath = [ [extrapMaxX , extraMaxY] ; roiPath.Position ; [extrapMinX , extraMinY]];
    end
    crossPoint = InterX(extrapRoiPath',roiBorder.Position'); % this is the intersection of EXTRAPOLATED measurement path and border
    if isempty(crossPoint)
        error('The measurement path -red- AND the border shape -blue- MUST CROSS for this algorithm to work. We have tried to extrapolate 10% the measurment path, but still no intersection. Border shape can be drawn *inside* the actual border, no need to be just in the border, we use only the shape.');
    end
end

%% master loop to sample the measurement path in deltaR increments
lastShape= roiBorder.Position + ([roiPath.Position(1,1);roiPath.Position(1,2)]' - crossPoint'); % this is the shape we start on, but it will morph in each deltaR point
keepTrackRelAngle = []; % keep tracking of relative angle with the border, to limit going further that angle
keepRsdVariation = []; % keep tracking of RSD variation (rsd of rsd) in each delta T point, as estimation of how dfined are the growth lines there
pointNumber = 1;
for i=1:size(roiPath.Position,1) - 1
    nextX = roiPath.Position(i+1,1);
    nextY = roiPath.Position(i+1,2);
    if size(seq_realPos,1) == 0 % first point, then prevX, prevY and prevR is already updated bellow
        prevX = roiPath.Position(1,1); % for the first point
        prevY = roiPath.Position(1,2);
        prevR = seq_firstAbsR; % the starting value for the R axis
    end
    if nextX>prevX  dirPathX=1; else dirPathX=-1; end
    if nextY>prevY  dirPathY=-1;else dirPathY=1;  end
    pathDeltaDistance =  sqrt((nextX-prevX)^2+(nextY-prevY)^2); % distance in mm
    pathDeltaArcTan = atan(abs(nextY-prevY)/abs(nextX-prevX));
    npoints = floor(pathDeltaDistance/ seq_deltaR);
    for j=1:npoints
        % we have here prevXY previous point and nextXY next point
        shift = [ seq_deltaR*cos(pathDeltaArcTan)*dirPathX; -seq_deltaR*sin(pathDeltaArcTan)*dirPathY]; % removed: crossPoint - [prevX;prevY];
        seq_realR = [ seq_realR ; prevR];
        seq_realPos = [ seq_realPos ; prevX  , prevY ];
        % now we have to sample but over the border shape
        % this is the starting shape
        %removed roiThisBorder= drawpolyline('Position',lastShape,'Color','white','LineWidth',0.4,'MarkerSize',0.1); % removed, better not shown, no, we need to know the shape against the current growth lines
        if rem(pointNumber,drawEveryNPoints)==0
            drawpoint('Position',[seq_realPos(end,1),seq_realPos(end,2)],'Color','white','MarkerSize',5);
        end
        drawnow;
        %%%%%%%%%%%%%%%%%
        %%%%% this is the inner loop to process along the morphed interpolated spline over the isochronous path  (MISOTIP)
        minX = min(lastShape(:,1)); maxX = max(lastShape(:,1)); minY = min(lastShape(:,2)); maxY = max(lastShape(:,2));
        sevenX = linspace(minX,maxX,7);
        sevenY = spline(lastShape(:,1),lastShape(:,2),sevenX); % spline interpolation of the seven knots
        % sevenDistance = sqrt((sevenY-sevenY(4)).^2 + (sevenX-sevenX(4)).^2); to the central point
        sevenDistance = sqrt((sevenY-prevY).^2 + (sevenX-prevX).^2); % to the current point in the sequence, better than the central point of the spline
        sp7 = spapi(4, sevenX, sevenY);         % to eval for new knots fnval(sp5,xx); % 
        %%%%%% optimization of the spline with change constraints to minimize RSD over the ratio 
        % the control points are in sp7 and the interpolated (delta) points are in interpSBx
        minimunRSD=99999;
        % new dev05 angle excursion could be limited by max difference with border
        currNegMaxShapeAngle = maxShapeAngle;
        currPosMaxShapeAngle = maxShapeAngle;
        if ~isempty(keepTrackRelAngle)% for first run, no problem with max angle excursion, it should be zero
            if (keepTrackRelAngle(end)-maxShapeAngle) < - maxAngleExcursion % max angle excursion exceeded!!
                currNegMaxShapeAngle = 0;
            end
            if (keepTrackRelAngle(end)+maxShapeAngle) > maxAngleExcursion % max angle excursion exceeded!!
                currPosMaxShapeAngle = 0;
            end
        end
        trials= -currNegMaxShapeAngle:0.1:currPosMaxShapeAngle;
        keepRSD = [];
        for t=trials%  trials
            % change randomly control points in x,y considering distance from current crosspoint
            currSevenX = sevenX + randn()*maxShapeXYShift*sevenDistance;
            currSevenY = sevenY + randn()*maxShapeXYShift*sevenDistance;
            % change angle using maxShapeAngle            
            % random change in angle: %por = rotate(po,2*(rand()-0.5)*maxShapeAngle,[prevX,prevY]).Vertices; % rotation about the current point in the deltaR sequence
            rotX = currSevenX; rotY = currSevenY;
            %new dev06 using rotation matrix instead of polyshape rotate
            for csr=1:size(currSevenX,2)
                rotX(csr) =  (currSevenX(csr)-prevX)*cosd(t) + (currSevenY(csr)-prevY)*sind(t) + prevX;
                rotY(csr) = -(currSevenX(csr)-prevX)*sind(t) + (currSevenY(csr)-prevY)*cosd(t) + prevY;
            end
            %por = sortrows([rotX , rotY],1); % order again
            currSevenX= rotX;  % get x,y points back from polyshape used for rotation
            currSevenY = rotY;
            sp7 = spapi(4, currSevenX, currSevenY); % get the spline from the current seven points
            % get the ratios along the current rotated+shifted shape
            interpSBX = linspace(currSevenX(1),currSevenX(end),sqrt((currSevenX(end)-currSevenX(1)).^2+(currSevenY(end)-currSevenY(1)).^2)/iso_deltaTBorder); %
            interpSBY = fnval(sp7,interpSBX);
            % plot trial shape
            %drawpolyline('Position',lastShape,'Color','magenta','LineWidth',0.4,'MarkerSize',0.1);drawnow;
            currRatios = interp2(ejeParam1,ejeParam2,ratio_filtered',interpSBX,interpSBY,'bilinear'); % bilinear is important as must be interpolated, got NaN in the borders and should be removed
            if any(isnan(currRatios))
                %error('NaN in ratios'); % only got NaN if shape is outside the entire figure, but not if outise the measured are due to the fillmissing2 of ratio_2D
            end
            rsd = std(currRatios,'omitnan') ./ mean(currRatios,'omitnan');
            %rsd = (max(currRatios)-min(currRatios))/max(currRatios); % with this approach, there are sudden jumps  with angle
            keepRSD = [ keepRSD , rsd];
            if prevX <11.5 && DEBUG>0
                figure(3);plot(por(:,1),por(:,2));xlim([10 12]);ylim([0.2 1.2]);drawnow;       
                disp(['t=' num2str(t) ' mean=' num2str(mean(currRatios,'omitnan')) ]);
                disp(currRatios);
                %pause;
            end
            if rsd<minimunRSD
                minimunRSD = rsd;
                bestX = currSevenX;
                bestY = currSevenY;
                bestAngle = t;
                if DEBUG disp(['   rsd: ' num2str(rsd)]); end
            end
        end % for t
        if prevX< 11.5 && DEBUG>0
            figure(4);plot(keepRSD);drawnow;    
        end
        %change shape for the better? one
        rsdVariation = std(keepRSD) / mean(keepRSD); % this is an estimation on variation along different shapes to detect defined growth lines versus areas of shapeless noise
        keepRsdVariation = [ keepRsdVariation , rsdVariation];
        if rsdVariation < minRsdVariation % if rsd variation is too low, better not to change anything
            bestX = sevenX;
            bestY = sevenY;
            bestAngle = 0;
        end
        sp7 = spapi(4, bestX, bestY); % new spline for the better points
        interpSBX = linspace(bestX(1),bestX(end),sqrt((bestX(end)-bestX(1)).^2+(bestY(end)-bestY(1)).^2)/iso_deltaTBorder); % get the interpolated values (deltaT)
        interpSBY = fnval(sp7,interpSBX);
        % track of angle
        if isempty(keepTrackRelAngle)
            keepTrackRelAngle = [ keepTrackRelAngle, bestAngle];
        else
            keepTrackRelAngle = [ keepTrackRelAngle, keepTrackRelAngle(end)+bestAngle];
        end
        % plot the final shape
        if rem(pointNumber,drawEveryNPoints)==0
            drawpolyline('Position',[interpSBX ;interpSBY]','Color','white','LineWidth',0.6,'MarkerSize',0.4); % removed, better not shown, no, we need to know the shape against the current growth lines
        end
        iso_realPos = [ interpSBX ; interpSBY]';
        %iso_realPos = [ interpSBX ; interpSBY]';
        if ishandle(h_fig_2d) && 0
            for k=1:size(iso_realPos,1)
                drawpoint('Position',iso_realPos(k,:),'Color','black','MarkerSize',3);
            end
        end
        % at this point, we have a isochronous path with the border shape
        iso_ratios = zeros(size(iso_realPos,1),1);
        iso_errors_abs= zeros(size(iso_realPos,1),1);
        iso_errors_rel= zeros(size(iso_realPos,1),1);
        for n=1:size(iso_realPos,1) % number of sampled points in *this* isochronous shape
            iso_ratios(n) = interp2(ejeParam1,ejeParam2,ratio_2D',iso_realPos(n,1),iso_realPos(n,2),'nearest'); % ojo!!! por defecto es interpolación bilineal, he probado 'nearest'
            iso_errors_abs(n) = interp2(ejeParam1,ejeParam2,errors_abs_2D',iso_realPos(n,1),iso_realPos(n,2),'nearest'); %
            iso_errors_rel(n) = interp2(ejeParam1,ejeParam2,errors_rel_2D',iso_realPos(n,1),iso_realPos(n,2),'nearest'); %
        end
        % and now we have the values of ratio and errors. errors of the "pooled" data can be obtained following  https://www.isixsigma.com/dictionary/pooled-standard-deviation/
        %figure; histogram(iso_ratios);
        avgRatio = mean(iso_ratios,'omitnan');
        minRatio = min(iso_ratios,[],'omitnan');
        maxRatio = max(iso_ratios,[],'omitnan');
        iea = sum((iso_ratios(~isnan(iso_ratios)) - avgRatio).^2);
        avgStdAbs = sqrt(iea / (actualNpulses*size(iso_ratios(~isnan(iso_ratios)),1) -size(iso_ratios(~isnan(iso_ratios)),1))) ; % this is the total number of degrees of freedom (total sample size minus the number of groups)
        avgStdRel = avgStdAbs / avgRatio;
        % original sequence without VS
        seq_ratios = [ seq_ratios ; interp2(ejeParam1,ejeParam2,ratio_2D',seq_realPos(end,1),seq_realPos(end,2),'nearest'); ];
        seq_errors_abs = [ seq_errors_abs ; interp2(ejeParam1,ejeParam2,errors_abs_2D',seq_realPos(end,1),seq_realPos(end,2),'nearest'); ];
        seq_errors_rel = [ seq_errors_rel ; interp2(ejeParam1,ejeParam2,errors_rel_2D',seq_realPos(end,1),seq_realPos(end,2),'nearest'); ];
        seqVS_ratios = [ seqVS_ratios; avgRatio]; % this is the VS value
        seqVS_errors_abs = [ seqVS_errors_abs; avgStdAbs ];
        seqVS_errors_rel = [ seqVS_errors_rel; avgStdRel ];
        seqVS_rsdEachPoint = [seqVS_rsdEachPoint ; minimunRSD ];
        seqVS_minEachPoint = [seqVS_minEachPoint ; minRatio ];
        seqVS_maxEachPoint = [seqVS_maxEachPoint ; maxRatio ];
        
        %update things 
        lastShape = iso_realPos + shift';
        % udpate AFTER recording position both x,y and R
        prevX = prevX + seq_deltaR*cos(pathDeltaArcTan)*dirPathX;
        prevY = prevY - seq_deltaR*sin(pathDeltaArcTan)*dirPathY;
        prevR = prevR + seq_deltaR;
        pointNumber = pointNumber +1;
    end % j npoints
    drawnow;
end % for i

sequence = zeros(1,size(seq_realR,1)); % output sequence with interpolated values
errors_abs = zeros(1,size(seq_realR,1)); % para calcular el error a partir de los valores 2D
errors_rel = zeros(1,size(seq_realR,1)); %


paramRFrom = seq_firstAbsR;
paramRTo = seq_deltaR*(size(seq_realR,1)-1);
ejeParamR = paramRFrom:seq_deltaR:paramRTo;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new values
fnormal=figure(26);
fnormal.Name = strcat('Sequence (no VS) from... ', input_filename);
% prueba de poner barras de error
plot (ejeParamR,seq_ratios, 'Color', 'black','linewidth', 2);
hold on;
errorbar(ejeParamR,seq_ratios,seq_errors_abs,'LineStyle','none', 'Color', '#AAAAAA','linewidth', 1);
%plot (ejeParamR,sequence);
hold off;
xlabel('deltaR on 1D path [mm]');
if isequal(fExt,'.fig')  % from our experiments
    ylabel('Mg/Ca ratio [mmol/mol]');
else
    ylabel('Mg/Ca line intensity ratio [a.u.]');
end
if ~all(fixedYaxis == [0 0])
    ylim(fixedYaxis);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% virtual sampling sequence
fvs=figure(27);
fvs.Name = strcat('Sequence (VirtualSamping) from... ', input_filename);
% prueba de poner barras de error
plot (ejeParamR,seqVS_ratios, 'Color', 'black','linewidth', 2);
hold on;
errorbar(ejeParamR,seqVS_ratios,seqVS_errors_abs,'LineStyle','none', 'Color', '#AAAAAA','linewidth', 1);
%plot (ejeParamR,sequence);
hold off;
xlabel('deltaR on 1D path [mm]');
if isequal(fExt,'.fig')  % from our experiments
    ylabel('VS Mg/Ca ratio [mmol/mol]');
else
    ylabel('VS Mg/Ca line intensity ratio [a.u.]');
end
if ~all(fixedYaxis == [0 0])
    ylim(fixedYaxis);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RSD in each point of the virtual sampling sequence
frsd=figure(28);
frsd.Name = strcat('RSD in each point of the sequence (VirtualSamping) from... ', input_filename);
% prueba de poner barras de error
plot (ejeParamR,100*seqVS_rsdEachPoint, 'Color', 'magenta','linewidth', 2);
xlabel('deltaR on 1D path [mm]');
ylabel('VS RSD of each sequence point [%]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tracking of angle
fangle=figure(29);
fangle.Name = strcat('Tracking of angle from... ', input_filename);
% prueba de poner barras de error
plot (ejeParamR,keepTrackRelAngle, 'Color', 'blue','linewidth', 2);
xlabel('deltaR on 1D path [mm]');
ylabel('VS tracking of angle [deg]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tracking of angle
fangle=figure(30);
fangle.Name = strcat('Tracking of rsd variations from... ', input_filename);
% prueba de poner barras de error
plot (ejeParamR,keepRsdVariation, 'Color', 'black','linewidth', 2);
xlabel('deltaR on 1D path [mm]');
ylabel('VS RSD variations in each deltaR point');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% min/max value of ratio
fminmax=figure(31);
clf;
fminmax.Name = strcat('Min/max values of ratio from... ', input_filename);
% prueba de poner barras de error
hold on;
plot (ejeParamR,seqVS_ratios, 'Color', 'black','linewidth', 2);
plot (ejeParamR,seqVS_minEachPoint, 'Color', [0.5 0.5 0.5], 'linewidth', 2);
plot (ejeParamR,seqVS_maxEachPoint, 'Color', [0.5 0.5 0.5], 'linewidth', 2);
hold off;
xlabel('deltaR on 1D path [mm]');
ylabel('Rango of ratio values in each deltaR point');
if ~all(fixedYaxis == [0 0])
    ylim(fixedYaxis);
end


%export file
data_to_write = [seq_realPos(:,1), seq_realPos(:,2), zeros(size(seq_realPos(:,1))) , seq_ratios ,  seq_errors_abs, seq_errors_rel ];
% Write to CSV file
processingTimeStr = datestr(now,'ddmmyyyy_HHMMSS'); % for plot titles and backup file
[ export_filename,export_path] = uiputfile('*.csv','Choose a name for the 1D CSV file to export', [fPath '\\ratio_1D_' processingTimeStr '.csv' ]);
if export_filename~=0
    [p,f,e]=fileparts(export_filename);
    csvwrite( [ export_path '\\' export_filename ], data_to_write);
    % and now the Virtually sampled sequence
    data_to_write = [seq_realPos(:,1), seq_realPos(:,2), zeros(size(seq_realPos(:,1))) , seqVS_ratios ,  seqVS_errors_abs, seqVS_errors_rel ];
    export_filename =  [f '_1D_VS.csv' ];
    csvwrite( [ export_path '\\' export_filename ], data_to_write);
    export_filename =  [f '_figure_2D_VS.fig' ];
    savefig( h_fig_2d, [ export_path '\\' export_filename ]); % save figure with all the transverse shapes
    export_filename =  [f '_figure_1D_normal.fig' ];
    savefig( fnormal, [ export_path '\\' export_filename ]); % save figure with all the transverse shapes
    export_filename =  [f '_figure_1D_VS.fig' ];
    savefig( fvs, [ export_path '\\' export_filename ]); % save figure with all the transverse shapes
    export_filename =  [f '_figure_rsd_VS.fig' ];
    savefig( frsd, [ export_path '\\' export_filename ]); % save figure with all the transverse shapes
    export_filename =  [f '_roiPath_VS.mat' ];
    save([ export_path '\\' export_filename ],'roiPath');% save roiPath, it could be read back using filpathFileToRead='thisfilename'
    export_filename =  [f '_roiBorder_VS.mat' ];
    save([ export_path '\\' export_filename ],'roiBorder');% save roiPath, it could be read back using filpathFileToRead='thisfilename'
    export_filename =  [f '_figure_minmaxRatio_VS.fig' ];
    savefig( fminmax, [ export_path '\\' export_filename ]); % save figure with the range of values
end

