%%% This programe requires two colored Z-stack of Airyscan image containing 40 zlices each. First stack should be from the colored images to be
%%% analyzed. It wirdks with one cell at a time. 
%%% It also requires two zip files from the ImageJ software using plug-in 'ROI manager'
%%% 1) Nodes ROI: Save x-y-z positions of nodes and the last one is polynomial area of the broad band of nodes from SUM projected images
%%% 2) background ROI: Save x-y-z positions of 8 circles from each tip and coordinates of cell length, cell width and cell bundaries from SUM projected image in the order as shown
%%% in Figure 2

%%% (To get the ROIs from ImageJ and find out the pixel intensities of the circle 13 pixel diamter from the sum images of 7
%%% zstacks where all the fluorescenc of a unitary nodes resides. Careful observation showed that the Matlab and ImageJ handles image differently Eg. pixel intensity of (94, 107) in imageJ will be
%%% smilar to that of (95 108) in the matlab image.) 

close all;
clear all;
% % [Image,PathName] = uigetfile('*.czi', 'Select image');
resultfolder='RESULTS\';
IMAGE= bfopen (' Open the image');
Image= IMAGE{1,1}(1:2:end,:);
[NodeROI,PathName] = uigetfile('*.zip', 'Select Nodes ROI');
Node=ReadImageJROI(NodeROI);
x=[];
y=[];
z=[];
[BackGroundROI,PathName] = uigetfile('*.zip', 'Select background ROI');
BackGroundROI=ReadImageJROI(BackGroundROI);
X=[];
Y=[];
%% getting all the cooridnates from the NodeROI
for i=1:1:size(Node,2)
   tf = strcmp(Node{1,i}.strType,'Oval'); % compairing if it is node or node broad band
     if tf 'true';
    A=Node{1,i}.strName;
    y=[y;str2num(A(6:9))+1];
    x=[x;str2num(A(11:14))+1];
    z=[z;str2num(A(1:4))];
    else
        A=Node{1,i}.strName;
    y=[y;str2num(A(1:4))+1];
    x=[x;str2num(A(6:9))+1];
    end
end
%% getting all the cooridnates from the BackgroundROI

for j=1:1:size(BackGroundROI,2)
    B=BackGroundROI{1,j}.strName;
    Y=[Y;str2num(B(1:4))+1];
    X=[X;str2num(B(6:9))+1];
end

%% Getting 7 images corresponing to z coordinates, sum them  and calculating pixel intensity in the 12 pixel diamter circle.
cmap = gray(256);A = {};NodeIntensity=[];AvgBackgroundIntensity=[];   BgremovedNode=[]; Background=[];backGround=[];NodeNumber=[];AvgBGInt=[];
    firstFrame=[];lastFrame=[]; NodeInt= [];TotalInt=[];
for k = 1:1:size (z)    %Go into each of the ROI 
      NodeNumber=[NodeNumber;k];
      fF = z (k)- 3; %first Frame
      firstFrame= [firstFrame;fF];
      lF  = z(k) + 3;
      lastFrame = [lastFrame;lF];
          for l = fF:1:lF
              A(l)= (Image(l, 1));
          end
              sumImage= sum(cat(3, A{:}), 3);
              figure ;
              h=imshow(sumImage,cmap);
      
              xc= x(k);yc= y(k);
              [imageSizeY imageSizeX] = size(sumImage);    
              [columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);  
              centerx = x(k);
              centery = y(k);
              radius = 6;
              mask = (rowsInImage - centery).^2 + (columnsInImage - centerx).^2 <= radius.^2;
% % %               figure;
% % %               bw=imshow(mask) ;
% % %               colormap([0 0 0; 1 1 1]);
% % %               title('Binary image of a circle');
              viscircles([centerx centery], radius,'color', ' r' )
              NodeIntensity =sum(sum(double(sumImage) .* mask)) ;
              NodeInt=[NodeInt;NodeIntensity];
               A={};
%% background calculation              
              for K = 1:1:size(X,1)
                  Tf = strcmp(BackGroundROI{1,K}.strType,'Oval'); % comparing if it is node or node broad band
                  if Tf 'true'; 
                    centerX = X(K);
                    centerY = Y(K);
                    Mask = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
                    viscircles([centerX centerY], radius,'color', ' w' )
                    Background= [Background;sum(sum(double(sumImage) .* Mask))];
                  end
% % %                 BackgroundIntensity=sum([BackgroundIntensity;sum(sum(double(sumImage) .* mask))])./16;
              end
              AvgBGInt = sum(Background)/16; %Average Background intensity
              backGround= [backGround,Background];% To save background
              BgremovedNode=[BgremovedNode;NodeIntensity-AvgBGInt];
              Background=[]; %Everytime background is also changing
              AvgBackgroundIntensity= [AvgBackgroundIntensity;AvgBGInt];
             
    end
     
   
%% Get the total intensity of the node broad band area
         for m = 1:1:40
             BB(m)= (Image (m, 1));
          end
              sumImage= sum(cat(3, BB{:}), 3);
              figure ;
              J=imshow(sumImage,cmap);
                 XX =Node{1,end}.mnCoordinates (:,1)+1;
                 YY =Node{1,end}.mnCoordinates (:,2)+1;
                 hold on;
                 plot (XX,YY);
                 BW = poly2mask(XX,YY,imageSizeY,imageSizeX);
                 BroadbandInt=sum(sum(sumImage.*BW));
                 BroadbandArea=polyarea(XX, YY);
                 TipArea=[];meanTipArea= [];TipInt=[];meanTipInt=[];
                 
                 for n=3:4
                   XXB= BackGroundROI{1,end-n}.mnCoordinates (:,1)+1;
                   YYB= BackGroundROI{1,end-n}.mnCoordinates (:,2)+1;
                   BBW = poly2mask(XXB,YYB,imageSizeY,imageSizeX);
                   TipInt=[TipInt;sum(sum(sumImage.*BBW))/sum(sum (BBW))]; %%avg tip intensity=sum of pixel intensity/total no. of pixels
                   hold on;
                   plot (XXB,YYB);
                 end
                    
       meanTipInt=sum(TipInt)/2;
    NdBdInt = BroadbandInt-(meanTipInt*BroadbandArea);
    NumberofNode = NdBdInt / mean(BgremovedNode(BgremovedNode<=10000));
    %%%Total intensity of cell background subtracted
                   XXB= BackGroundROI{1,end}.mnCoordinates (:,1)+1;
                   YYB= BackGroundROI{1,end}.mnCoordinates (:,2)+1;
                   BBW = poly2mask(XXB,YYB,imageSizeY,imageSizeX);
                   TotalInt= sum(sum(sumImage.*BBW))-(meanTipInt*(polyarea(XXB, YYB)));
                   plot (XXB,YYB);
    
     %% Get the size of the cell
            
                a1= BackGroundROI{1,19}.vnLinePoints (1,1);
                b1= BackGroundROI{1,19}.vnLinePoints (2,1);
                a2 = BackGroundROI{1,19}.vnLinePoints (3,1);
                b2= BackGroundROI{1,19}.vnLinePoints (4,1);
            L =sqrt((abs(a2-a1))^2+(abs(b2-b1))^2)*0.041 ; %%Length of the cell
            
             A1= BackGroundROI{1,20}.vnLinePoints (1,1);
                B1= BackGroundROI{1,20}.vnLinePoints (2,1);
                A2 = BackGroundROI{1,20}.vnLinePoints (3,1);
                B2= BackGroundROI{1,20}.vnLinePoints (4,1);
             D = sqrt((abs(A2-A1))^2+(abs(B2-B1))^2)*0.041;  %% Daimeter of the cells
             
             V = (pi * (D /2)^2 *L) +  ((4/3)*pi*(D/2)^3)*0.041;   %%Volume of the cell
             S = (2 *pi*(D/2)*L) + (4*pi* (D/2)^2)*0.041 ;        %% Surface area of the cell
       
%% Write everything to excel.
             o=max(numel(z)); %To make equal size vector
             NdBdInt(2:o, 1)=nan;
             NumberofNode(2:o, 1)=nan;
             L(2:o, 1)=nan;
             D(2:o, 1)=nan;
             
             V(2:o, 1)=nan;
             S(2:o, 1)=nan;
             TotalInt(2:o, 1)=nan;
T = table(NodeNumber,firstFrame,lastFrame, NodeInt,AvgBackgroundIntensity,BgremovedNode, NdBdInt,NumberofNode, L, D, V, S, TotalInt);
    filename = 'Image 12.xlsx';
writetable(T,filename,'Sheet',1,'Range','D1')         

