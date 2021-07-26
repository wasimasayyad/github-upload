%%% Similar to the NodaAnalysis.m file it also requies 1 z-stack image with two coloured 40 z-slices and two zip 'ROI manager' fils from ImageJ
%%% except it can work with more than one cells in the same image. 
%%% 1)CR ROI: Save Outline of CR in the cell, length, Diamter and outline area of the cell in order 
%%% 2)CRBackgroundROI: Save first two tip intensities of the cell, then 3rd and fourth for second cell and so on. 
 

%% use ImageJ to segment CR as follows Take max projection_ USE Process--FFt--Bandpassfilter
%% then FFT to filter 
close all;
clear all;
% % [Image,PathName] = uigetfile('*.czi', 'Select image');
resultfolder='RESULTS\';
IMAGE= bfopen (' Open the image');
Image= IMAGE{1,1}(1:2:end,:);
[CRROI,PathName] = uigetfile('*.zip', 'Select CR ROI');
CR=ReadImageJROI(CRROI);
x=[];
y=[];
z=[];
[CRBackGroundROI,PathName] = uigetfile('*.zip', 'Select CRbackground ROI');
CRBackGroundROI=ReadImageJROI(CRBackGroundROI);
X=[];
Y=[];
%% Get the total intensity of the CR area
cmap = gray(256);
         for m = 1:1:40
             BB(m)= (Image (m, 1));
          end
              sumImage= sum(cat(3, BB{:}), 3);
              figure ;
              J=imshow(sumImage,cmap);
              [imageSizeY imageSizeX] = size(sumImage);    
              [columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
              CRInt = []; TipArea=[];meanTipArea= [];TipInt=[];meanTipInt=[];MeanTipInt=[]; L=[];D=[];, V=[]; S=[];TotalInt=[]; d=[];s=[];,v=[];l=[];
              for k = 1:1:size (CR,2)/4  %Go into each of the CR ROI
                 XX =CR{1,k*4-3}.mnCoordinates (:,1)+1; %%The CR ROI will be repaeated after every 4th position starting with first position in ROI
                 YY =CR{1,k*4-3}.mnCoordinates (:,2)+1;
                 hold on;
                 plot (XX,YY);
                 BW = poly2mask(XX,YY,imageSizeY,imageSizeX);
                 crInt=sum(sum(sumImage.*BW));
                 CRArea=polyarea(XX, YY);
                 
                 
                   %%%Get the background roi of tip1 and 2    
                   XXB1= CRBackGroundROI{1,k*2-1}.mnCoordinates (:,1)+1;
                   YYB1= CRBackGroundROI{1,k*2-1}.mnCoordinates (:,2)+1;
                   BBW1 = poly2mask(XXB1,YYB1,imageSizeY,imageSizeX);
                   
                   TipInt1=[TipInt;sum(sum(sumImage.*BBW1))/sum(sum (BBW1))]; %%avg tip intensity=sum of pixel intensity/total no. of pixels
                   hold on;
                   plot (XXB1,YYB1);
                       
                   XXB2= CRBackGroundROI{1,k*2}.mnCoordinates (:,1)+1;
                   YYB2= CRBackGroundROI{1,k*2}.mnCoordinates (:,2)+1;
                   BBW2 = poly2mask(XXB2,YYB2,imageSizeY,imageSizeX);
                   
                   TipInt2=[TipInt;sum(sum(sumImage.*BBW2))/sum(sum (BBW2))]; %%avg tip intensity=sum of pixel intensity/total no. of pixels
                   hold on;
                   plot (XXB2,YYB2);
% %                    meanTipInt=(TipInt1+TipInt2)/2;
                    meanTipInt = (TipInt1+TipInt2)/2;
                    MeanTipInt=[MeanTipInt;meanTipInt];
                   CRInt = [CRInt;crInt-(meanTipInt*CRArea)];
                   
                   
                   %%%Get the size of the cell
                a1= CR{1,k.*4-2}.vnLinePoints (1,1); %%The length ROI will be repaeated after every 4th position starting with second position in ROI
                b1= CR{1,k.*4-2}.vnLinePoints (2,1);
                a2 = CR{1,k.*4-2}.vnLinePoints (3,1);
                b2= CR{1,k.*4-2}.vnLinePoints (4,1);
            l =sqrt((abs(a2-a1)).^2+(abs(b2-b1)).^2).*0.041 ; %%Length of the cell
            L=[L;l];
            
             A1= CR{1,k.*4-1}.vnLinePoints (1,1);  %%The dimeter ROI will be repaeated after every 4th position starting with third position in ROI
                B1= CR{1,k.*4-1}.vnLinePoints (2,1);
                A2 = CR{1,k.*4-1}.vnLinePoints (3,1);
                B2= CR{1,k.*4-1}.vnLinePoints (4,1);
             d = sqrt((abs(A2-A1)).^2+(abs(B2-B1)).^2).*0.041;  %% Daimeter of the cells
             D=[D;d];
             v = (pi.*(d./2).^2.*l) +  ((4/3)*pi.*(d/2).^3).*0.041;   %%Volume of the cell
             V=[V;v];
             s = (2. *pi.*(d/2).*l) + (4.*pi.* (d/2).^2).*0.041 ;        %% Surface area of the cell
             S=[S;s];
             
                    %%%Total intensity of cell background subtracted
                   XX= CR{1,k*4}.mnCoordinates (:,1)+1; %%The totalInt ROI will be repaeated after every 4th position starting with fourth position in ROI
                   YY= CR{1,k*4}.mnCoordinates (:,2)+1;
                  BBW = poly2mask(XX,YY,imageSizeY,imageSizeX);
                   TotalInt= [TotalInt;sum(sum(sumImage.*BBW))-(meanTipInt*(polyarea(XX, YY)))];
                   plot (XX,YY);
                   
              end
                 
              %% Write everything to excel.
% % %               o=k; %To make equal size vector
% % %                  CRInt(2:o, 1)=nan;
% % %                   meanTipInt(2:o, 1)=nan;
% % %                   L(2:o, 1)=nan;
% % %                   D(2:o, 1)=nan;
% % %              
% % %              V(2:o, 1)=nan;
% % %              S(2:o, 1)=nan;
% % %              TotalInt(2:o, 1)=nan;
              T = table(CRInt,MeanTipInt, L, D, V, S, TotalInt);
    filename = 'Image 04.xlsx';
writetable(T,filename,'Sheet',1,'Range','D1')   
     