# github-upload
This project is about automatically extracting information from images using MatLab code and by providing the coordinates of ROI saved in ROI-manager.zip files obtained using 
ImageJ software. Please see the following manuscript for more detailes
https://www.biorxiv.org/content/10.1101/2021.07.02.450929v1

It has two programs. Specifically, the programes were developed to work on the images obtained by using Airyscan microscopy to extract the intensities of unitary nodes and cell 
size parameters They requires two colored Z-stack of Airyscan images containing 40 z-slices each. First stack should be from the colored images to be analyzed.

A) NodeBroadbandAnalysis.m
It works with one cell at a time in an image and exrtacts information from a broad band of nodes 
Along with the image files it also requires two zip files from the ImageJ software using plug-in 'ROI manager'
1) Nodes ROI: Containing x-y-z positions of nodes and the last one should be polynomial area of the broad band of nodes from SUM projected images
2) background ROI: Containging x-y-z positions of 8 circles from each tip and coordinates of cell length, cell width and cell bundaries from SUM projected image in the order as
shown in Figure 2 of the manuscript mentioned above. 

B) CRAnalysis.m
It can work with more than one cells in the same image and exrtacts information from a CR of nodes.
Similar to the NodaAnalysis.m file it also requies two zip 'ROI manager' files from ImageJ. 
1)CR ROI: Save Outline of CR in the cell, length, Diamter and outline area of the cell in order 
2)CRBackgroundROI: Save first two tip intensities of the cell, then 3rd and fourth for second cell and so on. 
