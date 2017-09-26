# TimeseriesCalciumAnalysis
Software to analyse time-series neuron calcium signals from two-photon microscopy in-vivo

Requirements:

  Matlab + Image processing toolbox

Installation:

  Go into the folder, and add all folders to path

  To start: Run 'CalciumTimeSeriesAnalysisGUI.m' 

1. To start press 'Load mean image'- load an image that's a mean/ max projection from the time-series stack OR an single image of the same location.
It will then request to load the corresponding time-series TIFF stack.

2. To analyse the time-series cells press 'Analyse Ca2+'. Here the active cells will be analysed (i.e. the cells detected on the right image).

You have options to edit both the Cells and Active Cells. The difference is that active cells are detected from the standard deviation image, to show whether it's been active or now.

These are the following options for editing and display:
- If you want to remove/ add cells, click on the corresponding button, and press 'Enter' once done. Please note you can add/remove both cells (left figure) and active cells (right figure). 
For active cells, draw a box around the cell you want to add.
To display cells that were edited, click 'Plot Cells'.

- To automatically detect Active cells from the Cells selected, press 'Detect Active cells' - It will automatically take the cells on the left as seed points and will look whether these cells are active, and will display on the right image.

- To view traces, you can change the cell in the editable textbox at the bottom of the trace

- For better view of cells in the left or right image, press 'Change Contrast', and select the appropriate values. 
DO NOT CLICK "Adjust Data" as this will override the pixel values. Click X on the top right to exit.

3. To save, press 'Save'. This will save all data in the same folder as 'filename.mat'

4. To load data, press 'Load Data' to load previously analysed data
