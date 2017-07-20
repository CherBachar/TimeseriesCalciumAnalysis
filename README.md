# TimeseriesCalciumAnalysis
Software to analyse time-series neuron calcium signals from two-photon microscopy in-vivo

Requirements:

  Matlab + Image processing toolbox

Installation:

  Go into the folder, and add all folders to path

  To start: Run 'CalciumTimeSeriesAnalysisGUI.m' 

1. Press 'Load mean image'- load an image that's a mean/ max projection from the time-series stack OR an single image of the same location

2. Pres 'Segment' - Segmentation of cells, and detection of active cells will appear in the axes. Upload the corresponding time-series calcium image.

3. If you want to remove/ add cells, click on the corresponding button, and press 'Enter' once done. Please note you can add/remove both cells (left figure) and active cells (right figure). 
For active cells, draw a box around the cell you want to add.
To display cells that were edited, click 'Plot Cells'. 

4. Once you detected all the cells you want, press 'Analyse Ca2+' to analyse the time-series calcium signal. 

5. To view traces, you can change the cell in the editable textbox at the bottom of the trace

6. To save, press 'Save'. This will save all data in the same folder as 'Data.mat'

7. To load data, press 'Load Data' to load previously analysed data
