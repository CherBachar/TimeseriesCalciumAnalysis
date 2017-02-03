# TimeseriesCalciumAnalysis
Software to analyse time-series neuron calcium signals from two-photon microscopy in-vivo

Requirements:

Matlab + Image processing toolbox

Steps for analysis:

Run 'CalciumTimeSeriesAnalysisGUI.m' - make sure all folders are in path.

1. Press 'Load mean image'- load an image that's a mean/ max projection from the time-series stack OR an single image of the same location

2. Pres 'Segment' - Point detection and segmentation will appear in the axes

3. If you want to remove/ add points, click on the corresponding botton, and press Enter once done.

4. Once you detected all the cells you want, press 'Analyse Ca2+' to analyse the time-series calcium signal. Upload the corresponding time-series calcium image.

5. To view traces, you can change the cell in the editable textbox at the bottom of the trace

6. To save, press 'Save'. This will save all data in the same folder as 'Data.mat'

*Add/ remove peaks functions don't work yet
