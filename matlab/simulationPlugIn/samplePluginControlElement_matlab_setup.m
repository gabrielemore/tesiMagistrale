% To enable debug log functionality, declare "global debugLog" in the main script. 
function output = samplePluginControlElement_matlab_setup
output.requiresSqInsulinSupport = true;
output.requiresSqGlucoseSupport = false;
output.requiresInsulinDosingSupport = false;
output.requiresCPeptideSupport = false;
output.requiresSqGlucagonSupport = false;
output.numSignals = 1;
output.signalDescription = {'My CGM'};
output.numOutSignals = 0;
output.outSignalNames = {};
