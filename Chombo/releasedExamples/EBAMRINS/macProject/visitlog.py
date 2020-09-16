# Visit 3.0.1 log file
ScriptVersion = "3.0.1"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
ShowAllWindows()
SetWindowLayout(1)
OpenDatabase("localhost:/home/spencer/users/graves/_svnchombo/releasedExamples/EBAMRINS/macProject/divu.hdf5", 0)
metadata = GetMetaData("localhost:/home/spencer/users/graves/_svnchombo/releasedExamples/EBAMRINS/macProject/divu.hdf5", -1)
AddPlot("Pseudocolor", "component_0", 1, 0)
AddOperator("ThreeSlice", 0)
DeleteActivePlots()
AddPlot("Pseudocolor", "component_0", 1, 0)
DrawPlots()
