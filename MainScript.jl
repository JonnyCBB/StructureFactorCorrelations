using ElementInformation
using ParsePDB
using ReciprocalSpace
using ReflectionDatasetAnalysis
using Gadfly
using Debug
using DataFrames
using Color


############ Main script #############
# User defined params
pdbFileNames = ["test.pdb"]
λ = 1.0
resLim = 10.0
#minRefinResBin = 100
numberOfResBins = 12
numOfTrans = 30
transformType = "move"
progressiveTrans = true
numOfBoxPlots = 3
numOfRefPlots = 0
refToPlot = Array[[0,2,0]]
refsToCompare = Array[[0,2,0],[2,2,2]]
minNumPlotPoints = 2000
boxPlots = int(round(linspace(1,numOfTrans,numOfBoxPlots)))
refPlots = int(round(linspace(1,numOfTrans,numOfRefPlots)))
plotCorrelationBins = false

pdbFileName = pdbFileNames[1]
pdbCode = createDirectory(pdbFileName)
hklList, origMols, resbins, elementDict, pdbInfo, unitcell = setupMoleculesAndReflections(pdbFileName, λ, resLim, numberOfResBins)
intDiffMat, newMols = performMoleculeTrans(pdbCode, λ, hklList, origMols, resbins, elementDict, unitcell, numOfTrans, transformType,
                                           progressiveTrans, boxPlots, refPlots, numOfRefPlots, refToPlot, refsToCompare)
corArray, diffDist, diffRes = getRefDistancesAndCorrelations(intDiffMat, hklList, unitcell)
plotCorrelationGraphs(pdbCode, corArray, diffDist, diffRes, minNumPlotPoints, numberOfResBins, plotCorrelationBins)
writeCorrelationInfoToFile(pdbCode, corArray, diffDist, diffRes)

#resBinBoxPlot(pdbCode, datasetNum, resbins, hklChanged)

# ref1 = hklList[[3,1,-1]]
# ref2 = hklList[[1,-1,0]]
# plt1 = plotReflectionStrucFac(ref1, origMols, elementDict, false)
# plt2 = plotReflectionStrucFac(ref2, origMols, elementDict, false)
# draw(PDF("p1and2.pdf", 12inch, 6inch), hstack(plt1,plt2))

#writePDBFile(pdbCode, origMols, pdbInfo, "testpdbUnitCell.pdb")
#writePDBFile(pdbCode, unitcellMols, pdbInfo, "testpdbUnitCell2.pdb")
