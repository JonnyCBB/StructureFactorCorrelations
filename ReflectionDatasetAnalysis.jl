module ReflectionDatasetAnalysis

export
    summedRefProps,
    resBinBoxPlot,
    createDatasetInfoFile,
    writeDatasetInformation,
    createDirectory,
    setupMoleculesAndReflections,
    extractDatasetInfo!,
    addIntensityDiffColumn!,
    generateRandomTransformParams,
    chooseTransform,
    performMoleculeTrans,
    getRefDistancesAndCorrelations,
    plotCorrelationGraphs,
    writeCorrelationInfoToFile


using ReciprocalSpace
using ElementInformation
using ParsePDB
using DataFrames
using Gadfly
using LsqFit
using Color

#############################################################
#####        FUNCTIONS                                  #####
#############################################################

function performMoleculeTrans(pdbCode::SubString{ASCIIString}, waveLength::Float64, hklList::Dict{Vector{Int64},Reflection}, unitcellMols::MoleculeUnitCell,
                              resbins::Array{ResBin,1}, elementDict::Dict{ASCIIString,Element}, unitcell::Unitcell, numOfTrans::Int64,
                              transformType::ASCIIString, progressiveTrans::Bool, boxPlots::Vector{Int64}, refPlots::Vector{Int64}, numOfRefPlots::Int64,
                              refToPlot, refsToCompare)

    intValsforRefComparison = Array(Float64, numOfTrans, 2)

    datasetNum = 1
    hklOrig = calcScatFactors(hklList,unitcellMols,waveLength)
    extractDatasetInfo!(resbins, hklOrig, pdbCode, datasetNum)
    resBinBoxPlot(pdbCode, datasetNum, resbins, hklOrig)
    if numOfRefPlots ≥ 1
        for ref in refToPlot
            plotReflectionStrucFac(pdbCode, datasetNum, hklOrig[ref], unitcellMols, elementDict, false)
        end
    end

    hklPrev = hklOrig
    intDiffMat = Array(Float64,length(hklList),numOfTrans)
    for i = 1 : numOfTrans
        datasetNum += 1
        molsToTransform, serialNum, addAtomPosition, atomToAdd, transAtom, transMolecule, rotationAngle = generateRandomTransformParams(unitcellMols, unitcell)
        newAllMols, transformInfo = chooseTransform(transformType, unitcellMols, molsToTransform, serialNum, addAtomPosition, atomToAdd, transAtom, transMolecule,
                                                    rotationAngle, elementDict)
        hklChanged = calcScatFactors(hklList,newAllMols,waveLength)
        extractDatasetInfo!(resbins, hklChanged, pdbCode, datasetNum, transformInfo)
        addIntensityDiffColumn!(intDiffMat, hklPrev, hklChanged, i)
        if progressiveTrans
            unitcellMols = newAllMols
        end

        if any(boxPlots .== i) && i != 1
            resBinBoxPlot(pdbCode, datasetNum, resbins, hklChanged)
        end
        if any(refPlots .== i) && i != 1
            for ref in refToPlot
                plotReflectionStrucFac(pdbCode, datasetNum, hklChanged[ref], unitcellMols, elementDict, false)
            end
        end
        if length(refsToCompare) == 2
            intValsforRefComparison[i, 1] = hklPrev[refsToCompare[1]].squaredAmplitude - hklChanged[refsToCompare[1]].squaredAmplitude
            intValsforRefComparison[i, 2] = hklPrev[refsToCompare[2]].squaredAmplitude - hklChanged[refsToCompare[2]].squaredAmplitude
        end
        hklPrev = deepcopy(hklChanged)
    end
    if length(refsToCompare) == 2
        model(x, p) = p[1] * x + p[2]
        xdata = intValsforRefComparison[:,1]
        ydata = intValsforRefComparison[:,2]
        correlation = cor(xdata,ydata)
        fit = curve_fit(model, xdata, ydata, [1.0, 0.0])
        x = linspace(minimum(xdata), maximum(xdata), 100)
        y = fit.param[1] * x + fit.param[2]
        n = 2
        getColors = distinguishable_colors(n, ColorValue[LCHab(70, 60, 240)],
                       transform=c -> deuteranopic(c, 0.5),
                       lchoices=Float64[65, 70, 75, 80],
                       cchoices=Float64[0, 50, 60, 70],
                       hchoices=linspace(0, 330, 24))
        plt = plot(layer(x=xdata, y=ydata, Geom.point, Geom.line, Theme(default_color=color(getColors[1]))),
                   layer(x=x, y=y, Geom.line, Theme(default_color=color(getColors[2]))),
                   Guide.xlabel(@sprintf("Intensity differences for reflection %s",join(refsToCompare[1],","))),
                   Guide.ylabel(@sprintf("Intensity differences for reflection %s",join(refsToCompare[2],","))),
                   Guide.title(@sprintf("Correlation is: %.3f", correlation)),
                   Guide.manual_color_key("Colour Key",["Data", "Line of best fit"],[getColors[1],getColors[2]]))
        draw(PDF(@sprintf("%s/IntensityDiffComp_%s_and_%s.pdf",pdbCode,join(refsToCompare[1],"-"), join(refsToCompare[2],"-")), 8inch, 6inch), plt)
    end
    return intDiffMat, unitcellMols
end

function createDirectory(pdbFileName::ASCIIString)
    splitStrings = split(pdbFileName,".")
    if length(splitStrings) ≥ 2
        pdbCode = splitStrings[end - 1]
    else
        pdbCode = splitStrings[1]
    end

    if !isdir(pdbCode)
        mkdir(pdbCode)
    end
    createDatasetInfoFile(pdbCode)
    return pdbCode
end

function setupMoleculesAndReflections(pdbFileName::ASCIIString, xrayWavelength::Float64, resolutionLimit::Float64, numberOfResBins::Int64)
    elementDict = createElementDictionary()
    allMoleculesInUnitCell, molAsymUnit, unitcell, spacegroup, ncs, pdbInfo = generateUnitCell(pdbFileName, elementDict)
    reflectionList = getAllReflections(resolutionLimit,spacegroup,unitcell,xrayWavelength)
    allAngles = getAllScatteringAngles(reflectionList)
    minRefinResBin = int(round(length(reflectionList)/numberOfResBins))
    resolutionBins = sortResolutionBins(reflectionList, minRefinResBin)
    calcf0!(elementDict, allAngles, xrayWavelength)
    updateAtomicElements!(allMoleculesInUnitCell,elementDict)
    sortHKLIntoResBins!(resolutionBins,reflectionList)
    return reflectionList, allMoleculesInUnitCell, resolutionBins, elementDict, pdbInfo, unitcell
end

function generateRandomTransformParams(allMoleculesInUnitCell::MoleculeUnitCell, unitcell::Unitcell)
    numberOfMolsInUnitCell = length(allMoleculesInUnitCell.moleculeSequence)
    numberOfAtomsInMolecule = length(allMoleculesInUnitCell.moleculeSequence[1].atomSequence)
    numOfMolsToTransform = rand(1:numberOfMolsInUnitCell)
    molsToTransform = unique(vec(rand(1:numberOfMolsInUnitCell, numOfMolsToTransform, 1)))

    rotationAngle = float(rand(0:400)/10000)
    numMissingAtoms = allMoleculesInUnitCell.moleculeSequence[1].atomSequence[end].serialNum - numberOfAtomsInMolecule
    atom = allMoleculesInUnitCell.moleculeSequence[1].atomSequence[rand(1:numberOfAtomsInMolecule)]
    serialNum = atom.serialNum - numMissingAtoms ≤ 0 ? atom.serialNum : atom.serialNum - numMissingAtoms
    addAtomPosition = [float(rand(0:unitcell.a)), float(rand(0:unitcell.b)), float(rand(0:unitcell.c))]

    typicalElements = ["H", "O", "C", "N", "P", "S", "Se"]
    atomToAdd = typicalElements[rand(1:length(typicalElements))]

    meanDisplacement = sqrt(atom.bfactor/(8π^2))
    threeRandNums = rand(3,1)
    transAtom = 0.1*vec(sqrt(meanDisplacement^2 * threeRandNums/sum(threeRandNums)))

    #transMolecule = [0.005, 0, 0]
    transMolecule = transAtom

    return molsToTransform, serialNum, addAtomPosition, atomToAdd, transAtom, transMolecule, rotationAngle
end

function chooseTransform(transType::ASCIIString, molUnitCell::MoleculeUnitCell, moleculeNums::Vector{Int64}, serialID::Int64, position::Vector{Float64},
                         newElement::ASCIIString, translateAtom::Vector{Float64}, translateMolecule::Vector{Float64}, rotation::Float64,
                         elementDict::Dict{ASCIIString, Element})
    if transType == "move"
        transNum = rand(1:3)
        #transNum = 2
        if transNum == 1
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, serialID, translateAtom) #translate atom
        elseif transNum == 2
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, translateMolecule) #translate molecule
        else
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, rotation) #rotate molecule
        end
    elseif transType == "damage"
        transNum = rand(1:4)
        if transNum == 1
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, serialID, translateAtom) #translate atom
        elseif transNum == 2
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, translateMolecule) #translate molecule
        elseif transNum == 3
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, rotation) #rotate molecule
        else
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, serialID) #delete atom
        end
    else
        transNum = rand(1:6)
        if transNum == 1
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, serialID, translateAtom) #translate atom
        elseif transNum == 2
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, translateMolecule) #translate molecule
        elseif transNum == 3
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, rotation) #rotate molecule
        elseif transNum == 4
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, serialID) #delete atom
        elseif transNum == 5
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, moleculeNums, elementDict, serialID, newElement) #change atom
        elseif transNum == 6
            newMolUnitCell, transformInfo = transformMoleculeInUnitCell(molUnitCell, elementDict, position, newElement) #add atom
        end
    end
    return newMolUnitCell, transformInfo
end

function summedRefProps(reflectionList::Dict{Vector{Int64},Reflection})
    summedAmplitude = 0.0
    summedSqrdAmplitude = 0.0
    summedPhase = 0.0
    for key in keys(reflectionList)
        reflection = reflectionList[key]
        summedAmplitude += reflection.amplitude
        summedSqrdAmplitude += reflection.squaredAmplitude
        summedPhase += reflection.phase
    end
    modSummedPhase = mod(summedPhase,360)
    totalMeanIntensity = summedSqrdAmplitude/length(reflectionList)

    return summedAmplitude, summedSqrdAmplitude, summedPhase, modSummedPhase, totalMeanIntensity
end

function getRefDistancesAndCorrelations(intDiffMat::Array{Float64,2}, reflectionList::Dict{Vector{Int64},Reflection}, unitcell::Unitcell)
    corMat = cor(intDiffMat, vardim=2)
    n, m = size(corMat)
    numOfelements = int(round((n^2 - n)/2))
    corArray = Array(Float64, numOfelements,1)
    diffRes = similar(corArray)
    diffDist = similar(corArray)
    hkls = collect(keys(reflectionList))
    counter = 0
    recipUnitcell = convertDirectToReciprocal(unitcell)
    for col = 1:m
        for row = col+1:n
            if reflectionList[hkls[row]].resolution != Inf && reflectionList[hkls[col]].resolution != Inf
                counter += 1
                corArray[counter] = abs(corMat[row, col])
                ref1 = reflectionList[hkls[row]]
                ref2 = reflectionList[hkls[col]]
                diffRes[counter] = abs(ref1.resolution - ref2.resolution)
                point1 = [ref1.hkl[1] * recipUnitcell.a, ref1.hkl[2] * recipUnitcell.b, ref1.hkl[3] * recipUnitcell.c]
                point2 = [ref2.hkl[1] * recipUnitcell.a, ref2.hkl[2] * recipUnitcell.b, ref2.hkl[3] * recipUnitcell.c]
                diffDist[counter] = norm(point1 - point2)
            end
        end
    end
    #remove NaN values
    diffRes = diffRes[!isnan(corArray)]
    diffDist = diffDist[!isnan(corArray)]
    corArray = corArray[!isnan(corArray)]
    return corArray, diffDist, diffRes
end

function sortDistDiffBins(distanceArray::Vector{Float64}, correlationArray::Vector{Float64}, numOfDistsInBin::Int64)
    diffArray = [distanceArray correlationArray]
    diffArray = sortrows(diffArray)

    loopCounter = 1
    averageArray = Float64[]
    binDistArray = Vector{Float64}[]
    prevDist = diffArray[1,1]
    while loopCounter ≤ length(diffArray[:,2])
        totalCor = 0.0
        distCounter = 0
        minDist = prevDist
        while distCounter ≤ numOfDistsInBin || prevDist == diffArray[loopCounter,1]
            totalCor += diffArray[loopCounter, 2]
            prevDist = diffArray[loopCounter,1]
            distCounter += 1
            loopCounter += 1
            if loopCounter > length(diffArray[:,2])
                break
            end
        end
        push!(averageArray, totalCor/distCounter)
        push!(binDistArray,[minDist, prevDist, distCounter])
    end
    return averageArray, binDistArray
end

function plotCorrelationGraphs(pdbCode::SubString{ASCIIString}, corArray::Vector{Float64}, diffDist::Vector{Float64}, diffRes::Vector{Float64},
                               minNumPlotPoints::Int64, numberOfResBins::Int64, plotCorrelationBins::Bool=false)

    numOfelements = length(corArray)
    indicesToPlot = int(round(linspace(1,numOfelements,min(minNumPlotPoints,numOfelements))))
    if plotCorrelationBins
        n=2
        getColors = distinguishable_colors(n, ColorValue[LCHab(70, 60, 240)],
                                           transform=c -> deuteranopic(c, 0.5),
                                           lchoices=Float64[65, 70, 75, 80],
                                           cchoices=Float64[0, 50, 60, 70],
                                           hchoices=linspace(0, 330, 24))
        numberOfDistsinBin = int(round(numOfelements/numberOfResBins))
        averageArray, binDistArray = sortDistDiffBins(diffDist, corArray, numberOfDistsinBin)

        midAbsDistArray = Array(Float64, length(binDistArray), 1)
        counter = 0
        for vector in binDistArray
            counter += 1
            midAbsDistArray[counter] = (vector[2] + vector[1])/2
        end
        plt1 = plot(layer(x=diffDist[indicesToPlot],y=corArray[indicesToPlot],Geom.point, Theme(default_color=color(getColors[1])), order=1),
                    layer(x=midAbsDistArray[1:end],y=averageArray[1:end],Geom.point, Geom.line, Theme(default_color=color(getColors[2])), order=2),
                    Guide.xlabel("Distance between Reflections (Angstroms^-1)"), Guide.ylabel("Correlation"),
                    Guide.manual_color_key("Colour Key",["Individual Correlations", "Bin Averaged Correlations"],[getColors[1],getColors[2]]))

        averageArray, binDistArray = sortDistDiffBins(diffRes, corArray, numberOfDistsinBin)

        midResDistArray = Array(Float64, length(binDistArray), 1)
        counter = 0
        for vector in binDistArray
            counter += 1
            midResDistArray[counter] = (vector[2] + vector[1])/2
        end
        plt2 = plot(layer(x=diffRes[indicesToPlot],y=corArray[indicesToPlot],Geom.point, Theme(default_color=color(getColors[1])), order=1),
                    layer(x=midResDistArray[1:end],y=averageArray[1:end],Geom.point, Geom.line, Theme(default_color=color(getColors[2])), order=2),
                    Guide.xlabel("Difference in Reflection Resolution (Angstroms)"), Guide.ylabel("Correlation"),
                    Guide.manual_color_key("Colour Key",["Individual Correlations", "Bin Averaged Correlations"],[getColors[1],getColors[2]]))
        draw(PDF(@sprintf("%s/DiffResolutionCorrelation.pdf",pdbCode), 8inch, 6inch), plt2)
        draw(PDF(@sprintf("%s/DiffDistanceCorrelation.pdf",pdbCode), 8inch, 6inch), plt1)

    else
        plt1 = plot(x=diffDist[indicesToPlot],y=corArray[indicesToPlot],Geom.point,
                    Guide.xlabel("Difference in Reflection Resolution (Angstroms)"), Guide.ylabel("Correlation"))
        plt2 = plot(x=diffRes[indicesToPlot],y=corArray[indicesToPlot],Geom.point,
                    Guide.xlabel("Distance between Reflections (Angstroms^-1)"), Guide.ylabel("Correlation"))
        draw(PDF(@sprintf("%s/DiffResolutionCorrelation.pdf",pdbCode), 8inch, 6inch), plt2)
        draw(PDF(@sprintf("%s/DiffDistanceCorrelation.pdf",pdbCode), 8inch, 6inch), plt1)
    end
end

function resBinBoxPlot(dirPath::SubString{ASCIIString}, datasetNum::Int64, resolutionBins::Array{ResBin,1}, reflectionList::Dict{Vector{Int64},Reflection})
    cd(dirPath)
    csvfile = open("boxplotdata.csv","w")
    write(csvfile, join(["LogIntensity", "ResolutionRange"],","), "\n")
    for resbin in resolutionBins
        for key in resbin.hkls
            line = @sprintf("%f,%.2f-%.2f",log(reflectionList[key].squaredAmplitude), resbin.minRes, resbin.maxRes)
            write(csvfile,line,"\n")
        end
    end
    close(csvfile)
    df = readtable("boxplotdata.csv")
    plt = plot(df, x="ResolutionRange", y="LogIntensity", Geom.boxplot)
    draw(PDF(@sprintf("boxplot%d.pdf", datasetNum), 8inch, 6inch), plt)
    cd("..")
end

function createDatasetInfoFile(pdbCode::SubString{ASCIIString})
    cd(pdbCode)
    filename = @sprintf("%sdatasetInformation.txt",pdbCode)
    file = open(filename,"w")
    write(file,"#############################################################\n")
    write(file,"################# Dataset Information File ##################\n")
    write(file,"#############################################################\n")
    write(file,"\n\n")
    close(file)
    cd("..")
end

function writeDatasetInformation(pdbCode::SubString{ASCIIString}, datasetNum::Int64, sumAmp::Float64,
                                 sumSqrdAmp::Float64, sumPhase::Float64, modSumPhase::Float64,
                                 meanIntensity::Float64, resolutionBins::Array{ResBin,1}, transformInfo::ASCIIString)
    cd(pdbCode)
    filename = @sprintf("%sdatasetInformation.txt",pdbCode)
    file = open(filename,"a")
    write(file,"--------------------------------------------------\n")
    write(file,@sprintf("**************** Dataset Number %d ****************\n",datasetNum))
    write(file,"--------------------------------------------------\n")
    write(file,"\n")
    write(file,"********** Information about the transformation **********\n")
    write(file,"\n")
    write(file,@sprintf("Transformation - %s.\n", transformInfo))
    write(file,"\n")
    write(file,"********** Overall analysis of dataset **********\n")
    write(file,"\n")
    write(file,@sprintf("%-50s: %-12.3f\n", "Sum of all reflection amplitudes", sumAmp))
    write(file,@sprintf("%-50s: %-12.3f\n", "Sum of all reflection squared amplitudes", sumSqrdAmp))
    write(file,@sprintf("%-50s: %-12.3f\n", "Sum of all reflection phases (mod 360)", modSumPhase))
    write(file,@sprintf("%-50s: %-12.3f\n", "Average of all reflection squared amplitudes", meanIntensity))
    write(file,"\n")
    write(file,"********** Analysis of dataset by resolution bin **********\n")
    write(file,"\n")
    write(file,@sprintf("Total number of resolution bins: %d\n", length(resolutionBins)))
    write(file,"\n")
    write(file,@sprintf("%-2s%-22s%2s%-16s%2s%12s%2s%12s%2s%-12s%2s%12s%2s\n",
                        "|", "Resolution Range (Ang)", "|", "Number of Refs","|", "Mean Amp^2", "|", "Median Amp^2","|", "Std Dev","|","Std/Mean (\%)","|"))
    write(file,@sprintf("%-s%s%s%s%s%s%s%s%s%s%s%s%s\n",
                        "|", "------------------------", "|", "-----------------","|", "-------------", "|", "-------------","|", "-------------","|","-------------","|"))
    for resbin in resolutionBins
        write(file,@sprintf("%-2s%6.3f%3s%6.3f%7s%2s%16d%2s%12.3f%2s%12.3f%2s%12.3f%2s%12.3f%2s\n",
                            "|", resbin.minRes," - ", resbin.maxRes, "", "|",
                            resbin.numOfRef,"|", resbin.meanIntensity, "|", resbin.medianIntensity,"|", resbin.stdDev,"|",resbin.std2MeanRatio * 100,"|"))
    end
    write(file,"\n")
    write(file,"\n")
    close(file)
    cd("..")
end

function extractDatasetInfo!(resolutionBins::Array{ResBin,1}, reflectionList::Dict{Vector{Int64},Reflection}, pdbCode::SubString{ASCIIString},
                             datasetNum::Int64, transformInfo::ASCIIString="no transformation")
    sumAmp, sumSqrdAmp, sumPhase, modSumPhase, meanIntensity = summedRefProps(reflectionList)
    updateResBinInfo!(resolutionBins, reflectionList)
    writeDatasetInformation(pdbCode, datasetNum, sumAmp, sumSqrdAmp, sumPhase, modSumPhase, meanIntensity, resolutionBins, transformInfo)
end

function addIntensityDiffColumn!(intDiffMat::Array{Float64,2}, reflectionListOrig::Dict{Vector{Int64},Reflection},
                                 reflectionListChanged::Dict{Vector{Int64},Reflection}, transformationNum::Int64)
    intDiff = Array(Float64,length(reflectionListOrig),1)
    refCounter = 0
    for key in keys(reflectionListOrig)
        refCounter += 1
        intDiff[refCounter] = reflectionListOrig[key].squaredAmplitude - reflectionListChanged[key].squaredAmplitude
    end
    intDiffMat[:,transformationNum] = intDiff
end

function writeCorrelationInfoToFile(pdbCode::SubString{ASCIIString}, corArray::Vector{Float64}, diffDist::Vector{Float64}, diffRes::Vector{Float64})
    filename = @sprintf("%s/%sdatasetInformation.txt",pdbCode,pdbCode)
    file = open(filename,"a")
    write(file,"--------------------------------------------------\n")
    write(file,@sprintf("**************** Correlation Summary ****************\n"))
    write(file,"--------------------------------------------------\n")
    write(file,"\n")
    write(file,@sprintf("The correlation between the differences in the squared structure factor amplitudes
and the differences in the resolutions of the corresponding reflections is: %.3f\n\n", cor(corArray,diffRes)[1]))
    write(file,@sprintf("The correlation between the differences in the squared structure factor amplitudes
and the differences in (Euclidean) distances between the two corresponding reflections is: %.3f\n\n", cor(corArray,diffDist)[1]))
    close(file)
end


end
