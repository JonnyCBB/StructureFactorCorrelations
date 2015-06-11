module ReciprocalSpace

export
    Reflection,
    ResBin,
    AllReflections,
    convertDirectToReciprocal,
    calcResolution,
    calcScatteringAngle,
    getAllScatteringAngles,
    determineHKLClass,
    applyRecipSpaceTrans,
    getAllReflections,
    calcScatFactors,
    plotReflectionStrucFac,
    removeSysAbs!,
    sortResolutionBins,
    sortHKLIntoResBins!,
    groupSymEquivRefs,
    updateResBinInfo!

using ParsePDB
using ElementInformation
using Gadfly
using Compose
using Color

#############################################################
#####        CREATE TYPES                               #####
#############################################################
type Reflection
    hkl::Vector{Int64}
    resolution::Float64
    scatteringAngle::Float64
    epsilon::Int64
    isCentric::Bool
    isSysAbs::Bool
    amplitude::Float64
    phase::Float64
    squaredAmplitude::Float64
    atomScatFactors::Array{Float64,2}

    function Reflection(hkl::Vector{Int64}, spacegroup::Spacegroup, unitcell::Unitcell, xrayWavelength::Float64)
        reflection = new()
        reflection.hkl = hkl
        recipspaceunitcell = convertDirectToReciprocal(unitcell)
        reflection.resolution = calcResolution(recipspaceunitcell,hkl)
        reflection.scatteringAngle = calcScatteringAngle(reflection.resolution, xrayWavelength)
        determineHKLClass!(reflection, spacegroup, unitcell)
        reflection.amplitude = 0.0
        reflection.phase = 0.0
        reflection.squaredAmplitude = 0.0
        reflection.atomScatFactors = Array(Float64,1,2)
        return reflection
    end
end

type ResBin
    minRes::Float64
    maxRes::Float64
    numOfRef::Int64
    totalIntensity::Float64
    medianIntensity::Float64
    meanIntensity::Float64
    stdDev::Float64
    std2MeanRatio::Float64
    hkls::Array{Vector{Int64},1}
end
ResBin(minRes::Float64, maxRes::Float64, numOfRef::Int64) = ResBin(minRes, maxRes, numOfRef, 0.0, 0.0, 0.0, 0.0, 0.0, Array(Vector{Int64},1))

#############################################################
#####        FUNCTIONS                                  #####
#############################################################

function convertDirectToReciprocal(directUnitCell::Unitcell)
    a = directUnitCell.a
    b = directUnitCell.b
    c = directUnitCell.c
    α = deg2rad(directUnitCell.α)
    β = deg2rad(directUnitCell.β)
    γ = deg2rad(directUnitCell.γ)

    #Calculate volume of direct Unit cell
    V = a*b*c*sqrt(1 - (cos(α))^2 - (cos(β))^2 - (cos(γ))^2 + (2 * cos(α) * cos(β) * cos(γ)))

    #Calculate reciprocal lattice parameters
    aᵣ = (    b * c * sin(α))/V
    bᵣ = (a *     c * sin(β))/V
    cᵣ = (a * b *     sin(γ))/V
    αᵣ = rad2deg(acos((cos(β)  * cos(γ) - cos(α)) / (sin(β) * sin(γ))));
    βᵣ = rad2deg(acos((cos(α)  * cos(γ) - cos(β)) / (sin(α) * sin(γ))));
    γᵣ = rad2deg(acos((cos(α)  * cos(β) - cos(γ)) / (sin(α) * sin(β))));

    #Create reciprocal unit cell lattice
    reciprocalUnitcell = Unitcell(aᵣ, bᵣ, cᵣ, αᵣ, βᵣ, γᵣ)
    return reciprocalUnitcell
end

function calcResolution(reciprocalUnitcell::Unitcell, millerIndices::Vector{Int64})
    h = millerIndices[1]
    k = millerIndices[2]
    l = millerIndices[3]

    a = reciprocalUnitcell.a
    b = reciprocalUnitcell.b
    c = reciprocalUnitcell.c
    α = deg2rad(reciprocalUnitcell.α)
    β = deg2rad(reciprocalUnitcell.β)
    γ = deg2rad(reciprocalUnitcell.γ)

    resolution = 1/sqrt(h^2 * a^2 + k^2 * b^2 + l^2 * c^2 +
                            2 * h * k * a * b * cos(γ) +
                            2 * h * l * a * c * cos(β) +
                            2 * k * l * b * c * cos(α))

    return resolution
end

function calcScatteringAngle(resolution::Float64, xrayWavelength::Float64)
    if isinf(resolution)
        return 0
    else
        d = resolution
        λ = xrayWavelength
        θ = rad2deg(asin(λ/2d)) #calculate scattering angle from Bragg's Law
        return θ
    end
end

function getAllScatteringAngles(reflectionList::Dict{Vector{Int64},Reflection})
    scatteringAngles = Array(Float64,length(reflectionList))
    counter = 1
    for key in keys(reflectionList)
        scatteringAngles[counter] = reflectionList[key].scatteringAngle
        counter += 1
    end
    scatteringAngles = unique(scatteringAngles)
    return scatteringAngles
end

function determineHKLClass!(reflection::Reflection, spacegroup::Spacegroup, unitcell::Unitcell)
    tolerance = 1e-3
    reflection.epsilon = 0
    reflection.isSysAbs = false
    reflection.isCentric = false
    for key in keys(spacegroup.symops)
        symEquivReflection, phaseShift = applyRecipSpaceTrans(spacegroup.symops[key],reflection.hkl,unitcell)
        if(symEquivReflection == reflection.hkl)
            if (cos(phaseShift) > (1 - tolerance))
                reflection.epsilon += 1
            else
                reflection.epsilon, reflection.isSysAbs = 0, true
                break
            end
        elseif(symEquivReflection == -reflection.hkl)
            reflection.isCentric = true
        end
    end
end

function applyRecipSpaceTrans(transformation::Array{Float64,2}, point::Vector{Int64}, unitcell::Unitcell)
    realSpaceMatrix = transformation[1:3,1:3]
    recipSpaceMatrix = transpose(realSpaceMatrix)
    translation = transformation[1:3,4]
    if sort(translation)[end] > 1.0
        relativeTranslation = [translation[1]/unitcell.a, translation[2]/unitcell.b, translation[3]/unitcell.c]
    else
        relativeTranslation = translation
    end
    pointTransformation = recipSpaceMatrix * point
    pointTransformation = convert(Array{Int64,1},round(pointTransformation))
    phaseShift = 2π*dot(relativeTranslation, point)
    return pointTransformation, phaseShift
end

function removeSysAbs!(reflectionList::Dict{Vector{Int64},Reflection})
    for key in keys(reflectionList)
        if reflectionList[key].isSysAbs
            delete!(reflectionList, key)
        end
    end
end

function getAllReflections(resolutionLimit::Float64, spacegroup::Spacegroup, unitcell::Unitcell, xrayWavelength::Float64)
    h, k, l = 0 ,0, 0
    reflectionList = Dict{Vector{Int64},Reflection}()
    currentReflectionResolution = Inf

    while currentReflectionResolution ≥ resolutionLimit
        reflection = Reflection([h,k,l], spacegroup, unitcell, xrayWavelength)
        currentReflectionResolution = reflection.resolution
        if currentReflectionResolution ≤ resolutionLimit
            break
        end

        while currentReflectionResolution ≥ resolutionLimit
            reflection = Reflection([h,k,l], spacegroup, unitcell, xrayWavelength)
            currentReflectionResolution = reflection.resolution
            if currentReflectionResolution ≤ resolutionLimit
                break
            end

            while currentReflectionResolution ≥ resolutionLimit
                reflection = Reflection([h,k,l], spacegroup, unitcell, xrayWavelength)
                currentReflectionResolution = reflection.resolution
                if currentReflectionResolution ≥ resolutionLimit
                    reflectionList[[ h, k, l]] = Reflection([ h, k, l], spacegroup, unitcell, xrayWavelength)
                    reflectionList[[ h, k,-l]] = Reflection([ h, k,-l], spacegroup, unitcell, xrayWavelength)
                    reflectionList[[ h,-k, l]] = Reflection([ h,-k, l], spacegroup, unitcell, xrayWavelength)
                    reflectionList[[-h, k, l]] = Reflection([-h, k, l], spacegroup, unitcell, xrayWavelength)
                    reflectionList[[ h,-k,-l]] = Reflection([ h,-k,-l], spacegroup, unitcell, xrayWavelength)
                    reflectionList[[-h, k,-l]] = Reflection([-h, k,-l], spacegroup, unitcell, xrayWavelength)
                    reflectionList[[-h,-k, l]] = Reflection([-h,-k, l], spacegroup, unitcell, xrayWavelength)
                    reflectionList[[-h,-k,-l]] = Reflection([-h,-k,-l], spacegroup, unitcell, xrayWavelength)
                    l += 1
                end
            end
            l = 0
            k += 1
            currentReflectionResolution = Inf
        end
        k = 0
        h += 1
        currentReflectionResolution = Inf
    end
    removeSysAbs!(reflectionList)
    return reflectionList
end

function calcScatFactors(reflectionList::Dict{Vector{Int64},Reflection}, moleculeArray::MoleculeUnitCell, xrayWavelength::Float64)
    newReflectionList = deepcopy(reflectionList)
    atomicf0::Float64 = 0.0
    atomicPhase::Float64 = 0.0
    tempFactor::Float64 = 0.0
    fᵢ::Complex64 = 0.0
    f::Complex64 = 0.0
    structureFactor::Complex64 = 0.0 + 0.0im

    totalNumAtoms = 0
    for molecule in moleculeArray.moleculeSequence
        totalNumAtoms += length(molecule.atomSequence)
    end

    for key in keys(newReflectionList)
        reflection = deepcopy(newReflectionList[key])
        reflection.atomScatFactors = Array(Float64,totalNumAtoms,2)
        structureFactor = 0.0 + 0.0im
        counter::Int64 = 1
        for molecule in moleculeArray.moleculeSequence
            for atom in molecule.atomSequence
                atomicf0 = atom.element.f0[reflection.scatteringAngle]
                atomicPhase = 2π*dot(reflection.hkl, [atom.x, atom.y, atom.z])
                reflection.atomScatFactors[counter,1] = atomicf0
                reflection.atomScatFactors[counter,2] = mod(rad2deg(atomicPhase), 360)
                fᵢ = atomicf0 * exp(im * atomicPhase)
                tempFactor = exp(-atom.bfactor * sin(deg2rad(reflection.scatteringAngle) / xrayWavelength)^2)
                f = atom.occupancy * fᵢ * tempFactor
                structureFactor += f
                counter += 1
            end
        end
        reflection.amplitude = abs(structureFactor)
        reflection.phase = mod(rad2deg(angle(structureFactor)), 360)
        reflection.squaredAmplitude = abs2(structureFactor)
        newReflectionList[key] = reflection
    end
    return newReflectionList
end

function calcArrowPoints(amplitude::Float64, phase::Float64, start_pos::Vector{Float64}, scaleArrowBase::Float64, halfArrowHeadBaseLength::Float64)
    f_x = amplitude * cos(phase)
    f_y = amplitude * sin(phase)

    end_pos = [start_pos[1] + f_x, start_pos[2] + f_y]
    arrowHeadBase = start_pos + (1-scaleArrowBase)*[f_x, f_y]
    perpendicularAngle = atan(-f_x/f_y)
    arrow_x = halfArrowHeadBaseLength * cos(perpendicularAngle)
    arrow_y = halfArrowHeadBaseLength * sin(perpendicularAngle)

    arrowPoint1 = arrowHeadBase + [arrow_x, arrow_y]
    arrowPoint2 = arrowHeadBase - [arrow_x, arrow_y]

    return arrowPoint1, arrowPoint2, end_pos, arrowHeadBase
end

function matchElementToScatFac(ref::Reflection, allMols::MoleculeUnitCell, totalNumAtoms::Int64, sort::Bool)
    elementArray = Array(ASCIIString,totalNumAtoms)
    counter = 0
    for molecule in allMols.moleculeSequence
        for atom in molecule.atomSequence
            counter += 1
            elementArray[counter] = atom.element.symbol
        end
    end
    fullArray = [elementArray ref.atomScatFactors]
    if sort
        fullArray = sortrows(fullArray)
    end
    return fullArray
end

function plotReflectionStrucFac(dirPath::SubString{ASCIIString}, datasetNum::Int64, ref::Reflection, allMols::MoleculeUnitCell,
                                elementDict::Dict{ASCIIString,Element}, sortedElements::Bool=false)
    cd(dirPath)
    h, k, l = ref.hkl[1], ref.hkl[2], ref.hkl[3]
    arrowScale = 0.3
    arrowHeadBaseScale = 0.4
    start_pos = [0.0, 0.0]
    maxAmplitude = maximum(ref.atomScatFactors[:,1])
    totalNumAtoms = length(ref.atomScatFactors[:,1])
    maxArrowWidth = 5
    gadflyPlotLayers = Array(Layer,totalNumAtoms + 1)
    arrowHeadAnnotations = Array(ASCIIString, totalNumAtoms + 3)
    arrowHeadAnnotations[1] = "Guide.annotation(compose(context()"
    arrowHeadAnnotations[end] = "))"

    fullArray = matchElementToScatFac(ref, allMols, totalNumAtoms, sortedElements)

    for atmIndex = 1 : totalNumAtoms
        amp = fullArray[atmIndex,2]
        phase = deg2rad(fullArray[atmIndex,3])
        halfArrowHeadBaseLength = (amp * arrowHeadBaseScale)/2.0 #Length of half of the base of the triangle
        arrowWidth = (amp/maxAmplitude) * maxArrowWidth #Decide the width of the arrow base on it's relative atomic scattering contribution
        point1, point2, end_pos, arrowHeadBase = calcArrowPoints(amp, phase, start_pos, arrowScale, halfArrowHeadBaseLength)

        gadflyPlotLayers[atmIndex] = layer(x=[start_pos[1],arrowHeadBase[1]], y=[start_pos[2],arrowHeadBase[2]], Geom.line, Theme(default_color=color(elementDict[fullArray[atmIndex]].color),line_width=eval(arrowWidth)px))[1]
        arrowHeadAnnotations[atmIndex + 1] = "(context(),polygon([($point1[1], $point1[2]), ($point2[1], $point2[2]), ($end_pos[1], $end_pos[2])]),fill($elementDict[$fullArray[$atmIndex]].color))"

        start_pos = deepcopy(end_pos)
    end

    #Calculate total structure factor
    halfArrowHeadBaseLength = (maxAmplitude * 2 * arrowHeadBaseScale)/2.0
    arrowWidth = maxArrowWidth * 2
    start_pos = [0.0, 0.0]
    point1, point2, end_pos, arrowHeadBase = calcArrowPoints(ref.amplitude, deg2rad(ref.phase), start_pos, arrowScale, halfArrowHeadBaseLength)
    gadflyPlotLayers[end] = layer(x=[start_pos[1],arrowHeadBase[1]], y=[start_pos[2],arrowHeadBase[2]], Geom.line, Theme(default_color=color("#C44E52"),line_width=eval(arrowWidth)px))[1]
    arrowHeadAnnotations[end-1] = "(context(),polygon([($point1[1], $point1[2]), ($point2[1], $point2[2]), ($end_pos[1], $end_pos[2])]),fill(RGB(196/255, 78/255, 82/255)))"

    #Produce plot
    plt = plot(gadflyPlotLayers, eval(parse(join(arrowHeadAnnotations,",",""))),
         Guide.title("Structure Factor for reflection $h, $k, $l"),
         Guide.xlabel(nothing), Guide.ylabel(nothing),
         Guide.manual_color_key("Key of element colours",
                                ["Hydrogen","Carbon","Nitrogen","Oxygen","Phosphorous","Sulphur","Selenium","Total Structure Factor","All other elements"],
                                [elementDict["H"].color, elementDict["C"].color, elementDict["N"].color, elementDict["O"].color, elementDict["P"].color,
                                 elementDict["S"].color, elementDict["Se"].color, RGB(196/255, 78/255, 82/255), elementDict["He"].color]),

         Coord.cartesian(fixed=true), Theme(panel_fill=color("#D3D3D3"), grid_color=color("#A8A8A8")))
    draw(PDF(@sprintf("reflectionplot%d_hkl-%s.pdf", datasetNum, join(ref.hkl,"-")), 8inch, 6inch), plt)
    cd("..")
end

function sortResolutionBins(reflectionList::Dict{Vector{Int64},Reflection}, minRefPerBin::Int64)
    resbins = ResBin[]
    counter = 1
    numOfRef = length(reflectionList)
    resolutionArray = Array(Float64,numOfRef)

    #First we need to sort the resolution values into increasing resolution order
    #(i.e. from biggest resolution number to smallest)
    for key in keys(reflectionList)
        resolutionArray[counter] = reflectionList[key].resolution
        counter += 1
    end
    resolutionArray = sort(resolutionArray, rev=true)

    #Now we can sort the resolution values
    previousRes = Inf
    numOfRefInResBinCounter = 0
    minRes = previousRes
    for i = 1:numOfRef
        if numOfRefInResBinCounter >= minRefPerBin && resolutionArray[i] != previousRes
            maxRes = previousRes
            resbin = ResBin(minRes, maxRes, numOfRefInResBinCounter)
            push!(resbins, resbin)
            minRes = previousRes
            numOfRefInResBinCounter = 0
        end
        numOfRefInResBinCounter += 1
        previousRes = resolutionArray[i]
    end
    return resbins
end

function sortHKLIntoResBins!(resbins::Array{ResBin,1}, reflectionList::Dict{Vector{Int64},Reflection})
    for resbin in resbins
        resbin.hkls = Array(Vector{Int64}, resbin.numOfRef)
        refCounter = 1
        for key in keys(reflectionList)
            reflection = reflectionList[key]
            if reflection.hkl == [0,0,0] && resbin.minRes == Inf
                resbin.hkls[refCounter] = reflection.hkl
                if refCounter == resbin.numOfRef
                    break
                else
                    refCounter += 1
                end
            elseif resbin.minRes > reflection.resolution ≥ resbin.maxRes
                resbin.hkls[refCounter] = reflection.hkl
                if refCounter == resbin.numOfRef
                    break
                else
                    refCounter += 1
                end
            end
        end
    end
end

function groupSymEquivRefs(reflectionList::Dict{Vector{Int64},Reflection}, spacegroup::Spacegroup, unitcell::Unitcell)
    symEquivs = Array{Int64,2}[]
    hklArray = collect(keys(reflectionList))
    while !isempty(hklArray)
        reflection = reflectionList[hklArray[1]]
        symEquivRefs = Array(Int64,2*length(spacegroup.symops),3)
        symEquivCounter = 0
        for key in keys(spacegroup.symops)
            symEquivReflection, phaseShift = applyRecipSpaceTrans(spacegroup.symops[key],reflection.hkl,unitcell)
            f(x) = x == symEquivReflection
            indexToDelete = find(f,hklArray)
            if !isempty(indexToDelete)
                symEquivCounter += 1
                symEquivRefs[symEquivCounter,:] = symEquivReflection
                splice!(hklArray, indexToDelete[1])
            end

            symEquivReflectionFriedel, phaseShift = applyRecipSpaceTrans(-spacegroup.symops[key],reflection.hkl,unitcell)
            f(x) = x == symEquivReflectionFriedel
            indexToDelete = find(f,hklArray)
            if !isempty(indexToDelete)
                symEquivCounter += 1
                symEquivRefs[symEquivCounter,:] = symEquivReflectionFriedel
                splice!(hklArray, indexToDelete[1])
            end
        end
        symEquivRefs = symEquivRefs[1:symEquivCounter,:]
        push!(symEquivs, symEquivRefs)
    end
    return symEquivs
end

function updateResBinInfo!(resbins::Array{ResBin,1}, reflectionList::Dict{Vector{Int64},Reflection})
    binCounter = 1
    for resbin in resbins
        intensityArray = Array(Float64,resbin.numOfRef,1)
        refCounter = 0
        for key in resbin.hkls
            refCounter += 1
            intensityArray[refCounter] = reflectionList[key].squaredAmplitude
        end
        resbin.totalIntensity = sum(intensityArray)
        resbin.medianIntensity = median(intensityArray)
        resbin.meanIntensity = mean(intensityArray)
        resbin.stdDev = std(intensityArray)
        resbin.std2MeanRatio = resbin.stdDev/resbin.meanIntensity
        resbins[binCounter] = resbin
        binCounter += 1
    end
end

end
