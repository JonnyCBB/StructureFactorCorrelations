module ParsePDB
using ElementInformation

export
    Atom,
    Unitcell,
    Spacegroup,
    NCS,
    Molecule,
    MoleculeUnitCell,
    parseCrystLine,
    parsePDBFile,
    generateUnitCell,
    createFullAsymUnit,
    generateMolSymmCopies,
    applyTransformation,
    updateAtomicElements!,
    writePDBFile,
    addAtom,
    deleteAtom,
    translateAtom,
    translateMolecule,
    changeAtom,
    rotateMolecule,
    transformMoleculeInUnitCell


#############################################################
#####        CREATE TYPES                               #####
#############################################################

#Create unit cell type with the standard unit cell definitions
#Unit cell lenghts a, b and c with the angles alpha, beta and
# gamma
type Unitcell{T<:Float64}
    a::T
    b::T
    c::T
    α::T
    β::T
    γ::T
end
#Define an empty constructor for the unit cell
Unitcell() = Unitcell(0.0,0.0,0.0,90.0,90.0,90.0)

type Spacegroup
    num_symops::Int64
    spgrp_symbol::ASCIIString
    symops::Dict{Int64, Array{Float64,2}}
end
#Define an empty constructor for the spacegroup
Spacegroup() = Spacegroup(0,"Empty constructor used",[1=>[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0]])

type NCS
    num_ncsops::Int64
    ncsops::Dict{Int64, Array{Float64,2}}
end
#Define an empty constructor for the NCS operators
NCS() = NCS(0,[1=>[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0]])

type Atom
    serialNum::Int64
    atomName::ASCIIString
    element::Element
    residueName::ASCIIString
    chainID::Char
    resSeqNum::Int64
    x::Float64
    y::Float64
    z::Float64
    occupancy::Float64
    bfactor::Float64
    recordName::ASCIIString
    altLocation::Char
    iCode::Char
    charge::ASCIIString
end
Atom(atom::Atom, position::Vector{Float64}) = Atom(atom.serialNum, atom.atomName, atom.element,
                                                   atom.residueName, atom.chainID, atom.resSeqNum,
                                                   position[1], position[2], position[3],
                                                   atom.occupancy, atom.bfactor, atom.recordName,
                                                   atom.altLocation, atom.iCode, atom.charge)

type Molecule
    pdbCode::String
    atomSequence::Array{Atom,1}
end
Molecule() = Molecule("N/A",Atom[])
Molecule(molecule::Molecule) = Molecule(molecule.pdbCode, molecule.atomSequence)

type MoleculeUnitCell
    moleculeSequence::Array{Molecule,1}
end
MoleculeUnitCell() = MoleculeUnitCell(Molecule[])

type PDBInformation
   header::ASCIIString
end
PDBInformation() = PDBInformation("")

#############################################################
#####        FUNCTIONS                                  #####
#############################################################

######################
#Parse the cryst1 line
######################
function parseCrystLine(line::ASCIIString)
    cryst1 = split(line)
    unit_cell = Unitcell(parsefloat(cryst1[2]),parsefloat(cryst1[3]),parsefloat(cryst1[4]),
                         parsefloat(cryst1[5]),parsefloat(cryst1[6]),parsefloat(cryst1[7]))
    return unit_cell
end

######################
#Parse the PDB file
######################
function parsePDBFile(pdbFileName::ASCIIString, elementDict::Dict{ASCIIString, Element})
    #initialise empty objects to be created
    unitcell = Unitcell()
    spacegroup = Spacegroup()
    ncs = NCS()
    moleculeAsymUnit = Molecule()
    pdbInfo = PDBInformation()

    #define variables that are required for use outside the scope of
    #"for loop"
    spgrp_symbol = ""
    sym_op_dict = Dict{Int64, Array{Float64,2}}()
    ncs_op_dict = Dict{Int64, Array{Float64,2}}()
    number_of_symops = 0
    number_of_ncsops = 0
    mtrx_row1 = [0.0 0.0 0.0 0.0]
    mtrx_row2 = [0.0 0.0 0.0 0.0]
    mtrx_row3 = [0.0 0.0 0.0 0.0]
    mtrx_row4 = [0.0 0.0 0.0 0.0]
    headerInfo = true

    #Extract the PDB code.
    if length(split(pdbFileName,".")) > 1
        PDBcode = split(pdbFileName,".")[1]
    else
        PDBcode = pdbFileName
    end
    moleculeAsymUnit.pdbCode = PDBcode

    #Open the PDB file stream for reading
    pdbFile = open(pdbFileName)

    #Read PDB file line by line
    for line in eachline(pdbFile)
        #Extract the unit cell information
        if line[1:6] == "CRYST1"
            unitcell = parseCrystLine(line)
            #Extract the crystallographic symmetry operators
        elseif contains(line,"REMARK 290")
            if contains(line,"SPACE GROUP")
                index = search(line,':')
                spgrp_symbol = join(split(line[index+1:end]))
            elseif contains(line,"SMTRY")
                splitLine = split(line)
                mtrx_op_row = [parsefloat(splitLine[end-3]) parsefloat(splitLine[end-2]) parsefloat(splitLine[end-1]) parsefloat(splitLine[end])]
                if contains(line,"SMTRY1")
                    mtrx_row1 = mtrx_op_row
                elseif contains(line,"SMTRY2")
                    mtrx_row2 = mtrx_op_row
                else
                    mtrx_row3 = mtrx_op_row
                    number_of_symops = parseint(splitLine[end-4])
                    sym_op_dict[number_of_symops] = [mtrx_row1; mtrx_row2; mtrx_row3; mtrx_row4]
                end
            end
            #Extract the Non-crystallographic symmetry (NCS) operators
        elseif line[1:5] == "MTRIX"
            splitLine = split(line)
            mtrx_op_row = [parsefloat(splitLine[3]) parsefloat(splitLine[4]) parsefloat(splitLine[5]) parsefloat(splitLine[6])]
            if contains(line,"MTRIX1")
                mtrx_row1 = mtrx_op_row
            elseif contains(line,"MTRIX2")
                mtrx_row2 = mtrx_op_row
            else
                mtrx_row3 = mtrx_op_row
                number_of_ncsops = parseint(splitLine[2])
                ncs_op_dict[number_of_ncsops] = [mtrx_row1; mtrx_row2; mtrx_row3; mtrx_row4]
            end
            #Extract the atom coordinate information
        elseif contains(line[1:6],"ATOM") || contains(line[1:6],"HETATM")
            if headerInfo == true
                headerInfo = false
            end
            serialID = parseint(line[7:11])
            atomName = strip(line[13:16])
            elementSymbol = strip(line[77:78])
            residueName = strip(line[18:20])
            chainID = line[22]
            residueSequenceNumber = parseint(line[23:26])
            xcoord = parsefloat(line[31:38])
            ycoord = parsefloat(line[39:46])
            zcoord = parsefloat(line[47:54])
            occupancy = parsefloat(line[55:60])
            bfactor = parsefloat(line[61:66])
            recordName = line[1:6]
            altLocation = line[17]
            iCode = line[27]
            charge = strip(line[79:80])
            if haskey(elementDict,elementSymbol)
                element = elementDict[elementSymbol]
            else
                error("Element symbol does not exist. Please check your element symbols!\n
                If symbols are ok then contact the Garman Group: elspeth.garman@bioch.ox.ac.uk\n
                Perhaps the Cromer-Mann coefficient table needs updating.\n")
            end
            atom = Atom(serialID, atomName, element, residueName, chainID, residueSequenceNumber,
                        xcoord, ycoord, zcoord, occupancy, bfactor, recordName, altLocation,
                        iCode, charge)
            push!(moleculeAsymUnit.atomSequence,atom)

        end
        #Keep the PDB header information for when the PDB file needs to be rewritten
        #with new atomic coordinates.
        if headerInfo == true
            pdbInfo.header = pdbInfo.header*line
        end
        spacegroup = Spacegroup(number_of_symops, spgrp_symbol, sym_op_dict)
        ncs = NCS(number_of_ncsops, ncs_op_dict)
    end
    return moleculeAsymUnit, unitcell, spacegroup, ncs, pdbInfo
end

##################################
#Generate a unit cell of molecules
##################################

function generateUnitCell(pdbFileName::ASCIIString, elementDict::Dict{ASCIIString,Element})
    moleculeUnit, unitcell, spacegroup, ncs, pdbInfo = parsePDBFile(pdbFileName, elementDict)
    createFullAsymUnit!(moleculeUnit,ncs)
    moleculeUnitCell = generateMolSymmCopies(moleculeUnit,spacegroup)
    return moleculeUnitCell, moleculeUnit, unitcell, spacegroup, ncs, pdbInfo
end

function createFullAsymUnit!(moleculeUnit::Molecule, ncs::NCS)
    if ncs.num_ncsops != 0
        numberOfAtoms = length(moleculeUnit.atomSequence)
        for key in keys(ncs.ncsops)
            ncsOperation = ncs.ncsops[key]
            if !((ncsOperation[1:3,1:3] == eye(3)) && (ncsOperation[1:3,4] == [0,0,0]))
                for atomIndex in 1:numberOfAtoms
                    atomPosition = [moleculeUnit.atomSequence[atomIndex].x,
                                    moleculeUnit.atomSequence[atomIndex].y,
                                    moleculeUnit.atomSequence[atomIndex].z]
                    symAtomPos = applyTransformation(ncsOperation, atomPosition)
                    symAtom = Atom(moleculeUnit.atomSequence[atomIndex], symAtomPos)
                    push!(moleculeUnit.atomSequence,symAtom)
                end
            end
        end
    end
end

function generateMolSymmCopies(molAsymUnit::Molecule, spacegroup::Spacegroup)
    allMoleculeCopies = MoleculeUnitCell()
    for key in keys(spacegroup.symops)
        symRelMol = Molecule()
        symRelMol.pdbCode = molAsymUnit.pdbCode;
        symop = spacegroup.symops[key]
        if !((symop[1:3,1:3] == eye(3)) && (symop[1:3,4] == [0,0,0]))
            for atom in molAsymUnit.atomSequence
                atomPosition = [atom.x, atom.y, atom.z]
                symAtomPos = applyTransformation(symop, atomPosition)
                symAtom = Atom(atom, symAtomPos)
                push!(symRelMol.atomSequence,symAtom)
            end
        else
            symRelMol = molAsymUnit
        end
        push!(allMoleculeCopies.moleculeSequence, symRelMol)
    end
    return allMoleculeCopies
end

function applyTransformation(transformation::Array{Float64,2},point::Vector{Float64})
    matrix = transformation[1:3,1:3]
    translation = transformation[1:3,4]
    pointTransformation = matrix * point + translation
    return pointTransformation
end

function updateAtomicElements!(moleculeArray::MoleculeUnitCell, elementDict::Dict{ASCIIString,Element})
    for molecule in moleculeArray.moleculeSequence
        for atom in molecule.atomSequence
            atom.element = deepcopy(elementDict[atom.element.symbol])
        end
    end
end

function writeAtomLine(atom::Atom)
    line = @sprintf("%6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%11s%2s%2s\n",
                    atom.recordName, atom.serialNum, "", atom.atomName, atom.altLocation, atom.residueName, "",
                    atom.chainID, atom.resSeqNum, atom.iCode, "",
                    atom.x, atom.y, atom.z, atom.occupancy, atom.bfactor,"",
                    atom.element.symbol, atom.charge)
end

function writePDBFile(dirPath::SubString{ASCIIString}, molecule::Molecule, pdbHeader::PDBInformation, outputFilename::ASCIIString)
    cd(dirPath)
    f = open(outputFilename,"w")
    write(f,pdbHeader.header)
    for atom in molecule.atomSequence
        atomLine = writeAtomLine(atom)
        write(f,atomLine)
    end
    close(f)
    cd("..")
end

function writePDBFile(dirPath::SubString{ASCIIString}, moleculeArray::MoleculeUnitCell, pdbHeader::PDBInformation, outputFilename::ASCIIString)
    cd(dirPath)
    f = open(outputFilename,"w")
    write(f,pdbHeader.header)
    for molecule in moleculeArray.moleculeSequence
        for atom in molecule.atomSequence
            atomLine = writeAtomLine(atom)
            write(f,atomLine)
        end
    end
    close(f)
    cd("..")
end

function addAtom(molecule::Molecule, elementDict::Dict{ASCIIString, Element}, position::Vector{Float64},
                  elementSymbol::ASCIIString, occupancy::Float64=1.00, bfactor::Float64=0.0, chainID::Char='A',
                  altLocation::Char='A', resSeqNum::Int64=2000, iCode::Char='?', charge::ASCIIString="",
                  recordName::ASCIIString="ATOM  ")
    newMolecule = deepcopy(molecule)
    element = elementDict[elementSymbol]
    serialNum = length(newMolecule.atomSequence) + 1
    atom = Atom(serialNum, elementSymbol, element, elementSymbol, chainID, resSeqNum,
                position[1], position[2], position[3], occupancy, bfactor,
                recordName, altLocation, iCode, charge)

    push!(newMolecule.atomSequence, atom)
    return newMolecule
end

function transformMoleculeInUnitCell(molUnitCell::MoleculeUnitCell, elementDict::Dict{ASCIIString, Element},
                                     position::Vector{Float64}, elementSymbol::ASCIIString, occupancy::Float64=1.0, bfactor::Float64=0.0,
                                     chainID::Char='A', altLocation::Char='A', resSeqNum::Int64=2000, iCode::Char='?', charge::ASCIIString="",
                                     recordName::ASCIIString="ATOM  ")
    newMolUnitCell = deepcopy(molUnitCell)
    molecule = newMolUnitCell.moleculeSequence[1]
    newMolecule = addAtom(molecule, elementDict, position, elementSymbol, occupancy, bfactor,
                          chainID, altLocation, resSeqNum, iCode, charge, recordName)
    newMolUnitCell.moleculeSequence[1] = newMolecule
    transformInfo = @sprintf("addition of atom: '%s' at position (x, y, z) = (%.2f, %.2f, %.2f) in molecule number 1",
                             elementSymbol, position[1], position[2], position[3])
    return newMolUnitCell, transformInfo
end

function deleteAtom(molecule::Molecule, serialID::Int64)
    newMolecule = deepcopy(molecule)
    for i = 1:length(newMolecule.atomSequence)
        if serialID == newMolecule.atomSequence[i].serialNum
            splice!(newMolecule.atomSequence, i)
            break
        end
    end
    return newMolecule
end

function transformMoleculeInUnitCell(molUnitCell::MoleculeUnitCell, moleculeNums::Vector{Int64}, serialID::Int64)
    newMolUnitCell = deepcopy(molUnitCell)
    atom = newMolUnitCell.moleculeSequence[1].atomSequence[serialID]
    for molNum in moleculeNums
        molecule = newMolUnitCell.moleculeSequence[molNum]
        newMolecule = deleteAtom(molecule, serialID)
        newMolUnitCell.moleculeSequence[molNum] = newMolecule
    end
    transformInfo = @sprintf("deletion of atom: '%s' (atom number %d) at position (x, y, z) = (%.2f, %.2f, %.2f) in molecule numbers %s",
                             atom.atomName, serialID, atom.x, atom.y, atom.z, join(moleculeNums,","," and "))
    return newMolUnitCell, transformInfo
end

function translateAtom(molecule::Molecule, serialID::Int64, translation::Vector{Float64})
    newMolecule = deepcopy(molecule)
    for i = 1:length(newMolecule.atomSequence)
        if serialID == newMolecule.atomSequence[i].serialNum
            newMolecule.atomSequence[i].x = newMolecule.atomSequence[i].x + translation[1]
            newMolecule.atomSequence[i].y = newMolecule.atomSequence[i].y + translation[2]
            newMolecule.atomSequence[i].z = newMolecule.atomSequence[i].z + translation[3]
            break
        end
    end
    return newMolecule
end

function transformMoleculeInUnitCell(molUnitCell::MoleculeUnitCell, moleculeNums::Vector{Int64}, serialID::Int64, translation::Vector{Float64})
    newMolUnitCell = deepcopy(molUnitCell)
    atom = newMolUnitCell.moleculeSequence[1].atomSequence[serialID]
    for molNum in moleculeNums
        molecule = newMolUnitCell.moleculeSequence[molNum]
        newMolecule = translateAtom(molecule, serialID, translation)
        newMolUnitCell.moleculeSequence[molNum] = newMolecule
    end
    transformInfo = @sprintf("translation of atom: '%s' (atom number %d) by (x, y, z) = (%.2f, %.2f, %.2f) -
a total of %.2f Angstroms from original position (x, y, z) = (%.2f, %.2f, %.2f) in molecule numbers %s",
                             atom.atomName, serialID, translation[1], translation[2], translation[3], norm(translation), atom.x, atom.y, atom.z,
                             join(moleculeNums,","," and "))
    return newMolUnitCell, transformInfo
end

function translateMolecule(molecule::Molecule, translation::Vector{Float64})
    newMolecule = deepcopy(molecule)
    for atom in newMolecule.atomSequence
        atom.x = atom.x + translation[1]
        atom.y = atom.y + translation[2]
        atom.z = atom.z + translation[3]
    end
    return newMolecule
end

function transformMoleculeInUnitCell(molUnitCell::MoleculeUnitCell, moleculeNums::Vector{Int64}, translation::Vector{Float64})
    newMolUnitCell = deepcopy(molUnitCell)
    for molNum in moleculeNums
        molecule = newMolUnitCell.moleculeSequence[molNum]
        newMolecule = translateMolecule(molecule, translation)
        newMolUnitCell.moleculeSequence[molNum] = newMolecule
    end
    transformInfo = @sprintf("translation of molecules: '%s' by (x, y, z) = (%.2f, %.2f, %.2f) -
a total of %.2f Angstroms from their original positions",
                             join(moleculeNums,","," and "), translation[1], translation[2], translation[3], norm(translation))
    return newMolUnitCell, transformInfo
end

function changeAtom(molecule::Molecule, elementDict::Dict{ASCIIString, Element}, serialID::Int64, newElement::ASCIIString)
    newMolecule = deepcopy(molecule)
    element = elementDict[newElement]
    for i = 1:length(newMolecule.atomSequence)
        if serialID == newMolecule.atomSequence[i].serialNum
            newMolecule.atomSequence[i].element = element
            newMolecule.atomSequence[i].atomName = newElement
            break
        end
    end
    return newMolecule
end

function transformMoleculeInUnitCell(molUnitCell::MoleculeUnitCell, moleculeNums::Vector{Int64}, elementDict::Dict{ASCIIString, Element},
                                     serialID::Int64, newElement::ASCIIString)
    newMolUnitCell = deepcopy(molUnitCell)
    atom = newMolUnitCell.moleculeSequence[1].atomSequence[serialID]
    for molNum in moleculeNums
        molecule = newMolUnitCell.moleculeSequence[molNum]
        newMolecule = changeAtom(molecule, elementDict, serialID, newElement)
        newMolUnitCell.moleculeSequence[molNum] = newMolecule
    end
    transformInfo = @sprintf("change atom from '%s' to '%s' (atom number %d) at position (x, y, z) = (%.2f, %.2f, %.2f) in molecule numbers %s",
                             atom.atomName, newElement, serialID, atom.x, atom.y, atom.z, join(moleculeNums,","," and "))
    return newMolUnitCell, transformInfo
end

function calcMoleculeCentroid(molecule::Molecule)
    x_pos_counter = 0.0
    y_pos_counter = 0.0
    z_pos_counter = 0.0
    counter = 0
    for atom in molecule.atomSequence
        x_pos_counter += atom.x
        y_pos_counter += atom.y
        z_pos_counter += atom.z
        counter += 1
    end
    xCen = x_pos_counter / counter
    yCen = y_pos_counter / counter
    zCen = z_pos_counter / counter

    return [xCen, yCen, zCen]
end

function rotateAtom(atom::Atom, point::Vector{Float64}, dirVec::Vector{Float64}, dirVecLength::Float64, dirVecLengthsqrd::Float64, rotationAngle::Float64,)

    x, y, z = atom.x, atom.y, atom.z
    a, b, c = point[1], point[2], point[3]
    u, v, w = dirVec[1], dirVec[2], dirVec[3]
    L = dirVecLengthsqrd
    sqrtL = dirVecLength
    θ = rotationAngle

    x_new = ((a*(v^2 + w^2) - u*(b*v + c*w - u*x - v*y - w*z))*(1 - cos(θ)) +
        L*x*cos(θ) +
        sqrtL * (-c*v + b*w - w*y + v*z) * sin(θ))/L

    y_new = ((b*(u^2 + w^2) - v*(a*u + c*w - u*x - v*y - w*z))*(1 - cos(θ)) +
        L*y*cos(θ) +
        sqrtL * (c*u - a*w + w*x - u*z) * sin(θ))/L

    z_new = ((c*(v^2 + u^2) - w*(a*u + b*v - u*x - v*y - w*z))*(1 - cos(θ)) +
        L*z*cos(θ) +
        sqrtL * (-b*u + a*v - v*x + u*y) * sin(θ))/L

    return [x_new, y_new, z_new]
end

function rotateMolecule(molecule::Molecule, rotationAngle::Float64, point1::Vector{Float64}=[0.0,0.0,0.0], point2::Vector{Float64}=[0.0,0.0,0.0])
    newMolecule = deepcopy(molecule)
    counter = 1
    p1 = point1
    p2 = point2
    θ  = deg2rad(rotationAngle)
    if p1 == [0.0,0.0,0.0]
        p1 = calcMoleculeCentroid(molecule)
    end
    directionVector = p2 - p1
    u, v, w = directionVector[1], directionVector[2], directionVector[3]
    L = u^2 + v^2 + w^2
    sqrtL = sqrt(L)
    for atom in newMolecule.atomSequence
        new_pos = rotateAtom(atom, point1, directionVector, sqrtL, L, θ)
        newMolecule.atomSequence[counter] = Atom(atom, new_pos)
        counter += 1
    end
    return newMolecule
end

function transformMoleculeInUnitCell(molUnitCell::MoleculeUnitCell, moleculeNums::Vector{Int64}, rotationAngle::Float64,
                                     point1::Vector{Float64}=[0.0,0.0,0.0], point2::Vector{Float64}=[0.0,0.0,0.0])
    newMolUnitCell = deepcopy(molUnitCell)
    for molNum in moleculeNums
        molecule = newMolUnitCell.moleculeSequence[molNum]
        newMolecule = rotateMolecule(molecule, rotationAngle, point1, point2)
        newMolUnitCell.moleculeSequence[molNum] = newMolecule
    end
    if point1 == [0.0,0.0,0.0]
        transformInfo = @sprintf("rotation of molecules: '%s' by %.1f degrees about line passing through the centroid of each molecule and (%s)",
                                 join(moleculeNums,","), rotationAngle, join(point2,","))
    else
        transformInfo = @sprintf("rotation of molecules: '%s' by %.1f degrees about line passing through (%s) and (%s)",
                                 join(moleculeNums,","," and "), rotationAngle, join(point1,","), join(point2,","))
    end
    return newMolUnitCell, transformInfo
end

end
