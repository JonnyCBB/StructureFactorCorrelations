module ElementInformation

using Color

export
    Element,
    ElementDatabase,
    createElementDictionary,
    calcf0!

#############################################################
#####        CREATE TYPES                               #####
#############################################################

type Element
    symbol::ASCIIString
    atomicNumber::Int64
    cm_a1::Float64
    cm_a2::Float64
    cm_a3::Float64
    cm_a4::Float64
    cm_b1::Float64
    cm_b2::Float64
    cm_b3::Float64
    cm_b4::Float64
    cm_c ::Float64
    f0::Dict{Float64,Float64}
    color::RGB
end

#############################################################
#####        FUNCTIONS                                  #####
#############################################################

function createElementDictionary()
    maxColorValue = 255
    elementDictionary = Dict{ASCIIString, Element}()
    atmScatFacDict = Dict{Float64, Float64}()
    cmFile = open("Cromer-Mann_Coefficients.txt")
    for line in eachline(cmFile)
        if line[1:2] == "#S" && length(split(line)) > 3
            symbol = split(line)[3]
            atomicNumber = parseint(split(line)[2])
            a1 = parsefloat(split(line)[4])
            a2 = parsefloat(split(line)[5])
            a3 = parsefloat(split(line)[6])
            a4 = parsefloat(split(line)[7])
            c  = parsefloat(split(line)[8])
            b1 = parsefloat(split(line)[9])
            b2 = parsefloat(split(line)[10])
            b3 = parsefloat(split(line)[11])
            b4 = parsefloat(split(line)[12])
            red   = parseint(split(line)[13])
            green = parseint(split(line)[14])
            blue  = parseint(split(line)[15])
            rgb = RGB(red/maxColorValue, green/maxColorValue, blue/maxColorValue)
            elementInfo = Element(symbol, atomicNumber, a1, a2, a3, a4,
                                  b1, b2, b3, b4, c, atmScatFacDict, rgb)
            elementDictionary[symbol] = elementInfo
        end
    end
    return elementDictionary
end


function calcf0!(elementDict::Dict{ASCIIString,Element}, angles::Vector{Float64}, wavelength::Float64)
    k_sqrd = (sin(deg2rad(angles))/wavelength).^2
    for key in keys(elementDict)
        element = deepcopy(elementDict[key])
        a1 = element.cm_a1
        a2 = element.cm_a2
        a3 = element.cm_a3
        a4 = element.cm_a4
        b1 = element.cm_b1
        b2 = element.cm_b2
        b3 = element.cm_b3
        b4 = element.cm_b4
        c  = element.cm_c
        f0_values = c + a1*exp(-(b1 * k_sqrd)) + a2*exp(-(b2 * k_sqrd)) + a3*exp(-(b3 * k_sqrd)) + a4*exp(-(b4 * k_sqrd))
        for i = 1:length(angles)
           element.f0[angles[i]] = f0_values[i]
        end
        elementDict[key] = element
    end
end

end
