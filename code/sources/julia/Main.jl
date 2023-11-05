include("Args.jl") # Parsing the command-line
##################################################
# General parameters
##################################################
"""
    PARALLEL

Determining if the code is run in parallel.
"""
const PARALLEL = parsed_args["parallel"]
if ((PARALLEL != "yes") && (PARALLEL != "no"))
    error("ERROR: UNKNOWN PARALLEL") # Unknown parallel procedure
end


"""
    nbK_default

Default number of sampling points for the 3D integrals.
"""
const nbK_default = parsed_args["nbK"]

"""
    nbAvr_default

Default number of sampling points for the orbit-averaging integrals.
"""
const nbAvr_default = parsed_args["nbAvr"]

"""
    qCalc

Anisotropy parameter q for the Plummer model.
"""
const qCalc = parsed_args["q"]

"""
    alpha

Rotation parameter for the Lynden-Bell Demon.
Frot(E,L,Lz) = F(E,L) * (1 + alpha * g(Lz/L))
"""
const alphaRot = parsed_args["a"]

"""
    g(x)

Rotation fonction g(cos I)=g(Lz/L)
Frot(E,L,Lz) = F(E,L) * (1 + alpha * g(Lz/L))
"""
const gfunction = parsed_args["g"]

const nbAvrTh_default = parsed_args["nbAvrTh"]
const nbu0 = parsed_args["nbu"]
const a_err = parsed_args["a_err"]

const nbw_default = parsed_args["nbw"]
const nbphi_default = parsed_args["nbphi"]
const nbvarphi_default = parsed_args["nbvarphi"]

##################################################
##################################################

include("Constants.jl")
include("Mean.jl")
include("EffectiveAnomaly.jl")

include("Grad_Jr.jl")
include("AlphaBeta.jl")
include("Inversion.jl")
include("OrbitParameters.jl")
include("Bath.jl")
include("DistributionFunctionLz.jl")
include("LocalDeflectionAngleAverage.jl")
include("OrbitAverage.jl")
include("Flux.jl")
