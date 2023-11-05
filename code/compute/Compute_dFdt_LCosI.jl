include("../sources/julia/Main.jl")

using HDF5

const LminMeasure, LmaxMeasure = _L0*0.02, _L0*1.0 # Jr range
const nbLMeasure = 25 # Number of Jr sampling points

const epsRef = 0.01

const nbJr = 50
const Jrmax = 3.0


# cosI-sampling for cos I>0
const cosIminPos = 0.02
const cosImaxPos = 0.98
const nbcosIPos = 20

const tabcosIPos = [cosIminPos + (cosImaxPos-cosIminPos)*(i-0)/nbcosIPos for i=1:nbcosIPos]
const tabFluxIPos = zeros(Float64,nbcosIPos)


# cosI-sampling for cos I<0
const cosIminNeg = -0.98
const cosImaxNeg = -0.02
const nbcosINeg = 20

const tabcosINeg = [cosIminNeg + (cosImaxNeg-cosIminNeg)*(i-1)/nbcosINeg for i=1:nbcosINeg]

const tabCosIMeasure = vcat(tabcosINeg,tabcosIPos)
const nbCosIMeasure = nbcosINeg+nbcosIPos


########################################

const tabLMeasure = collect(range(LminMeasure,length=nbLMeasure,LmaxMeasure))

const nbLCosIGrid = nbLMeasure*nbCosIMeasure

const tabCosILGrid  = zeros(Float64,2,nbLCosIGrid) # Location (L,Jr) of the grid points where the diffusion coefficients are computed
const tabdFdt  = zeros(Float64,nbLCosIGrid) # Values of dF/dt on the (L,Jr)-grid


########################################



########################################
# Functions to fill the arrays of coefficients
########################################

function tabCosILGrid!()
    index = 1
    for iL=1:nbLMeasure
    LMeasure = tabLMeasure[iL]
        for iCosI=1:nbCosIMeasure
            CosIMeasure = tabCosIMeasure[iCosI]
            tabCosILGrid[1,index], tabCosILGrid[2,index] = CosIMeasure, LMeasure
            index += 1
        end
    end
end

function tabdFdt!()

    println("Nb Threads = ",Threads.nthreads())

    Threads.@threads for iGrid=1:nbLCosIGrid
        CosIMeasure, LMeasure = tabCosILGrid[1,iGrid], tabCosILGrid[2,iGrid]
        dfdt = dFdt2D_LcosI(LMeasure,CosIMeasure,m_field,alphaRot,nbJr,Jrmax,nbAvr_default,nbw_default,nbvarphi_default,nbphi_default,nbu0,epsRef)

        tabdFdt[iGrid] = dfdt
    end

end


########################################
namefile = "../data/Dump_dFdt_LCosI_Map_q_"*string(qCalc)*"_a_"*string(alphaRot)*".hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"q",qCalc)
    write(file,"a",alphaRot)
    write(file,"eps",epsRef)

    write(file,"tabL",tabLMeasure)
    write(file,"tabCosI",tabCosIMeasure)
    write(file,"tabCosIL",tabCosILGrid)
    write(file,"tabdFdt",tabdFdt)

    write(file,"nbLMeasure",nbLMeasure)
    write(file,"nbCosIMeasure",nbCosIMeasure)


    write(file,"Mtot",_M)
    write(file,"Npart",nbGlobularCluster)

    close(file) # Closing the file
end
########################################

@time tabCosILGrid!()
@time tabdFdt!()

########################################
writedump!(namefile) # Dumping the computed table
