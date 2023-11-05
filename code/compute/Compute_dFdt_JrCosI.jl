include("../sources/julia/Main.jl")

using HDF5

const JrminMeasure, JrmaxMeasure = _L0*0.02, _L0*0.5 # Jr range
const nbJrMeasure = 25#50 # Number of Jr sampling points

const epsRef = 0.01

const nbL = 50
const Lmax = 3.0


# cosI-sampling for cos I>0
const cosIminPos = 0.02
const cosImaxPos = 0.98
const nbcosIPos = 20#25

const tabcosIPos = [cosIminPos + (cosImaxPos-cosIminPos)*(i-0)/nbcosIPos for i=1:nbcosIPos]
const tabFluxIPos = zeros(Float64,nbcosIPos)


# cosI-sampling for cos I<0
const cosIminNeg = -0.98
const cosImaxNeg = -0.02
const nbcosINeg = 20#25

const tabcosINeg = [cosIminNeg + (cosImaxNeg-cosIminNeg)*(i-1)/nbcosINeg for i=1:nbcosINeg]

const tabCosIMeasure = vcat(tabcosINeg,tabcosIPos)
const nbCosIMeasure = nbcosINeg+nbcosIPos

#println(tabCosIMeasure)

########################################

const tabJrMeasure = collect(range(JrminMeasure,length=nbJrMeasure,JrmaxMeasure))

const nbJrCosIGrid = nbJrMeasure*nbCosIMeasure

const tabCosIJrGrid  = zeros(Float64,2,nbJrCosIGrid) # Location (L,Jr) of the grid points where the diffusion coefficients are computed
const tabdFdt  = zeros(Float64,nbJrCosIGrid) # Values of dF/dt on the (L,Jr)-grid


########################################



########################################
# Functions to fill the arrays of coefficients
########################################

function tabCosIJrGrid!()
    index = 1
    for iJr=1:nbJrMeasure
    JrMeasure = tabJrMeasure[iJr]
        for iCosI=1:nbCosIMeasure
            CosIMeasure = tabCosIMeasure[iCosI]
            tabCosIJrGrid[1,index], tabCosIJrGrid[2,index] = CosIMeasure, JrMeasure
            index += 1
        end
    end
end

function tabdFdt!()

    Threads.@threads for iGrid=1:nbJrCosIGrid
        CosIMeasure, JrMeasure = tabCosIJrGrid[1,iGrid], tabCosIJrGrid[2,iGrid]
        dfdt = dFdtOptiExactSign_2D_JrcosI(JrMeasure,CosIMeasure,m_field,alphaRot,nbL,Lmax,nbAvr_default,nbw_default,nbvarphi_default,nbphi_default,nbu0,epsRef)

        tabdFdt[iGrid] = dfdt
    end

end


########################################
namefile = "../data/Dump_dFdt_JrCosI_Map_q_"*string(qCalc)*"_a_"*string(alphaRot)*".hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"q",qCalc)
    write(file,"a",alphaRot)
    write(file,"eps",epsRef)

    write(file,"tabJr",tabJrMeasure)
    write(file,"tabCosI",tabCosIMeasure)
    write(file,"tabCosIJr",tabCosIJrGrid)
    write(file,"tabdFdt",tabdFdt)

    write(file,"nbJrMeasure",nbJrMeasure)
    write(file,"nbCosIMeasure",nbCosIMeasure)


    write(file,"Mtot",_M)
    write(file,"Npart",nbGlobularCluster)

    close(file) # Closing the file
end
########################################

@time tabCosIJrGrid!()
@time tabdFdt!()

########################################
writedump!(namefile) # Dumping the computed table
