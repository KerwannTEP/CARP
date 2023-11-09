include("../sources/julia/Main.jl")

using HDF5

const JrminMeasure, JrmaxMeasure = _L0*0.02, _L0*0.5 # Jr range
const LminMeasure, LmaxMeasure = _L0*0.02, _L0*1.0 # L range
const nbJrMeasure = 5 # Number of Jr sampling points
const nbLMeasure = 10 # Number of L sampling points

const epsRef = 0.01

const nbCosI = 10

########################################

const tabJrMeasure = collect(range(JrminMeasure,length=nbJrMeasure,JrmaxMeasure))
const tabLMeasure = collect(range(LminMeasure,length=nbLMeasure,LmaxMeasure))
const nbJrLGrid = nbJrMeasure*nbLMeasure

const tabLJrGrid  = zeros(Float64,2,nbJrLGrid) # Location (L,Jr) of the grid points where the diffusion coefficients are computed
const tabdFdt  = zeros(Float64,nbJrLGrid) # Values of dF/dt on the (L,Jr)-grid


########################################



########################################
# Functions to fill the arrays of coefficients
########################################

function tabLJrGrid!()
    index = 1
    for iJr=1:nbJrMeasure
    JrMeasure = tabJrMeasure[iJr]
        for iL=1:nbLMeasure
            LMeasure = tabLMeasure[iL]
            tabLJrGrid[1,index], tabLJrGrid[2,index] = LMeasure, JrMeasure
            index += 1
        end
    end
end

function tabdFdt!()
    

    countbin = Threads.Atomic{Int}(0);
    println("Nb Threads = ",Threads.nthreads())

    Threads.@threads for iGrid=1:nbJrLGrid
        LMeasure, JrMeasure = tabLJrGrid[1,iGrid], tabLJrGrid[2,iGrid]
        dfdt = dFdt2D_JrL(JrMeasure,LMeasure,m_field,alphaRot,nbCosI,nbAvr_default,nbw_default,nbvarphi_default,nbphi_default,nbu0,epsRef)

        tabdFdt[iGrid] = dfdt

        Threads.atomic_add!(countbin, 1)
        println("Progress = ",countbin[],"/",nbJrLGrid)
    end

end


########################################
namefile = "../data/Dump_dFdt_JrL_Map_q_"*string(qCalc)*"_a_"*string(alphaRot)*".hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"q",qCalc)
    write(file,"a",alphaRot)
    write(file,"eps",epsRef)

    write(file,"tabJr",tabJrMeasure)
    write(file,"tabL",tabLMeasure)
    write(file,"tabLJr",tabLJrGrid)
    write(file,"tabdFdt",tabdFdt)

    write(file,"nbJrMeasure",nbJrMeasure)
    write(file,"nbLMeasure",nbLMeasure)


    write(file,"Mtot",_M)
    write(file,"Npart",nbGlobularCluster)

    close(file) # Closing the file
end
########################################

println("Setting up the action-space grid ... ")

@time tabLJrGrid!()

println("Computing the relaxation rate in action space ... ")

@time tabdFdt!()

########################################
writedump!(namefile) # Dumping the computed table
