include("Main.jl")

using LinearAlgebra
using HDF5

const nbJr = 50
const nbDeltaL = 50
const nbLzHalf = 50

const Jrmin = 0.01
const Jrmax = 0.5
const DeltaLmin = 0.01
const DeltaLmax = 1.5
const Lzmin = 0.001
const Lzmax = 1.0

const deltaJr = (Jrmax-Jrmin)/nbJr
const deltaL = (DeltaLmax-DeltaLmin)/nbDeltaL
const deltaLz = (Lzmax-Lzmin)/nbLzHalf

const matrixDiffCoeffs_Jr_DeltaL_LzPos = zeros(Float64,9,nbJr,nbDeltaL,nbLzHalf)
const matrixDiffCoeffs_Jr_DeltaL_LzNeg = zeros(Float64,9,nbJr,nbDeltaL,nbLzHalf)

const tabJr = zeros(Float64,nbJr)
const tabDeltaL = zeros(Float64,nbDeltaL)
const tabLzPos = zeros(Float64,nbLzHalf)
const tabLzNeg = zeros(Float64,nbLzHalf)

const nbGridJrL = nbJr*nbDeltaL

println("a")

function filltab!()

    for iJr=1:nbJr
        tabJr[iJr] = Jrmin + deltaJr*(iJr-0.5)
    end

    for iDeltaL=1:nbDeltaL
        tabDeltaL[iDeltaL] = DeltaLmin + deltaL*(iDeltaL-0.5)
    end

    for iLz=1:nbLzHalf
        tabLzPos[iLz] = Lzmin + (Lzmax-Lzmin)/nbLzHalf*(iLz-0.5)
        tabLzNeg[iLz] = -Lzmax + (-Lzmin+Lzmax)/nbLzHalf*(iLz-0.5)
    end

end

@time filltab!()

# println("tabLzPos")
# println(tabLzPos)
# println("tabLzNeg")
# println(tabLzNeg)

function fillMatrixDiff!()

    # Positive Lz
    for iLz=1:nbLzHalf
        Lz = tabLzPos[iLz]

        # println("Lz=",Lz)

        Threads.@threads for iGrid=1:nbGridJrL

            iDeltaL = floor(Int64,iGrid/nbJr)+1
            iJr = iGrid - nbJr*(iDeltaL-1)

            # println("igrid=",iGrid)



            if (iJr == 0)
                iJr = nbJr
                iDeltaL -= 1
            end

            #println("(iJr,iL) = ",(iJr,iDeltaL))

            Jr = tabJr[iJr]
            DeltaL = tabDeltaL[iDeltaL]
            L = abs(Lz) + DeltaL


            dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExact(Jr,L,Lz/L,m_field,alphaRot,nbAvr_default,nbK_default,nbu0,m_field)


            # println((dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz))

            E = E_from_Jr_L(Jr,L,nbu0)


            Frot = _Frot(E,L,Lz,alphaRot)

            matrixDiffCoeffs_Jr_DeltaL_LzPos[1,iJr,iDeltaL,iLz] = dJr*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzPos[2,iJr,iDeltaL,iLz] = dL*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzPos[3,iJr,iDeltaL,iLz] = dLz*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzPos[4,iJr,iDeltaL,iLz] = dJrJr*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzPos[5,iJr,iDeltaL,iLz] = dLL*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr,iDeltaL,iLz] = dLzLz*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzPos[7,iJr,iDeltaL,iLz] = dJrL*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr,iDeltaL,iLz] = dJrLz*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr,iDeltaL,iLz] = dLLz*Frot

        end

    end

    # println("PAUSE")


    # Negative Lz
    for iLz=1:nbLzHalf
        Lz = tabLzNeg[iLz]

        # println("Lz=",Lz)

        Threads.@threads for iGrid=1:nbGridJrL

            # println("igrid=",iGrid)

            iDeltaL = floor(Int64,iGrid/nbJr)+1
            iJr = iGrid - nbJr*(iDeltaL-1)

            if (iJr == 0)
                iJr = nbJr
                iDeltaL -= 1
            end

            Jr = tabJr[iJr]
            DeltaL = tabDeltaL[iDeltaL]
            L = abs(Lz) + DeltaL

            # if (abs(Lz/L)>1.0)
            #     println("halp")
            # end

            dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExact(Jr,L,Lz/L,m_field,alphaRot,nbAvr_default,nbK_default,nbu0,m_field)


            # println((dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz))

            E = E_from_Jr_L(Jr,L,nbu0)


            Frot = _Frot(E,L,Lz,alphaRot)

            matrixDiffCoeffs_Jr_DeltaL_LzNeg[1,iJr,iDeltaL,iLz] = dJr*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzNeg[2,iJr,iDeltaL,iLz] = dL*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzNeg[3,iJr,iDeltaL,iLz] = dLz*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzNeg[4,iJr,iDeltaL,iLz] = dJrJr*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzNeg[5,iJr,iDeltaL,iLz] = dLL*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr,iDeltaL,iLz] = dLzLz*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzNeg[7,iJr,iDeltaL,iLz] = dJrL*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr,iDeltaL,iLz] = dJrLz*Frot
            matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr,iDeltaL,iLz] = dLLz*Frot



        end

    end

end

@time fillMatrixDiff!()

#
# println("matrixDiffCoeffs_Jr_DeltaL_LzPos")
# println(matrixDiffCoeffs_Jr_DeltaL_LzPos)
# println("matrixDiffCoeffs_Jr_DeltaL_LzNeg")
# println(matrixDiffCoeffs_Jr_DeltaL_LzNeg)


# Remove extreme points

const nbJrNew = nbJr-2
const nbDeltaLNew = nbDeltaL-2
const nbLzHalfNew = nbLzHalf-2


const tabJrNew = zeros(Float64,nbJrNew)
const tabDeltaLNew = zeros(Float64,nbDeltaLNew)
const tabLzPosNew = zeros(Float64,nbLzHalfNew)
const tabLzNegNew = zeros(Float64,nbLzHalfNew)

function filltabNew!()

    for iJr=1:nbJrNew
        tabJrNew[iJr] = tabJr[iJr+1]
    end

    for iDeltaL=1:nbDeltaLNew
        tabDeltaLNew[iDeltaL] = tabDeltaL[iDeltaL+1]
    end

    for iLz=1:nbLzHalfNew
        tabLzPosNew[iLz] = tabLzPos[iLz+1]
        tabLzNegNew[iLz] = tabLzNeg[iLz+1]
    end

end

@time filltabNew!()

println("b")

# df(Jr,L,Lz)/dLz
# = d{f(Jr,|Lz|+DeltaL,Lz)}/dLz
# = sign(Lz) * df(Jr,|Lz|+DeltaL,Lz)/dL + df(Jr,|Lz|+DeltaL,Lz)/dLz
# = sign(Lz) * df(Jr,|Lz|+DeltaL,Lz)/dDeltaL + df(Jr,|Lz|+DeltaL,Lz)/dLz

########################################################################################
# Jr-component
########################################################################################


const matrix_GradFluxJr_LzPos = zeros(Float64,nbJrNew,nbDeltaLNew,nbLzHalfNew)
const matrix_GradFluxJr_LzNeg = zeros(Float64,nbJrNew,nbDeltaLNew,nbLzHalfNew)
const matrix_integratedGradFluxJr_LzPos = zeros(Float64,nbJrNew,nbLzHalfNew)
const matrix_integratedGradFluxJr_LzNeg = zeros(Float64,nbJrNew,nbLzHalfNew)


function fillGradientFluxJr!()

    # FJr = dJr F - 0.5*[ d(dJrJr F)/dJr +  d(dJrL F)/dL +  d(dJrLz F)/dLz]
    # d(FJr)/dJr = d(dJr F)/dJr - 0.5*[ d2(dJrJr F)/dJr2 +  d2(dJrL F)/dJrdL +  d2(dJrLz F)/dJrdLz]

    # Positive Lz
    for iLz=1:nbLzHalfNew

        for iJr=1:nbJrNew

            for iDeltaL=1:nbDeltaLNew
                DeltaL = tabDeltaLNew[iDeltaL]

                # CAUTION : ALL INDEX ARE TRANSLATED BY +1

                # d(dJr F)/dJr
                dJrF_p = matrixDiffCoeffs_Jr_DeltaL_LzPos[1,iJr+2,iDeltaL+1,iLz+1]
                dJrF_m = matrixDiffCoeffs_Jr_DeltaL_LzPos[1,iJr,iDeltaL+1,iLz+1]
                dJrF =  (dJrF_p-dJrF_m)/(2.0*deltaJr)

                # println("dJrF_p=",dJrF_p)
                # println("dJrF_m=",dJrF_m)
                # println("dJrF=",dJrF)
                # println("deltaJr=",deltaJr)

                # d2(dJrJr F)/dJrdJr
                dJrJrF_p = matrixDiffCoeffs_Jr_DeltaL_LzPos[4,iJr+2,iDeltaL+1,iLz+1]
                dJrJrF_m = matrixDiffCoeffs_Jr_DeltaL_LzPos[4,iJr,iDeltaL+1,iLz+1]
                dJrJrF = matrixDiffCoeffs_Jr_DeltaL_LzPos[4,iJr+1,iDeltaL+1,iLz+1]

                dJrJrF =  (dJrJrF_p+dJrJrF_m-2.0*dJrJrF)/(deltaJr)^2

                #println("dJrJrF =",dJrJrF )

                # d2(dJrL F)/dJrdL
                dJrLF_pp = matrixDiffCoeffs_Jr_DeltaL_LzPos[7,iJr+2,iDeltaL+2,iLz+1]
                dJrLF_mm = matrixDiffCoeffs_Jr_DeltaL_LzPos[7,iJr,iDeltaL,iLz+1]
                dJrLF_pm = matrixDiffCoeffs_Jr_DeltaL_LzPos[7,iJr+2,iDeltaL,iLz+1]
                dJrLF_mp = matrixDiffCoeffs_Jr_DeltaL_LzPos[7,iJr,iDeltaL+2,iLz+1]

                dJrLF =  (dJrLF_pp+dJrLF_mm-dJrLF_pm-dJrLF_mp)/(4.0*deltaJr*deltaL)

                #println("dJrLF=",dJrLF)
                # d2(dJrLz F)/dJrdLz = sign(Lz) * d2(dJrLz F)/dJrdDeltaL + d2(dJrLz F)/dJrdLz

                # d2(dJrLz F)/dJrdDeltaL

                dJrLzF_pp_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr+2,iDeltaL+2,iLz+1]
                dJrLzF_mm_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr,iDeltaL,iLz+1]
                dJrLzF_pm_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr+2,iDeltaL,iLz+1]
                dJrLzF_mp_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr,iDeltaL+2,iLz+1]

                # d2(dJrLz F)/dJrdLz
                dJrLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr+2,iDeltaL+1,iLz+2]
                dJrLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr,iDeltaL+1,iLz]
                dJrLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr+2,iDeltaL+1,iLz]
                dJrLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr,iDeltaL+1,iLz+2]

                # sign(Lz) = 1

                dJrLzF =  (dJrLzF_pp_1+dJrLzF_mm_1-dJrLzF_pm_1-dJrLzF_mp_1)/(4.0*deltaJr*deltaL) + (dJrLzF_pp_2+dJrLzF_mm_2-dJrLzF_pm_2-dJrLzF_mp_2)/(4.0*deltaJr*deltaLz)

                #println("dJrLzF=",dJrLzF)

                matrix_GradFluxJr_LzPos[iJr,iDeltaL,iLz] = dJrF - 0.5*(dJrJrF+dJrLF+dJrLzF)

                # println("-------")

            end
        end
    end

    # println("matrix_GradFluxJr_LzPos")
    # println(matrix_GradFluxJr_LzPos)
    #
    # println("PAUSE")


    # Negative Lz
    for iLz=1:nbLzHalfNew

        for iJr=1:nbJrNew

            for iDeltaL=1:nbDeltaLNew

                # CAUTION : ALL INDEX ARE TRANSLATED BY +1

                # d(dJr F)/dJr
                dJrF_p = matrixDiffCoeffs_Jr_DeltaL_LzNeg[1,iJr+2,iDeltaL+1,iLz+1]
                dJrF_m = matrixDiffCoeffs_Jr_DeltaL_LzNeg[1,iJr,iDeltaL+1,iLz+1]
                dJrF =  (dJrF_p-dJrF_m)/(2.0*deltaJr)

                # d2(dJrJr F)/dJrdJr
                dJrJrF_p = matrixDiffCoeffs_Jr_DeltaL_LzNeg[4,iJr+2,iDeltaL+1,iLz+1]
                dJrJrF_m = matrixDiffCoeffs_Jr_DeltaL_LzNeg[4,iJr,iDeltaL+1,iLz+1]
                dJrJrF = matrixDiffCoeffs_Jr_DeltaL_LzNeg[4,iJr+1,iDeltaL+1,iLz+1]

                dJrJrF =  (dJrJrF_p+dJrJrF_m-2.0*dJrJrF)/(deltaJr)^2

                # d2(dJrL F)/dJrdL
                dJrLF_pp = matrixDiffCoeffs_Jr_DeltaL_LzNeg[7,iJr+2,iDeltaL+2,iLz+1]
                dJrLF_mm = matrixDiffCoeffs_Jr_DeltaL_LzNeg[7,iJr,iDeltaL,iLz+1]
                dJrLF_pm = matrixDiffCoeffs_Jr_DeltaL_LzNeg[7,iJr+2,iDeltaL,iLz+1]
                dJrLF_mp = matrixDiffCoeffs_Jr_DeltaL_LzNeg[7,iJr,iDeltaL+2,iLz+1]

                dJrLF =  (dJrLF_pp+dJrLF_mm-dJrLF_pm-dJrLF_mp)/(4.0*deltaJr*deltaL)

                # d2(dJrLz F)/dJrdLz = sign(Lz) * d2(dJrLz F)/dJrdDeltaL + d2(dJrLz F)/dJrdLz

                # d2(dJrLz F)/dJrdDeltaL

                dJrLzF_pp_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr+2,iDeltaL+2,iLz+1]
                dJrLzF_mm_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr,iDeltaL,iLz+1]
                dJrLzF_pm_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr+2,iDeltaL,iLz+1]
                dJrLzF_mp_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr,iDeltaL+2,iLz+1]

                # d2(dJrLz F)/dJrdLz
                dJrLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr+2,iDeltaL+1,iLz+2]
                dJrLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr,iDeltaL+1,iLz]
                dJrLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr+2,iDeltaL+1,iLz]
                dJrLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr,iDeltaL+1,iLz+2]

                # sign(Lz) = -1

                dJrLzF =  -(dJrLzF_pp_1+dJrLzF_mm_1-dJrLzF_pm_1-dJrLzF_mp_1)/(4.0*deltaJr*deltaL) + (dJrLzF_pp_2+dJrLzF_mm_2-dJrLzF_pm_2-dJrLzF_mp_2)/(4.0*deltaJr*deltaLz)


                matrix_GradFluxJr_LzNeg[iJr,iDeltaL,iLz] = dJrF - 0.5*(dJrJrF+dJrLF+dJrLzF)

            end
        end
    end

    # println("matrix_GradFluxJr_LzNeg")
    # println(matrix_GradFluxJr_LzNeg)
end

@time fillGradientFluxJr!()



function fillintegratedGradientFluxJr!()

    for iLz=1:nbLzHalfNew
        for iJr=1:nbJrNew
            for iDeltaL=1:nbDeltaLNew

                # Positive Lz
                matrix_integratedGradFluxJr_LzPos[iJr,iLz] += deltaL*matrix_GradFluxJr_LzPos[iJr,iDeltaL,iLz]

                # Negative Lz
                matrix_integratedGradFluxJr_LzNeg[iJr,iLz] += deltaL*matrix_GradFluxJr_LzNeg[iJr,iDeltaL,iLz]
            end
        end
    end
end

@time fillintegratedGradientFluxJr!()

# println("matrix_integratedGradFluxJr_LzPos")
# println(matrix_integratedGradFluxJr_LzPos)
# println("matrix_integratedGradFluxJr_LzNeg")
# println(matrix_integratedGradFluxJr_LzNeg)


########################################################################################
# L-component
########################################################################################


const matrix_GradFluxL_LzPos = zeros(Float64,nbJrNew,nbDeltaLNew,nbLzHalfNew)
const matrix_GradFluxL_LzNeg = zeros(Float64,nbJrNew,nbDeltaLNew,nbLzHalfNew)
const matrix_integratedGradFluxL_LzPos = zeros(Float64,nbJrNew,nbLzHalfNew)
const matrix_integratedGradFluxL_LzNeg = zeros(Float64,nbJrNew,nbLzHalfNew)

println("c")

function fillGradientFluxL!()

    # FL = dL F - 0.5*[ d(dJrL F)/dJr +  d(dLL F)/dL +  d(dLLz F)/dLz]
    # d(FL)/dL = d(dL F)/dJr - 0.5*[ d2(dJrL F)/dJrdL +  d2(dLL F)/dL2 +  d2(dLLz F)/dLdLz]

    # Positive Lz
    for iLz=1:nbLzHalfNew

        for iJr=1:nbJrNew

            for iDeltaL=1:nbDeltaLNew
                DeltaL = tabDeltaLNew[iDeltaL]

                # CAUTION : ALL INDEX ARE TRANSLATED BY +1

                # d(dL F)/dL
                dLF_p = matrixDiffCoeffs_Jr_DeltaL_LzPos[2,iJr+1,iDeltaL+2,iLz+1]
                dLF_m = matrixDiffCoeffs_Jr_DeltaL_LzPos[2,iJr+1,iDeltaL,iLz+1]
                dLF =  (dLF_p-dLF_m)/(2.0*deltaL)

                # d2(dJrL F)/dJrdL
                dJrLF_pp = matrixDiffCoeffs_Jr_DeltaL_LzPos[7,iJr+2,iDeltaL+2,iLz+1]
                dJrLF_mm = matrixDiffCoeffs_Jr_DeltaL_LzPos[7,iJr,iDeltaL,iLz+1]
                dJrLF_pm = matrixDiffCoeffs_Jr_DeltaL_LzPos[7,iJr+2,iDeltaL,iLz+1]
                dJrLF_mp = matrixDiffCoeffs_Jr_DeltaL_LzPos[7,iJr,iDeltaL+2,iLz+1]

                dJrLF =  (dJrLF_pp+dJrLF_mm-dJrLF_pm-dJrLF_mp)/(4.0*deltaJr*deltaL)

                # d2(dLL F)/dL2
                dLLF_p = matrixDiffCoeffs_Jr_DeltaL_LzPos[5,iJr+1,iDeltaL+2,iLz+1]
                dLLF_m = matrixDiffCoeffs_Jr_DeltaL_LzPos[5,iJr+1,iDeltaL,iLz+1]
                dLLF = matrixDiffCoeffs_Jr_DeltaL_LzPos[5,iJr+1,iDeltaL+1,iLz+1]


                dLLF =  (dLLF_p+dLLF_m-2.0*dLLF)/(deltaL)^2

                # d2(dLLz F)/dLdLz = sign(Lz) * d2(dLLz F)/dDeltaL2 + d2(dLLz F)/dLdLz

                # d2(dJrLz F)/dDeltaL2
                dLLzF_p_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL+2,iLz+1]
                dLLzF_m_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL,iLz+1]
                dLLzF_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL+1,iLz+1]

                # d2(dJrLz F)/dLdLz
                dLLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL+2,iLz+2]
                dLLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL,iLz]
                dLLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL+2,iLz]
                dLLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL,iLz+2]

                # sign(Lz) = +1
                dLLzF =  (dLLzF_p_1+dLLzF_m_1-2.0*dLLzF_1)/(deltaL)^2 + (dLLzF_pp_2+dLLzF_mm_2-dLLzF_pm_2-dLLzF_mp_2)/(4.0*deltaJr*deltaLz)

                matrix_GradFluxL_LzPos[iJr,iDeltaL,iLz] = dLF - 0.5*(dJrLF+dLLF+dLLzF)

            end
        end
    end


    # Negative Lz
    for iLz=1:nbLzHalfNew

        for iJr=1:nbJrNew

            for iDeltaL=1:nbDeltaLNew

                # CAUTION : ALL INDEX ARE TRANSLATED BY +1

                # d(dL F)/dL
                dLF_p = matrixDiffCoeffs_Jr_DeltaL_LzNeg[2,iJr+1,iDeltaL+2,iLz+1]
                dLF_m = matrixDiffCoeffs_Jr_DeltaL_LzNeg[2,iJr+1,iDeltaL,iLz+1]
                dLF =  (dLF_p-dLF_m)/(2.0*deltaL)

                # d2(dJrL F)/dJrdL
                dJrLF_pp = matrixDiffCoeffs_Jr_DeltaL_LzNeg[7,iJr+2,iDeltaL+2,iLz+1]
                dJrLF_mm = matrixDiffCoeffs_Jr_DeltaL_LzNeg[7,iJr,iDeltaL,iLz+1]
                dJrLF_pm = matrixDiffCoeffs_Jr_DeltaL_LzNeg[7,iJr+2,iDeltaL,iLz+1]
                dJrLF_mp = matrixDiffCoeffs_Jr_DeltaL_LzNeg[7,iJr,iDeltaL+2,iLz+1]

                dJrLF =  (dJrLF_pp+dJrLF_mm-dJrLF_pm-dJrLF_mp)/(4.0*deltaJr*deltaL)

                # d2(dLL F)/dL2
                dLLF_p = matrixDiffCoeffs_Jr_DeltaL_LzNeg[5,iJr+1,iDeltaL+2,iLz+1]
                dLLF_m = matrixDiffCoeffs_Jr_DeltaL_LzNeg[5,iJr+1,iDeltaL,iLz+1]
                dLLF = matrixDiffCoeffs_Jr_DeltaL_LzNeg[5,iJr+1,iDeltaL+1,iLz+1]


                dLLF =  (dLLF_p+dLLF_m-2.0*dLLF)/(deltaL)^2

                # d2(dLLz F)/dLdLz = sign(Lz) * d2(dLLz F)/dDeltaL2 + d2(dLLz F)/dLdLz

                # d2(dJrLz F)/dDeltaL2
                dLLzF_p_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL+2,iLz+1]
                dLLzF_m_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL,iLz+1]
                dLLzF_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL+1,iLz+1]

                # d2(dJrLz F)/dLdLz
                dLLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL+2,iLz+2]
                dLLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL,iLz]
                dLLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL+2,iLz]
                dLLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL,iLz+2]


                # sign(Lz) = -1
                dLLzF =  -(dLLzF_p_1+dLLzF_m_1-2.0*dLLzF_1)/(deltaL)^2 + (dLLzF_pp_2+dLLzF_mm_2-dLLzF_pm_2-dLLzF_mp_2)/(4.0*deltaJr*deltaLz)

                matrix_GradFluxL_LzNeg[iJr,iDeltaL,iLz] = dLF - 0.5*(dJrLF+dLLF+dLLzF)

            end
        end
    end
end

@time fillGradientFluxL!()

function fillintegratedGradientFluxL!()

    for iLz=1:nbLzHalfNew
        for iJr=1:nbJrNew
            for iDeltaL=1:nbDeltaLNew

                # Positive Lz
                matrix_integratedGradFluxL_LzPos[iJr,iLz] += deltaL*matrix_GradFluxL_LzPos[iJr,iDeltaL,iLz]

                # Negative Lz
                matrix_integratedGradFluxL_LzNeg[iJr,iLz] += deltaL*matrix_GradFluxL_LzNeg[iJr,iDeltaL,iLz]
            end
        end
    end
end

@time fillintegratedGradientFluxL!()

# println("matrix_integratedGradFluxL_LzPos")
# println(matrix_integratedGradFluxL_LzPos)
# println("matrix_integratedGradFluxL_LzNeg")
# println(matrix_integratedGradFluxL_LzNeg)

########################################################################################
# Lz-component
########################################################################################


const matrix_GradFluxLz_LzPos = zeros(Float64,nbJrNew,nbDeltaLNew,nbLzHalfNew)
const matrix_GradFluxLz_LzNeg = zeros(Float64,nbJrNew,nbDeltaLNew,nbLzHalfNew)
const matrix_integratedGradFluxLz_LzPos = zeros(Float64,nbJrNew,nbLzHalfNew)
const matrix_integratedGradFluxLz_LzNeg = zeros(Float64,nbJrNew,nbLzHalfNew)

println("d")


function fillGradientFluxLz!()

    # FLz = dLz F - 0.5*[ d(dJrLz F)/dJr +  d(dLLz F)/dL +  d(dLzLz F)/dLz]
    # d(FLz)/dLz = d(dLz F)/dLz - 0.5*[ d2(dJrLz F)/dJrdLz +  d2(dLLz F)/dLdLz +  d2(dLzLz F)/dLz2]

    # Positive Lz


    for iJr=1:nbJrNew

        for iDeltaL=1:nbDeltaLNew

            # Non LzLz components
            for iLz=1:nbLzHalfNew

                # CAUTION : ALL INDEX ARE TRANSLATED BY +1

                # d(dLz F)/dLz = sign(Lz) * d(dLz F)/dDeltaL + d(dLz F)/dLz

                # d(dLz F)/dDeltaL
                dLzF_p_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[3,iJr+1,iDeltaL+2,iLz+1]
                dLzF_m_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[3,iJr+1,iDeltaL,iLz+1]

                # d(dLz F)/dLz
                dLzF_p_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[3,iJr+1,iDeltaL+1,iLz+2]
                dLzF_m_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[3,iJr+1,iDeltaL+1,iLz]

                # sign(Lz) = 1
                dLzF =  (dLzF_p_1-dLzF_m_1)/(2.0*deltaL) + (dLzF_p_2-dLzF_m_2)/(2.0*deltaLz)




                # d2(dJrLz F)/dJrdLz = sign(Lz) * d2(dJrLz F)/dJrdDeltaL + d2(dJrLz F)/dJrdLz

                # d2(dJrLz F)/dJrdDeltaL

                dJrLzF_pp_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr+2,iDeltaL+2,iLz+1]
                dJrLzF_mm_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr,iDeltaL,iLz+1]
                dJrLzF_pm_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr+2,iDeltaL,iLz+1]
                dJrLzF_mp_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr,iDeltaL+2,iLz+1]

                # d2(dJrLz F)/dJrdLz
                dJrLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr+2,iDeltaL+1,iLz+2]
                dJrLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr,iDeltaL+1,iLz]
                dJrLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr+2,iDeltaL+1,iLz]
                dJrLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[8,iJr,iDeltaL+1,iLz+2]

                # sign(Lz) = 1

                dJrLzF =  (dJrLzF_pp_1+dJrLzF_mm_1-dJrLzF_pm_1-dJrLzF_mp_1)/(4.0*deltaJr*deltaL) + (dJrLzF_pp_2+dJrLzF_mm_2-dJrLzF_pm_2-dJrLzF_mp_2)/(4.0*deltaJr*deltaLz)


                # d2(dLLz F)/dLdLz = sign(Lz) * d2(dLLz F)/dDeltaL2 + d2(dLLz F)/dLdLz

                # d2(dJrLz F)/dDeltaL2
                dLLzF_p_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL+2,iLz+1]
                dLLzF_m_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL,iLz+1]
                dLLzF_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL+1,iLz+1]

                # d2(dJrLz F)/dLdLz
                dLLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL+2,iLz+2]
                dLLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL,iLz]
                dLLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL+2,iLz]
                dLLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[9,iJr+1,iDeltaL,iLz+2]

                # sign(Lz) = +1
                dLLzF =  (dLLzF_p_1+dLLzF_m_1-2.0*dLLzF_1)/(deltaL)^2 + (dLLzF_pp_2+dLLzF_mm_2-dLLzF_pm_2-dLLzF_mp_2)/(4.0*deltaJr*deltaLz)


                # d2(dLzLz F)/dLz2 = d/dLz (sign(Lz) * d(dLzLz F)/dDeltaL + d(dLzLz F)/dLz)
                #  = sign(Lz) * [d/dLz d(dLzLz F)/dDeltaL] + d/dLz d(dLzLz F)/dLz
                #  = sign(Lz) * [sign(Lz) * d2(dLzLz F)/dDeltaL2 + d2(dLzLz F)/dDeltaLdLz] + sign(Lz) * d2(dLzLz F)/dDeltaLdLz + d2(dLzLz F)/dLz2
                #  = d2(dLzLz F)/dDeltaL2  + 2 sign(Lz) * d2(dLzLz F)/dDeltaLdLz + d2(dLzLz F)/dLz2

                # d2(dLzLz F)/dDeltaL2
                dLzLzF_p_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr+1,iDeltaL+2,iLz+1]
                dLzLzF_m_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr+1,iDeltaL,iLz+1]
                dLzLzF_1 = matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr+1,iDeltaL+1,iLz+1]

                # d2(dLzLz F)/dDeltaLdLz
                dLzLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr+1,iDeltaL+2,iLz+2]
                dLzLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr+1,iDeltaL,iLz]
                dLzLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr+1,iDeltaL+2,iLz]
                dLzLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr+1,iDeltaL,iLz+2]

                # sign(Lz) = 1
                dLzLzF = (dLzLzF_p_1+dLzLzF_m_1-2.0*dLzLzF_1)/(deltaL)^2 + 2.0* (dLzLzF_pp_2+dLzLzF_mm_2-dLzLzF_pm_2-dLzLzF_mp_2)/(4.0*deltaL*deltaLz)


                # Fill partially
                matrix_GradFluxLz_LzPos[iJr,iDeltaL,iLz] = dLzF - 0.5*(dJrLzF+dLLzF+dLzLzF)


            end

            # Consider LzLz coefficients d2(dLzLz F)/dLz2
            tabGradDLzLzF = zeros(Float64,nbLzHalfNew)
            for iLz=1:nbLzHalfNew
                DLzLzF_p = matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr+1,iDeltaL+1,iLz+2]
                DLzLzF_m = matrixDiffCoeffs_Jr_DeltaL_LzPos[6,iJr+1,iDeltaL+1,iLz]
                tabGradDLzLzF[iLz] = (DLzLzF_p-DLzLzF_m)/(2.0*deltaLz)

            end

            # Interpolation least-square fit : degree-4 polynomial

            matX = [tabLzPosNew[iLz]^j for iLz=1:nbLzHalfNew,j=0:5]
            tmatX = transpose(matX)
            S = inv(tmatX*matX)
            beta = S*tmatX*tabGradDLzLzF

            for iLz=1:nbLzHalfNew
                Lz = tabLzPosNew[iLz]
                matrix_GradFluxLz_LzPos[iJr,iDeltaL,iLz] += - 0.5*(beta[2]+2.0*beta[3]*Lz+3.0*beta[4]*Lz^2+4.0*beta[5]*Lz^3+5.0*beta[6]*Lz^4)
            end



        end
    end





    # Negative Lz


    for iJr=1:nbJrNew

        for iDeltaL=1:nbDeltaLNew

            # Non LzLz components
            for iLz=1:nbLzHalfNew

                # CAUTION : ALL INDEX ARE TRANSLATED BY +1

                # d(dLz F)/dLz = sign(Lz) * d(dLz F)/dDeltaL + d(dLz F)/dLz

                # d(dLz F)/dDeltaL
                dLzF_p_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[3,iJr+1,iDeltaL+2,iLz+1]
                dLzF_m_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[3,iJr+1,iDeltaL,iLz+1]

                # d(dLz F)/dLz
                dLzF_p_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[3,iJr+1,iDeltaL+1,iLz+2]
                dLzF_m_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[3,iJr+1,iDeltaL+1,iLz]

                # sign(Lz) = -1
                dLzF =  -(dLzF_p_1-dLzF_m_1)/(2.0*deltaL) + (dLzF_p_2-dLzF_m_2)/(2.0*deltaLz)




                # d2(dJrLz F)/dJrdLz = sign(Lz) * d2(dJrLz F)/dJrdDeltaL + d2(dJrLz F)/dJrdLz

                # d2(dJrLz F)/dJrdDeltaL

                dJrLzF_pp_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr+2,iDeltaL+2,iLz+1]
                dJrLzF_mm_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr,iDeltaL,iLz+1]
                dJrLzF_pm_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr+2,iDeltaL,iLz+1]
                dJrLzF_mp_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr,iDeltaL+2,iLz+1]

                # d2(dJrLz F)/dJrdLz
                dJrLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr+2,iDeltaL+1,iLz+2]
                dJrLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr,iDeltaL+1,iLz]
                dJrLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr+2,iDeltaL+1,iLz]
                dJrLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[8,iJr,iDeltaL+1,iLz+2]

                # sign(Lz) = -1

                dJrLzF =  -(dJrLzF_pp_1+dJrLzF_mm_1-dJrLzF_pm_1-dJrLzF_mp_1)/(4.0*deltaJr*deltaL) + (dJrLzF_pp_2+dJrLzF_mm_2-dJrLzF_pm_2-dJrLzF_mp_2)/(4.0*deltaJr*deltaLz)


                # d2(dLLz F)/dLdLz = sign(Lz) * d2(dLLz F)/dDeltaL2 + d2(dLLz F)/dLdLz

                # d2(dJrLz F)/dDeltaL2
                dLLzF_p_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL+2,iLz+1]
                dLLzF_m_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL,iLz+1]
                dLLzF_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL+1,iLz+1]

                # d2(dJrLz F)/dLdLz
                dLLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL+2,iLz+2]
                dLLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL,iLz]
                dLLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL+2,iLz]
                dLLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[9,iJr+1,iDeltaL,iLz+2]

                # sign(Lz) = -1
                dLLzF =  -(dLLzF_p_1+dLLzF_m_1-2.0*dLLzF_1)/(deltaL)^2 + (dLLzF_pp_2+dLLzF_mm_2-dLLzF_pm_2-dLLzF_mp_2)/(4.0*deltaJr*deltaLz)


                # d2(dLzLz F)/dLz2 = d/dLz (sign(Lz) * d(dLzLz F)/dDeltaL + d(dLzLz F)/dLz)
                #  = sign(Lz) * [d/dLz d(dLzLz F)/dDeltaL] + d/dLz d(dLzLz F)/dLz
                #  = sign(Lz) * [sign(Lz) * d2(dLzLz F)/dDeltaL2 + d2(dLzLz F)/dDeltaLdLz] + sign(Lz) * d2(dLzLz F)/dDeltaLdLz + d2(dLzLz F)/dLz2
                #  = d2(dLzLz F)/dDeltaL2  + 2 sign(Lz) * d2(dLzLz F)/dDeltaLdLz + d2(dLzLz F)/dLz2

                # d2(dLzLz F)/dDeltaL2
                dLzLzF_p_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr+1,iDeltaL+2,iLz+1]
                dLzLzF_m_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr+1,iDeltaL,iLz+1]
                dLzLzF_1 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr+1,iDeltaL+1,iLz+1]

                # d2(dLzLz F)/dDeltaLdLz
                dLzLzF_pp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr+1,iDeltaL+2,iLz+2]
                dLzLzF_mm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr+1,iDeltaL,iLz]
                dLzLzF_pm_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr+1,iDeltaL+2,iLz]
                dLzLzF_mp_2 = matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr+1,iDeltaL,iLz+2]

                # sign(Lz) = -1
                dLzLzF = (dLzLzF_p_1+dLzLzF_m_1-2.0*dLzLzF_1)/(deltaL)^2 - 2.0* (dLzLzF_pp_2+dLzLzF_mm_2-dLzLzF_pm_2-dLzLzF_mp_2)/(4.0*deltaL*deltaLz)


                # Fill partially
                matrix_GradFluxLz_LzPos[iJr,iDeltaL,iLz] = dLzF - 0.5*(dJrLzF+dLLzF+dLzLzF)


            end

            # Consider LzLz coefficients d2(dLzLz F)/dLz2
            tabGradDLzLzF = zeros(Float64,nbLzHalfNew)
            for iLz=1:nbLzHalfNew
                DLzLzF_p = matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr+1,iDeltaL+1,iLz+2]
                DLzLzF_m = matrixDiffCoeffs_Jr_DeltaL_LzNeg[6,iJr+1,iDeltaL+1,iLz]
                tabGradDLzLzF[iLz] = (DLzLzF_p-DLzLzF_m)/(2.0*deltaLz)

            end

            # Interpolation least-square fit : degree-4 polynomial

            matX = [tabLzPosNew[iLz]^j for iLz=1:nbLzHalfNew,j=0:5]
            tmatX = transpose(matX)
            S = inv(tmatX*matX)
            beta = S*tmatX*tabGradDLzLzF

            for iLz=1:nbLzHalfNew
                Lz = tabLzNegNew[iLz]
                matrix_GradFluxLz_LzNeg[iJr,iDeltaL,iLz] += - 0.5*(beta[2]+2.0*beta[3]*Lz+3.0*beta[4]*Lz^2+4.0*beta[5]*Lz^3+5.0*beta[6]*Lz^4)
            end



        end
    end

end


@time fillGradientFluxLz!()


# println("matrix_GradFluxLz_LzPos")
# println(matrix_GradFluxLz_LzPos)
# println("matrix_GradFluxLz_LzNeg")
# println(matrix_GradFluxLz_LzNeg)


function fillintegratedGradientFluxLz!()

    for iLz=1:nbLzHalfNew
        for iJr=1:nbJrNew
            for iDeltaL=1:nbDeltaLNew

                # Positive Lz
                matrix_integratedGradFluxLz_LzPos[iJr,iLz] += deltaL*matrix_GradFluxLz_LzPos[iJr,iDeltaL,iLz]

                # Negative Lz
                matrix_integratedGradFluxLz_LzNeg[iJr,iLz] += deltaL*matrix_GradFluxLz_LzNeg[iJr,iDeltaL,iLz]
            end
        end
    end
end

@time fillintegratedGradientFluxLz!()


########################################################################################
# Put together
########################################################################################

const matrix_dFdt_JrLz = zeros(Float64,nbJrNew,2*nbLzHalfNew)


function filldFdt_JrLz!()

    for iJr=1:nbJrNew
        for iLz=1:nbLzHalfNew

            # Lz positive
            matrix_dFdt_JrLz[iJr,nbLzHalfNew+iLz] = -(matrix_integratedGradFluxJr_LzPos[iJr,iLz]+matrix_integratedGradFluxL_LzPos[iJr,iLz]+matrix_integratedGradFluxLz_LzPos[iJr,iLz])

            # Lz negative
            matrix_dFdt_JrLz[iJr,iLz] = -(matrix_integratedGradFluxJr_LzNeg[iJr,iLz]+matrix_integratedGradFluxL_LzNeg[iJr,iLz]+matrix_integratedGradFluxLz_LzNeg[iJr,iLz])
        end
    end

end


println("e")


@time filldFdt_JrLz!()

println("tabdFdt")
println(matrix_dFdt_JrLz)



const tabLzJrGrid = zeros(Float64,2,nbJrNew*2*nbLzHalfNew)
const tabdFdt = zeros(Float64,nbJrNew*2*nbLzHalfNew)

function tabLzJr!()
    index = 1
    for iJr=1:nbJrNew
        JrMeasure = tabJrNew[iJr]
        for iLz=1:nbLzHalfNew
            # Neg
            LzMeasure = tabLzNegNew[iLz]
            tabLzJrGrid[1,index], tabLzJrGrid[2,index] = LzMeasure, JrMeasure
            index += 1
        end

        for iLz=1:nbLzHalfNew
            # Pos
            LzMeasure = tabLzPosNew[iLz]
            tabLzJrGrid[1,index], tabLzJrGrid[2,index] = LzMeasure, JrMeasure
            index += 1
        end
    end
end

@time tabLzJr!()


function tabdFdt!()

    index = 1
    for iJr=1:nbJrNew

        for iLz=1:2*nbLzHalfNew

            tabdFdt[index] = matrix_dFdt_JrLz[iJr,iLz]
            # if (isnan(tabdFdt[index]))
            #     println("(Lz,Jr) = ",(tabLzJrGrid[1,index], tabLzJrGrid[2,index]))
            # end
            index += 1


        end

    end

end

@time tabdFdt!()

#println("f")

########################################################################################
# Dump in file
########################################################################################


########################################
namefile = "../../data/Dump_dFdt_JrLz_InterpolateLzLz_q_"*string(qCalc)*"_alpha_"*string(alphaRot)*".hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"q",qCalc)
    write(file,"alpha",alphaRot)

    write(file,"tabLzJr",tabLzJrGrid)
    write(file,"tabdFdt",tabdFdt)




    close(file) # Closing the file
end
########################################


########################################
@time writedump!(namefile) # Dumping the computed table
