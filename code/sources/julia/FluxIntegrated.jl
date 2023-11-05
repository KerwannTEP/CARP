include("Main.jl")

# function dFdtIntegratedLz(Jr::Float64, L::Float64, m_field::Float64,
#             nbAvr::Int64=nbAvr_default, nbK::Int64=nbK_default, nbu::Int64=nbu0,
#             eps::Float64=10^(-2), m_test::Float64=m_field)
#
#     # Case Jr,L > eps
#


function LeftoverNoRot(Jr::Float64, L::Float64, m_field::Float64,
            nbAvr::Int64=nbAvr_default, nbK::Int64=nbK_default, nbu::Int64=nbu0,
            eps::Float64=10^(-2), m_test::Float64=m_field)

    # Non-derivative


    _, _, dI_1, _, _, dII_1, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    _, _, dI_m1, _, _, dII_m1, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,-1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)

    E = E_from_Jr_L(Jr,L)
    Ftot = L*_F(E,L)

    term0 = dI_1*Ftot - dI_m1*Ftot

    # (Jr,L) derivatives

    Jr_p = Jr + _L0*eps
    Jr_m = Jr - _L0*eps
    L_p = L + _L0*eps
    L_m = L - _L0*eps
    E_rp = E_from_Jr_L(Jr_p,L)
    E_rm = E_from_Jr_L(Jr_m,L)
    E_Lp = E_from_Jr_L(Jr,L_p)
    E_Lm = E_from_Jr_L(Jr,L_m)

    _, _, _, _, _, _, _, dJrI_rp1, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_p,L,1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    Ftot_rp1 = L*_F(E_rp,L)
    _, _, _, _, _, _, _, dJrI_rm1, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_m,L,1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    Ftot_rm1 = L*_F(E_rm,L)
    dJr1 = (dJrI_rp1*Ftot_rp1-dJrI_rm1*Ftot_rm1)/(2.0*_L0*eps)

    _, _, _, _, _, _, _, dJrI_rpm1, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_p,L,-1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    Ftot_rpm1 = L*_F(E_rp,L)
    _, _, _, _, _, _, _, dJrI_rmm1, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_m,L,-1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    Ftot_rmm1 = L*_F(E_rm,L)
    dJrm1 = (dJrI_rpm1*Ftot_rpm1-dJrI_rmm1*Ftot_rmm1)/(2.0*_L0*eps)




    _, _, _, _, _, _, _, _, dLI_Lp1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    Ftot_Lp1 = L_p*_F(E_Lp,L_p)
    _, _, _, _, _, _, _, _, dLI_Lm1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_m,1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    Ftot_Lm1 = L_m*_F(E_Lm,L_m)

    dL1 = (dLI_Lp1*Ftot_Lp1-dLI_Lm1*Ftot_Lm1)/(2.0*_L0*eps)

    _, _, _, _, _, _, _, _, dLI_Lpm1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,-1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    Ftot_Lpm1 = L_p*_F(E_Lp,L_p)
    _, _, _, _, _, _, _, _, dLI_Lmm1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_m,-1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    Ftot_Lmm1 = L_m*_F(E_Lm,L_m)

    dLm1 = (dLI_Lpm1*Ftot_Lpm1-dLI_Lmm1*Ftot_Lmm1)/(2.0*_L0*eps)



    term1 = dJr1+dL1-dJrm1-dLm1

    # cosI derivative
    # use forward finite diff
    # 3-point stencil
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient

    cosI_1_m = 1.0 - eps
    cosI_1_mm = 1.0 - 2.0*eps

    cosI_m1_p = -1.0 + eps
    cosI_m1_pp = -1.0 + 2.0*eps


    _, _, _, _, _, dII_1_m, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_1_m,m_field,0.0,nbAvr,nbK,nbu,m_test)
    _, _, _, _, _, dII_1_mm, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_1_mm,m_field,0.0,nbAvr,nbK,nbu,m_test)
    _, _, _, _, _, dII_m1_p, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_m1_p,m_field,0.0,nbAvr,nbK,nbu,m_test)
    _, _, _, _, _, dII_m1_pp, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_m1_pp,m_field,0.0,nbAvr,nbK,nbu,m_test)

    # f'(1) = [1.5*f(1)-2.0*f(1-eps)+0.5*f(1-2*eps)/eps
    dcosI_1 = (1.5*dII_1 - 2.0*dII_1_m + 0.5*dII_1_mm)*Ftot/(eps)

    # f'(-1) = [-1.5*f(-1)+2.0*f(-1+eps)-0.5*f(-1+2*eps)/eps
    dcosI_m1 = (-1.5*dII_m1 + 2.0*dII_m1_p - 0.5*dII_m1_pp)*Ftot/(eps)

    println("d/dcosI at cosI = 1 :")
    println("naive  : ",(dII_1-dII_1_m)/eps*Ftot)
    println("precise: ",dcosI_1)
    println("------------------------")
    println("d/dcosI at cosI = -1 :")
    println("naive  : ",(dII_m1_p-dII_m1)/eps*Ftot)
    println("precise: ",dcosI_m1)
    println("------------------------")

    term2 = 0.5*(dcosI_1-dcosI_m1)

    println("(0,1,2) = ",(term0,term1,term2))

    return term0 - term1 - term2
end



function Test_FD_NoRot(Jr::Float64, L::Float64, m_field::Float64,
            nbAvr::Int64=nbAvr_default, nbK::Int64=nbK_default, nbu::Int64=nbu0,
            eps::Float64=10^(-2), m_test::Float64=m_field)

    # Non-derivative


    _, _, dI_1, _, _, dII_1, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)
    _, _, dI_m1, _, _, dII_m1, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,-1.0,m_field,0.0,nbAvr,nbK,nbu,m_test)

    E = E_from_Jr_L(Jr,L)
    Ftot = L*_F(E,L)

    term0 = dI_1*Ftot - dI_m1*Ftot



    # cosI derivative
    # use forward finite diff
    # 3-point stencil
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient

    cosI_1_m = 1.0 - eps
    cosI_1_mm = 1.0 - 2.0*eps

    cosI_m1_p = -1.0 + eps
    cosI_m1_pp = -1.0 + 2.0*eps


    _, _, _, _, _, dII_1_m, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_1_m,m_field,0.0,nbAvr,nbK,nbu,m_test)
    _, _, _, _, _, dII_1_mm, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_1_mm,m_field,0.0,nbAvr,nbK,nbu,m_test)
    _, _, _, _, _, dII_m1_p, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_m1_p,m_field,0.0,nbAvr,nbK,nbu,m_test)
    _, _, _, _, _, dII_m1_pp, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_m1_pp,m_field,0.0,nbAvr,nbK,nbu,m_test)

    # f'(1) = [1.5*f(1)-2.0*f(1-eps)+0.5*f(1-2*eps)/eps
    dcosI_1 = (1.5*dII_1 - 2.0*dII_1_m + 0.5*dII_1_mm)*Ftot/(eps)

    # f'(-1) = [-1.5*f(-1)+2.0*f(-1+eps)-0.5*f(-1+2*eps)/eps
    dcosI_m1 = (-1.5*dII_m1 + 2.0*dII_m1_p - 0.5*dII_m1_pp)*Ftot/(eps)


    term2 = 0.5*(dcosI_1-dcosI_m1)

    println("(DI,grad DII) = ",(term0,term2))
    println("FD relation   = ",term0-term2)

    return term0, term2


end





function Leftover(Jr::Float64, L::Float64, m_field::Float64, alpha::Float64=alphaRot,
            nbAvr::Int64=nbAvr_default, nbK::Int64=nbK_default, nbu::Int64=nbu0,
            eps::Float64=10^(-2), m_test::Float64=m_field)

    # Non-derivative


    _, _, dI_1, _, _, dII_1, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
    _, _, dI_m1, _, _, dII_m1, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)

    E = E_from_Jr_L(Jr,L)
    Ftot_1 = _Frot_cosI(E,L,1.0)
    Ftot_m1 = _Frot_cosI(E,L,-1.0)

    term0 = dI_1*Ftot_1 - dI_m1*Ftot_m1

    # (Jr,L) derivatives

    Jr_p = Jr + _L0*eps
    Jr_m = Jr - _L0*eps
    L_p = L + _L0*eps
    L_m = L - _L0*eps
    E_rp = E_from_Jr_L(Jr_p,L)
    E_rm = E_from_Jr_L(Jr_m,L)
    E_Lp = E_from_Jr_L(Jr,L_p)
    E_Lm = E_from_Jr_L(Jr,L_m)

    _, _, _, _, _, _, _, dJrI_rp1, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_p,L,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_rp1 = _Frot_cosI(E_rp,L,1.0)
    _, _, _, _, _, _, _, dJrI_rm1, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_m,L,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_rm1 = _Frot_cosI(E_rm,L,1.0)
    dJr1 = (dJrI_rp1*Ftot_rp1-dJrI_rm1*Ftot_rm1)/(2.0*_L0*eps)

    _, _, _, _, _, _, _, dJrI_rpm1, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_p,L,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_rpm1 = _Frot_cosI(E_rp,L,-1.0)
    _, _, _, _, _, _, _, dJrI_rmm1, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_m,L,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_rmm1 = _Frot_cosI(E_rm,L,-1.0)
    dJrm1 = (dJrI_rpm1*Ftot_rpm1-dJrI_rmm1*Ftot_rmm1)/(2.0*_L0*eps)


    if (L > eps)



        _, _, _, _, _, _, _, _, dLI_Lp1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lp1 = _Frot_cosI(E_Lp,L_p,1.0)
        _, _, _, _, _, _, _, _, dLI_Lm1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_m,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lm1 = _Frot_cosI(E_Lm,L_m,1.0)

        dL1 = (dLI_Lp1*Ftot_Lp1-dLI_Lm1*Ftot_Lm1)/(2.0*_L0*eps)

        _, _, _, _, _, _, _, _, dLI_Lpm1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lpm1 = _Frot_cosI(E_Lp,L_p,-1.0)
        _, _, _, _, _, _, _, _, dLI_Lmm1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_m,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lmm1 = _Frot_cosI(E_Lm,L_m,-1.0)

        dLm1 = (dLI_Lpm1*Ftot_Lpm1-dLI_Lmm1*Ftot_Lmm1)/(2.0*_L0*eps)

    else # use forward diff

        L_pp = L + 2.0*_L0*eps
        E_Lpp = E_from_Jr_L(Jr,L_pp)

        _, _, _, _, _, _, _, _, dLI_L_m1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_L_m1 = _Frot_cosI(E,L,-1.0)
        _, _, _, _, _, _, _, _, dLI_Lp_m1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lp_m1 = _Frot_cosI(E_Lp,L_p,-1.0)
        _, _, _, _, _, _, _, _, dLI_Lpp_m1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lpp_m1 = _Frot_cosI(E_Lpp,L_pp,-1.0)

        dLm1 = (-1.5*dLI_L_m1*Ftot_L_m1 + 2.0*dLI_Lp_m1*Ftot_Lp_m1 - 0.5*dLI_Lpp_m1*Ftot_Lpp_m1)/(eps)

        _, _, _, _, _, _, _, _, dLI_L_1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_L_1 = _Frot_cosI(E,L,1.0)
        _, _, _, _, _, _, _, _, dLI_Lp_1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lp_1 = _Frot_cosI(E_Lp,L_p,1.0)
        _, _, _, _, _, _, _, _, dLI_Lpp_1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lpp_1 = _Frot_cosI(E_Lpp,L_pp,1.0)

        dL1 = (-1.5*dLI_L_m1*Ftot_L_m1 + 2.0*dLI_Lp_m1*Ftot_Lp_m1 - 0.5*dLI_Lpp_m1*Ftot_Lpp_m1)/(eps)


    end



    term1 = dJr1+dL1-dJrm1-dLm1

    # cosI derivative
    # use forward finite diff
    # 3-point stencil
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient

    cosI_1_m = 1.0 - eps
    cosI_1_mm = 1.0 - 2.0*eps

    Ftot_1_m = _Frot_cosI(E,L,cosI_1_m)
    Ftot_1_mm = _Frot_cosI(E,L,cosI_1_mm)

    cosI_m1_p = -1.0 + eps
    cosI_m1_pp = -1.0 + 2.0*eps

    Ftot_m1_p = _Frot_cosI(E,L,cosI_m1_p)
    Ftot_m1_pp = _Frot_cosI(E,L,cosI_m1_pp)


    _, _, _, _, _, dII_1_m, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_1_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    _, _, _, _, _, dII_1_mm, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_1_mm,m_field,alpha,nbAvr,nbK,nbu,m_test)
    _, _, _, _, _, dII_m1_p, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_m1_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    _, _, _, _, _, dII_m1_pp, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_m1_pp,m_field,alpha,nbAvr,nbK,nbu,m_test)

    # f'(1) = [1.5*f(1)-2.0*f(1-eps)+0.5*f(1-2*eps)/eps
    dcosI_1 = (1.5*dII_1*Ftot_1 - 2.0*dII_1_m*Ftot_1_m + 0.5*dII_1_mm*Ftot_1_mm)/(eps)

    # f'(-1) = [-1.5*f(-1)+2.0*f(-1+eps)-0.5*f(-1+2*eps)/eps
    dcosI_m1 = (-1.5*dII_m1*Ftot_m1 + 2.0*dII_m1_p*Ftot_m1_p - 0.5*dII_m1_pp*Ftot_m1_pp)/(eps)

    # println("d/dcosI at cosI = 1 :")
    # println("naive  : ",(dII_1-dII_1_m)/eps*Ftot)
    # println("precise: ",dcosI_1)
    # println("------------------------")
    # println("d/dcosI at cosI = -1 :")
    # println("naive  : ",(dII_m1_p-dII_m1)/eps*Ftot)
    # println("precise: ",dcosI_m1)
    # println("------------------------")

    term2 = 0.5*(dcosI_1-dcosI_m1)

    #println("(0,1,2) = ",(term0,term1,term2))

    return term0 - term1 - term2
end



###########################################################
# (Jr,L)-Integrated diffusion coefficients
###########################################################

function dJrF_IntegratedPar(Lz::Float64,m_field::Float64, DeltaLmax::Float64=3.0*_L0, nbL::Int64=30,
                        alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                        nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    sumJr = Threads.Atomic{Float64}(0.0)
    sumJrLz_1 = Threads.Atomic{Float64}(0.0)


    Threads.@threads for iL=1:nbL
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        E = E_from_Jr_L(0.0,L)
        Frot = _Frot(E,L,Lz,alpha)
        avrDJr, _, _, _, _, _, _, avrDJrLz, _ = orbitAverageActionCoeffsOptiExact(0.0,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)



        Threads.atomic_add!(sumJr,avrDJr*Frot)
        Threads.atomic_add!(sumJrLz_1,avrDJrLz*Frot)


    end

    sumJr[] *= -DeltaLmax/nbL
    sumJrLz_1[] *= -DeltaLmax/nbL

    _, _, _, _, _, _, _, avrDJrLz, _ = orbitAverageActionCoeffsOptiExact(0.0,abs(Lz),sign(Lz),m_field,alpha,nbAvr,nbK,nbu,m_test)

    sumJrLz_2 = -sign(Lz)*avrDJrLz

    return sumJr[], sumJrLz_1[], sumJrLz_2
end


function dLF_IntegratedPar(Lz::Float64,m_field::Float64, Jrmax::Float64=3.0*_L0, nbJr::Int64=30,
                        alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, eps::Float64=10^(-2),
                        nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    sumL = Threads.Atomic{Float64}(0.0)
    sumLLz_1 = Threads.Atomic{Float64}(0.0)
    sumLLz_2 = Threads.Atomic{Float64}(0.0)

    Threads.@threads for iJr=1:nbJr
        Jr = Jrmax*(iJr-0.5)/nbJr
        E = E_from_Jr_L(Jr,abs(Lz))
        E_p = E_from_Jr_L(Jr,abs(Lz)+eps*_L0)
        E_pp = E_from_Jr_L(Jr,abs(Lz)+2.0*eps*_L0)

        Frot = _Frot(E,abs(Lz),Lz,alpha)
        Frot_p = _Frot(E_p,abs(Lz)+eps*_L0,Lz,alpha)
        Frot_pp = _Frot(E_pp,abs(Lz)+2.0*eps*_L0,Lz,alpha)

        _, avrDL, _, _, _, _, _, _, avrDLLz = orbitAverageActionCoeffsOptiExact(Jr,abs(Lz),sign(Lz),m_field,alpha,nbAvr,nbK,nbu,m_test)
        _, _, _, _, _, _, _, _, avrDLLz_p = orbitAverageActionCoeffsOptiExact(Jr,abs(Lz)+eps*_L0,Lz/(abs(Lz)+eps*_L0),m_field,alpha,nbAvr,nbK,nbu,m_test)
        _, _, _, _, _, _, _, _, avrDLLz_pp = orbitAverageActionCoeffsOptiExact(Jr,abs(Lz)+2.0*eps*_L0,Lz/(abs(Lz)+2.0*eps*_L0),m_field,alpha,nbAvr,nbK,nbu,m_test)

        dLLz_dL = (-1.5*avrDLLz*Frot + 2.0*avrDLLz_p*Frot_p - 0.5*avrDLLz_pp*Frot_pp)/(eps*_L0)


        Threads.atomic_add!(sumL,avrDL*Frot)
        Threads.atomic_add!(sumLLz_1,avrDLLz*Frot)
        Threads.atomic_add!(sumLLz_2,dLLz_dL)


    end

    sumL[] *= -Jrmax/nbJr
    sumLLz_1[] *= -Jrmax/nbJr
    sumLLz_2[] *= sign(Lz)*Jrmax/nbJr

    return sumL[], sumLLz_1[], sumLLz_2[]
end

function dJrLF_Integrated(Lz::Float64,m_field::Float64,
                        alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                        nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)


    E = E_from_Jr_L(0.0,abs(Lz))
    Frot = _Frot(E,abs(Lz),Lz,alpha)
    _, _, _, _, _, _, dJrL, _ = orbitAverageActionCoeffsOptiExact(0.0,abs(Lz),sign(Lz),m_field,alpha,nbAvr,nbK,nbu,m_test)


    return dJrL*Frot
end

function dLzF_IntegratedPar(Lz::Float64,m_field::Float64, Jrmax::Float64=3.0*_L0, nbJr::Int64=30,
                        DeltaLmax::Float64=3.0*_L0, nbL::Int64=30, eps::Float64=10^(-2),
                        alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                        nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    #sum = Threads.Atomic{Float64}(0.0)


    nbJrL = nbJr*nbL
    tab_JrL = zeros(Float64,2,nbJrL)

    sumLz_1 = Threads.Atomic{Float64}(0.0)
    sumLzLz_1 = Threads.Atomic{Float64}(0.0)

    sumLz_2 = Threads.Atomic{Float64}(0.0)
    sumLzLz_2 = Threads.Atomic{Float64}(0.0)
    sumLzLz_3 = Threads.Atomic{Float64}(0.0)

    iGrid = 1
    for iJr=1:nbJr, iL=1:nbL
        Jr = Jrmax* (iJr-0.5)/nbJr
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        tab_JrL[1,iGrid], tab_JrL[2,iGrid] = Jr, L
        iGrid += 1
    end



    Threads.@threads for iGrid=1:nbJrL
        Jr = tab_JrL[1,iGrid]
        L = tab_JrL[2,iGrid]
        E = E_from_Jr_L(Jr,L)


        Frot = _Frot(E,L,Lz,alpha)
        _, _, avrDLz, _, _, avrDLzLz, _ = orbitAverageActionCoeffsOptiExact(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)



        Threads.atomic_add!(sumLz_1,avrDLz*Frot)
        Threads.atomic_add!(sumLzLz_1,avrDLzLz*Frot)

    end

    sumLz_1[] *= Jrmax/nbJr*DeltaLmax/nbL
    sumLzLz_1[] *= Jrmax/nbJr*DeltaLmax/nbL




    Threads.@threads for iJr=1:nbJr
        Jr = Jrmax*(iJr-0.5)/nbJr
        E = E_from_Jr_L(Jr,abs(Lz))
        E_p = E_from_Jr_L(Jr,abs(Lz)+eps*_L0)
        E_pp = E_from_Jr_L(Jr,abs(Lz)+2.0*eps*_L0)

        Frot = _Frot(E,abs(Lz),Lz,alpha)
        Frot_p = _Frot(E_p,abs(Lz)+eps*_L0,Lz,alpha)
        Frot_pp = _Frot(E_pp,abs(Lz)+2.0*eps*_L0,Lz,alpha)

        _, _, avrDLz, _, _, avrDLzLz, _ = orbitAverageActionCoeffsOptiExact(Jr,abs(Lz),sign(Lz),m_field,alpha,nbAvr,nbK,nbu,m_test)
        _, _, _, _, _, avrDLzLz_p, _ = orbitAverageActionCoeffsOptiExact(Jr,abs(Lz)+eps*_L0,Lz/(abs(Lz)+eps*_L0),m_field,alpha,nbAvr,nbK,nbu,m_test)
        _, _, _, _, _, avrDLzLz_pp, _ = orbitAverageActionCoeffsOptiExact(Jr,abs(Lz)+2.0*eps*_L0,Lz/(abs(Lz)+2.0*eps*_L0),m_field,alpha,nbAvr,nbK,nbu,m_test)

        dLzLz_dL = (-1.5*avrDLzLz*Frot + 2.0*avrDLzLz_p*Frot_p - 0.5*avrDLzLz_pp*Frot_pp)/(eps*_L0)


        Threads.atomic_add!(sumLz_2,avrDLz*Frot)
        Threads.atomic_add!(sumLzLz_2,dLzLz_dL)
        Threads.atomic_add!(sumLzLz_3,avrDLzLz*Frot)



    end

    sumLz_2[] *= sign(Lz)*Jrmax/nbJr
    sumLzLz_2[] *=-Jrmax/nbJr
    sumLzLz_3[] *= 2.0*sign(Lz)*Jrmax/nbJr


    return sumLz_1[], sumLzLz_1[], sumLz_2[], sumLzLz_2[], sumLzLz_3[]
end

function dJrJrF_IntegratedPar(Lz::Float64,m_field::Float64, DeltaLmax::Float64=3.0*_L0, nbL::Int64=30,
                        alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, eps::Float64=10^(-2),
                        nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    sumJrJr = Threads.Atomic{Float64}(0.0)


    Threads.@threads for iL=1:nbL
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL

        E = E_from_Jr_L(0.0,L)
        E_p = E_from_Jr_L(eps*_L0,L)
        E_pp = E_from_Jr_L(2.0*eps*_L0,L)

        Frot = _Frot(E,L,Lz,alpha)
        Frot_p = _Frot(E_p,L,Lz,alpha)
        Frot_pp = _Frot(E_pp,L,Lz,alpha)

        _, _, _, avrDJrJr, _ = orbitAverageActionCoeffsOptiExact(0.0,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        _, _, _, avrDJrJr_p, _ = orbitAverageActionCoeffsOptiExact(eps*_L0,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        _, _, _, avrDJrJr_pp, _ = orbitAverageActionCoeffsOptiExact(2.0*eps*_L0,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

        dJrJr_dJr = (-1.5*avrDJrJr*Frot + 2.0*avrDJrJr_p*Frot_p - 0.5*avrDJrJr_pp*Frot_pp)/(eps*_L0)


        Threads.atomic_add!(sumJrJr,dJrJr_dJr)


    end

    sumJrJr[] *= -DeltaLmax/nbL

    return sumJrJr[]
end

function dLLF_IntegratedPar(Lz::Float64,m_field::Float64, Jrmax::Float64=3.0*_L0, nbJr::Int64=30,
                        alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, eps::Float64=10^(-2),
                        nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    sumLL = Threads.Atomic{Float64}(0.0)


    Threads.@threads for iJr=1:nbJr
        Jr = Jrmax*(iJr-0.5)/nbJr

        E = E_from_Jr_L(Jr,abs(Lz))
        E_p = E_from_Jr_L(Jr,abs(Lz)+eps*_L0)
        E_pp = E_from_Jr_L(Jr,abs(Lz)+2.0*eps*_L0)

        Frot = _Frot(E,abs(Lz),Lz,alpha)
        Frot_p = _Frot(E_p,abs(Lz)+eps*_L0,Lz,alpha)
        Frot_pp = _Frot(E_pp,abs(Lz)+2.0*eps*_L0,Lz,alpha)

        _, _, _, _, avrDLL, _ = orbitAverageActionCoeffsOptiExact(Jr,abs(Lz),Lz/abs(Lz),m_field,alpha,nbAvr,nbK,nbu,m_test)
        _, _, _, _, avrDLL_p, _ = orbitAverageActionCoeffsOptiExact(Jr,abs(Lz)+eps*_L0,Lz/(abs(Lz)+eps*_L0),m_field,alpha,nbAvr,nbK,nbu,m_test)
        _, _, _, _, avrDLL_pp, _ = orbitAverageActionCoeffsOptiExact(Jr,abs(Lz)+2.0*eps*_L0,Lz/(abs(Lz)+2.0*eps*_L0),m_field,alpha,nbAvr,nbK,nbu,m_test)

        dLL_dL = (-1.5*avrDLL*Frot + 2.0*avrDLL_p*Frot_p - 0.5*avrDLL_pp*Frot_pp)/(eps*_L0)


        Threads.atomic_add!(sumLL,dLL_dL)


    end

    sumLL[] *= -Jrmax/nbJr

    return sumLL[]
end
