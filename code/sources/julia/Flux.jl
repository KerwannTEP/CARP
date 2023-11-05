function flux(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    # F_i = D_i Frot - 0.5 * d/dJk [D_ik Frot]

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffs(Jr,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot = _Frot(E,L,Lz,alpha)


    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffs(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffs(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffs(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffs(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffs(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffs(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzm =_Frot(E,L,Lz_m,alpha)


    # Jr-component

    DJrF = dJr*Frot

    dJrJrF_Jr = (dJrJr_Jrp*Frot_Jrp-dJrJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dJrLF_L = (dJrL_Lp*Frot_Lp-dJrL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dJrLzF_Lz = (dJrLz_Lzp*Frot_Lzp-dJrLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxJr = DJrF - 0.5*(dJrJrF_Jr+dJrLF_L+dJrLzF_Lz)

    # L-component

    DLF = dL*Frot

    dJrLF_Jr = (dJrL_Jrp*Frot_Jrp-dJrL_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLF_L = (dLL_Lp*Frot_Lp-dLL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLLzF_Lz = (dLLz_Lzp*Frot_Lzp-dLLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxL = DLF - 0.5*(dJrLF_Jr+dLLF_L+dLLzF_Lz)

    # L-component

    DLzF = dLz*Frot

    dJrLzF_Jr = (dJrLz_Jrp*Frot_Jrp-dJrLz_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLzF_L = (dLLz_Lp*Frot_Lp-dLLz_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLzLzF_Lz = (dLzLz_Lzp*Frot_Lzp-dLzLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxLz = DLzF - 0.5*(dJrLzF_Jr+dLLzF_L+dLzLzF_Lz)

    return fluxJr, fluxL, fluxLz
end

function fluxPar(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    # F_i = D_i Frot - 0.5 * d/dJk [D_ik Frot]

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot = _Frot(E,L,Lz,alpha)


    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzm =_Frot(E,L,Lz_m,alpha)


    # Jr-component

    DJrF = dJr*Frot

    dJrJrF_Jr = (dJrJr_Jrp*Frot_Jrp-dJrJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dJrLF_L = (dJrL_Lp*Frot_Lp-dJrL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dJrLzF_Lz = (dJrLz_Lzp*Frot_Lzp-dJrLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxJr = DJrF - 0.5*(dJrJrF_Jr+dJrLF_L+dJrLzF_Lz)

    # L-component

    DLF = dL*Frot

    dJrLF_Jr = (dJrL_Jrp*Frot_Jrp-dJrL_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLF_L = (dLL_Lp*Frot_Lp-dLL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLLzF_Lz = (dLLz_Lzp*Frot_Lzp-dLLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxL = DLF - 0.5*(dJrLF_Jr+dLLF_L+dLLzF_Lz)

    # L-component

    DLzF = dLz*Frot

    dJrLzF_Jr = (dJrLz_Jrp*Frot_Jrp-dJrLz_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLzF_L = (dLLz_Lp*Frot_Lp-dLLz_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLzLzF_Lz = (dLzLz_Lzp*Frot_Lzp-dLzLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxLz = DLzF - 0.5*(dJrLzF_Jr+dLLzF_L+dLzLzF_Lz)

    return fluxJr, fluxL, fluxLz
end

function fluxOpti(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    # F_i = D_i Frot - 0.5 * d/dJk [D_ik Frot]

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOpti(Jr,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot = _Frot(E,L,Lz,alpha)


    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOpti(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOpti(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOpti(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOpti(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOpti(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOpti(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzm =_Frot(E,L,Lz_m,alpha)


    # Jr-component

    DJrF = dJr*Frot

    dJrJrF_Jr = (dJrJr_Jrp*Frot_Jrp-dJrJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dJrLF_L = (dJrL_Lp*Frot_Lp-dJrL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dJrLzF_Lz = (dJrLz_Lzp*Frot_Lzp-dJrLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxJr = DJrF - 0.5*(dJrJrF_Jr+dJrLF_L+dJrLzF_Lz)

    # L-component

    DLF = dL*Frot

    dJrLF_Jr = (dJrL_Jrp*Frot_Jrp-dJrL_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLF_L = (dLL_Lp*Frot_Lp-dLL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLLzF_Lz = (dLLz_Lzp*Frot_Lzp-dLLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxL = DLF - 0.5*(dJrLF_Jr+dLLF_L+dLLzF_Lz)

    # L-component

    DLzF = dLz*Frot

    dJrLzF_Jr = (dJrLz_Jrp*Frot_Jrp-dJrLz_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLzF_L = (dLLz_Lp*Frot_Lp-dLLz_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLzLzF_Lz = (dLzLz_Lzp*Frot_Lzp-dLzLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxLz = DLzF - 0.5*(dJrLzF_Jr+dLLzF_L+dLzLzF_Lz)

    return fluxJr, fluxL, fluxLz
end

function fluxOptiPar(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    # F_i = D_i Frot - 0.5 * d/dJk [D_ik Frot]

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot = _Frot(E,L,Lz,alpha)


    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzm =_Frot(E,L,Lz_m,alpha)


    # Jr-component

    DJrF = dJr*Frot

    dJrJrF_Jr = (dJrJr_Jrp*Frot_Jrp-dJrJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dJrLF_L = (dJrL_Lp*Frot_Lp-dJrL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dJrLzF_Lz = (dJrLz_Lzp*Frot_Lzp-dJrLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxJr = DJrF - 0.5*(dJrJrF_Jr+dJrLF_L+dJrLzF_Lz)

    # L-component

    DLF = dL*Frot

    dJrLF_Jr = (dJrL_Jrp*Frot_Jrp-dJrL_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLF_L = (dLL_Lp*Frot_Lp-dLL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLLzF_Lz = (dLLz_Lzp*Frot_Lzp-dLLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxL = DLF - 0.5*(dJrLF_Jr+dLLF_L+dLLzF_Lz)

    # L-component

    DLzF = dLz*Frot

    dJrLzF_Jr = (dJrLz_Jrp*Frot_Jrp-dJrLz_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLzF_L = (dLLz_Lp*Frot_Lp-dLLz_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLzLzF_Lz = (dLzLz_Lzp*Frot_Lzp-dLzLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxLz = DLzF - 0.5*(dJrLzF_Jr+dLLzF_L+dLzLzF_Lz)

    return fluxJr, fluxL, fluxLz
end


function fluxOptiExactSign(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    # F_i = D_i Frot - 0.5 * d/dJk [D_ik Frot]

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExact(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    Frot = _Frot(E,L,Lz,alpha)


    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lzm =_Frot(E,L,Lz_m,alpha)


    # Jr-component

    DJrF = dJr*Frot

    dJrJrF_Jr = (dJrJr_Jrp*Frot_Jrp-dJrJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dJrLF_L = (dJrL_Lp*Frot_Lp-dJrL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dJrLzF_Lz = (dJrLz_Lzp*Frot_Lzp-dJrLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxJr = DJrF - 0.5*(dJrJrF_Jr+dJrLF_L+dJrLzF_Lz)

    # L-component

    DLF = dL*Frot

    dJrLF_Jr = (dJrL_Jrp*Frot_Jrp-dJrL_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLF_L = (dLL_Lp*Frot_Lp-dLL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLLzF_Lz = (dLLz_Lzp*Frot_Lzp-dLLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxL = DLF - 0.5*(dJrLF_Jr+dLLF_L+dLLzF_Lz)

    # L-component

    DLzF = dLz*Frot

    dJrLzF_Jr = (dJrLz_Jrp*Frot_Jrp-dJrLz_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLzF_L = (dLLz_Lp*Frot_Lp-dLLz_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLzLzF_Lz = (dLzLz_Lzp*Frot_Lzp-dLzLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxLz = DLzF - 0.5*(dJrLzF_Jr+dLLzF_L+dLzLzF_Lz)

    return fluxJr, fluxL, fluxLz
end


function fluxOptiExactSignPar(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    # F_i = D_i Frot - 0.5 * d/dJk [D_ik Frot]

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    Frot = _Frot(E,L,Lz,alpha)


    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lzm =_Frot(E,L,Lz_m,alpha)


    # Jr-component

    DJrF = dJr*Frot

    dJrJrF_Jr = (dJrJr_Jrp*Frot_Jrp-dJrJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dJrLF_L = (dJrL_Lp*Frot_Lp-dJrL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dJrLzF_Lz = (dJrLz_Lzp*Frot_Lzp-dJrLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxJr = DJrF - 0.5*(dJrJrF_Jr+dJrLF_L+dJrLzF_Lz)

    # L-component

    DLF = dL*Frot

    dJrLF_Jr = (dJrL_Jrp*Frot_Jrp-dJrL_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLF_L = (dLL_Lp*Frot_Lp-dLL_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLLzF_Lz = (dLLz_Lzp*Frot_Lzp-dLLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxL = DLF - 0.5*(dJrLF_Jr+dLLF_L+dLLzF_Lz)

    # L-component

    DLzF = dLz*Frot

    dJrLzF_Jr = (dJrLz_Jrp*Frot_Jrp-dJrLz_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLzF_L = (dLLz_Lp*Frot_Lp-dLLz_Lm*Frot_Lm)/(2.0*eps*_L0)
    dLzLzF_Lz = (dLzLz_Lzp*Frot_Lzp-dLzLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    fluxLz = DLzF - 0.5*(dJrLzF_Jr+dLLzF_L+dLzLzF_Lz)

    return fluxJr, fluxL, fluxLz
end



# USED
function fluxOptiExactSign_JrL(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default,
            nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
            nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    # F_i = D_i Frot - 0.5 * d/dJk [D_ik Frot]

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0



    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExact_nbKseparate(Jr,L,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)


    Frot = _Frot_cosI(E,L,cosI,alpha)


    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExact_nbKseparate(Jr_p,L,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)
    Frot_Jrp = _Frot_cosI(E_Jrp,L,cosI,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExact_nbKseparate(Jr_m,L,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)
    Frot_Jrm = _Frot_cosI(E_Jrm,L,cosI,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExact_nbKseparate(Jr,L_p,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)
    Frot_Lp = _Frot_cosI(E_Lp,L_p,cosI,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExact_nbKseparate(Jr,L_m,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)
    Frot_Lm = _Frot_cosI(E_Lm,L_m,cosI,alpha)

    # Jr-component

    DJrF = dJr*Frot

    dJrJrF_Jr = (dJrJr_Jrp*Frot_Jrp-dJrJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dJrLF_L = (dJrL_Lp*Frot_Lp-dJrL_Lm*Frot_Lm)/(2.0*eps*_L0)

    fluxJr = DJrF - 0.5*(dJrJrF_Jr+dJrLF_L)

    # L-component

    DLF = dL*Frot

    dJrLF_Jr = (dJrL_Jrp*Frot_Jrp-dJrL_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    dLLF_L = (dLL_Lp*Frot_Lp-dLL_Lm*Frot_Lm)/(2.0*eps*_L0)

    fluxL = DLF - 0.5*(dJrLF_Jr+dLLF_L)

    return fluxJr, fluxL
end


function dFdtOptiExactSign_2D_JrL(Jr::Float64, L::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbCosI::Int64=50, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default,
            nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
            nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    sumJr = 0.0
    sumL = 0.0

    for i=1:nbCosI
        cosI = -1.0 + 2.0/nbCosI*(i-0.5)

        fJr_p, _ = fluxOptiExactSign_JrL(Jr+eps,L,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,eps,m_test)
        fJr_m, _ = fluxOptiExactSign_JrL(Jr-eps,L,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,eps,m_test)

        dJr = (fJr_p-fJr_m)/(2.0*eps)

        _, fL_p  = fluxOptiExactSign_JrL(Jr,L+eps,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,eps,m_test)
        _, fL_m  = fluxOptiExactSign_JrL(Jr,L-eps,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,eps,m_test)

        dL = (fL_p-fL_m)/(2.0*eps)

        sumJr += dJr
        sumL  += dL
    end

    sumJr *= 2.0/nbCosI
    sumL *= 2.0/nbCosI

    return -(sumJr+sumL)
end

# USED
function dFdtOptiExactSign_2D_JrcosI(Jr::Float64, cosI::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbL::Int64=50, Lmax::Float64=3.0, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default,
            nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
            nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    sumJr = 0.0
    sumcosI = 0.0

    for i=1:nbL
        L = Lmax/nbL*(i-0.5)

        

        fJr_p, _ = fluxOptiExactSign_JrL(Jr+eps,L,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,eps,m_test)
        fJr_m, _ = fluxOptiExactSign_JrL(Jr-eps,L,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,eps,m_test)

        dJr = (fJr_p-fJr_m)/(2.0*eps)

        fcosI_p  = FluxCosIOptiExact_nbKseparate(Jr,L,cosI+eps,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)
        fcosI_m  = FluxCosIOptiExact_nbKseparate(Jr,L,cosI-eps,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)


        dcosI = (fcosI_p-fcosI_m)/(2.0*eps)

        sumJr += dJr
        sumcosI  += dcosI
    end


    sumJr *= Lmax/nbL
    sumcosI *= Lmax/nbL

    return -(sumJr+sumcosI)
end

function dFdtOptiExactSign_2D_LcosI(L::Float64, cosI::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbJr::Int64=50, Jrmax::Float64=3.0, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default,
            nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
            nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    sumL = 0.0
    sumcosI = 0.0

    for i=1:nbJr
        Jr = Jrmax/nbJr*(i-0.5)

        _ , fL_p = fluxOptiExactSign_JrL(Jr,L+eps,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,eps,m_test)
        _ , fL_m = fluxOptiExactSign_JrL(Jr,L-eps,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,eps,m_test)

        dL = (fL_p-fL_m)/(2.0*eps)

        fcosI_p  = FluxCosIOptiExact_nbKseparate(Jr,L,cosI+eps,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)
        fcosI_m  = FluxCosIOptiExact_nbKseparate(Jr,L,cosI-eps,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)


        dcosI = (fcosI_p-fcosI_m)/(2.0*eps)

        sumL += dL
        sumcosI  += dcosI
    end

    sumL *= Jrmax/nbJr
    sumcosI *= Jrmax/nbJr

    return -(sumL+sumcosI)
end

###########################################################
###########################################################


function dFdt(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    # Value at point

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffs(Jr,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot= _Frot(E,L,Lz,alpha)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffs(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffs(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffs(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffs(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffs(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffs(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzm = _Frot(E,L,Lz_m,alpha)


    DJrF_Jr = (dJr_Jrp*Frot_Jrp-dJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    DLF_L = (dL_Lp*Frot_Lp-dL_Lm*Frot_Lm)/(2.0*eps*_L0)
    DLzF_Lz = (dLz_Lzp*Frot_Lzp-dLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    DJrJr_F_JrJr = (dJrJr_Jrp*Frot_Jrp + dJrJr_Jrm*Frot_Jrm - 2.0*dJrJr*Frot)/(eps*_L0)^2
    DLL_F_LL = (dLL_Lp*Frot_Lp + dLL_Lm*Frot_Lm - 2.0*dLL*Frot)/(eps*_L0)^2
    DLzLz_F_LzLz = (dLzLz_Lzp*Frot_Lzp + dLzLz_Lzm*Frot_Lzm - 2.0*dLzLz*Frot)/(eps*_L0)^2





    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffs(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLp = _Frot(E_Jrp_Lp,L_p,Lz,alpha)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffs(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLm = _Frot(E_Jrm_Lm,L_m,Lz,alpha)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffs(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLm = _Frot(E_Jrp_Lm,L_m,Lz,alpha)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffs(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLp = _Frot(E_Jrm_Lp,L_p,Lz,alpha)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffs(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLzp = _Frot(E_Jrp,L,Lz_p,alpha)
    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffs(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLzm = _Frot(E_Jrm,L,Lz_m,alpha)

    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffs(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLzm = _Frot(E_Jrp,L,Lz_m,alpha)
    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffs(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLzp = _Frot(E_Jrm,L,Lz_p,alpha)

    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffs(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LpLzp = _Frot(E_Lp,L_p,Lz_p,alpha)
    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffs(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LmLzm = _Frot(E_Lm,L_m,Lz_m,alpha)

    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffs(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LpLzm = _Frot(E_Lp,L_p,Lz_m,alpha)
    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffs(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LmLzp =_Frot(E_Lm,L_m,Lz_p,alpha)




    DJrL_F_JrL = (dJrL_JrpLp*Frot_JrpLp+dJrL_JrmLm*Frot_JrmLm-dJrL_JrpLm*Frot_JrpLm-dJrL_JrmLp*Frot_JrmLp)/(2.0*eps*_L0)^2

    DJrLz_F_JrLz = (dJrLz_JrpLzp*Frot_JrpLzp+dJrLz_JrmLzm*Frot_JrmLzm-dJrLz_JrpLzm*Frot_JrpLzm-dJrLz_JrmLzp*Frot_JrmLzp)/(2.0*eps*_L0)^2

    DLLz_F_LLz = (dLLz_LpLzp*Frot_LpLzp+dLLz_LmLzm*Frot_LmLzm-dLLz_LpLzm*Frot_LpLzm-dLLz_LmLzp*Frot_LmLzp)/(2.0*eps*_L0)^2


    fluxJr_Jr = DJrF_Jr - 0.5*(DJrJr_F_JrJr + DJrL_F_JrL + DJrLz_F_JrLz)
    fluxL_L = DLF_L - 0.5*(DJrL_F_JrL + DLL_F_LL + DLLz_F_LLz)
    fluxLz_Lz = DLzF_Lz - 0.5*(DJrLz_F_JrLz + DLLz_F_LLz + DLzLz_F_LzLz)

    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz)

end

function dFdtOpti(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    # Value at point

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOpti(Jr,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot= _Frot(E,L,Lz,alpha)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOpti(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOpti(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOpti(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOpti(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOpti(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOpti(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzm = _Frot(E,L,Lz_m,alpha)


    DJrF_Jr = (dJr_Jrp*Frot_Jrp-dJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    DLF_L = (dL_Lp*Frot_Lp-dL_Lm*Frot_Lm)/(2.0*eps*_L0)
    DLzF_Lz = (dLz_Lzp*Frot_Lzp-dLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    DJrJr_F_JrJr = (dJrJr_Jrp*Frot_Jrp + dJrJr_Jrm*Frot_Jrm - 2.0*dJrJr*Frot)/(eps*_L0)^2
    DLL_F_LL = (dLL_Lp*Frot_Lp + dLL_Lm*Frot_Lm - 2.0*dLL*Frot)/(eps*_L0)^2
    DLzLz_F_LzLz = (dLzLz_Lzp*Frot_Lzp + dLzLz_Lzm*Frot_Lzm - 2.0*dLzLz*Frot)/(eps*_L0)^2





    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsOpti(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLp = _Frot(E_Jrp_Lp,L_p,Lz,alpha)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsOpti(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLm = _Frot(E_Jrm_Lm,L_m,Lz,alpha)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsOpti(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLm = _Frot(E_Jrp_Lm,L_m,Lz,alpha)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsOpti(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLp = _Frot(E_Jrm_Lp,L_p,Lz,alpha)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsOpti(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLzp = _Frot(E_Jrp,L,Lz_p,alpha)
    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsOpti(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLzm = _Frot(E_Jrm,L,Lz_m,alpha)

    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsOpti(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLzm = _Frot(E_Jrp,L,Lz_m,alpha)
    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsOpti(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLzp = _Frot(E_Jrm,L,Lz_p,alpha)

    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsOpti(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LpLzp = _Frot(E_Lp,L_p,Lz_p,alpha)
    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsOpti(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LmLzm = _Frot(E_Lm,L_m,Lz_m,alpha)

    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsOpti(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LpLzm = _Frot(E_Lp,L_p,Lz_m,alpha)
    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsOpti(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LmLzp =_Frot(E_Lm,L_m,Lz_p,alpha)




    DJrL_F_JrL = (dJrL_JrpLp*Frot_JrpLp+dJrL_JrmLm*Frot_JrmLm-dJrL_JrpLm*Frot_JrpLm-dJrL_JrmLp*Frot_JrmLp)/(2.0*eps*_L0)^2

    DJrLz_F_JrLz = (dJrLz_JrpLzp*Frot_JrpLzp+dJrLz_JrmLzm*Frot_JrmLzm-dJrLz_JrpLzm*Frot_JrpLzm-dJrLz_JrmLzp*Frot_JrmLzp)/(2.0*eps*_L0)^2

    DLLz_F_LLz = (dLLz_LpLzp*Frot_LpLzp+dLLz_LmLzm*Frot_LmLzm-dLLz_LpLzm*Frot_LpLzm-dLLz_LmLzp*Frot_LmLzp)/(2.0*eps*_L0)^2


    fluxJr_Jr = DJrF_Jr - 0.5*(DJrJr_F_JrJr + DJrL_F_JrL + DJrLz_F_JrLz)
    fluxL_L = DLF_L - 0.5*(DJrL_F_JrL + DLL_F_LL + DLLz_F_LLz)
    fluxLz_Lz = DLzF_Lz - 0.5*(DJrLz_F_JrLz + DLLz_F_LLz + DLzLz_F_LzLz)

    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz)

end

function dFdtOptiExact(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    # Value at point

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExact(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot= _Frot(E,L,Lz,alpha)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lzm = _Frot(E,L,Lz_m,alpha)


    DJrF_Jr = (dJr_Jrp*Frot_Jrp-dJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    DLF_L = (dL_Lp*Frot_Lp-dL_Lm*Frot_Lm)/(2.0*eps*_L0)
    DLzF_Lz = (dLz_Lzp*Frot_Lzp-dLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    DJrJr_F_JrJr = (dJrJr_Jrp*Frot_Jrp + dJrJr_Jrm*Frot_Jrm - 2.0*dJrJr*Frot)/(eps*_L0)^2
    DLL_F_LL = (dLL_Lp*Frot_Lp + dLL_Lm*Frot_Lm - 2.0*dLL*Frot)/(eps*_L0)^2
    DLzLz_F_LzLz = (dLzLz_Lzp*Frot_Lzp + dLzLz_Lzm*Frot_Lzm - 2.0*dLzLz*Frot)/(eps*_L0)^2





    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsOptiExact(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrpLp = _Frot(E_Jrp_Lp,L_p,Lz,alpha)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsOptiExact(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrmLm = _Frot(E_Jrm_Lm,L_m,Lz,alpha)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsOptiExact(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrpLm = _Frot(E_Jrp_Lm,L_m,Lz,alpha)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsOptiExact(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrmLp = _Frot(E_Jrm_Lp,L_p,Lz,alpha)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrpLzp = _Frot(E_Jrp,L,Lz_p,alpha)
    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrmLzm = _Frot(E_Jrm,L,Lz_m,alpha)

    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrpLzm = _Frot(E_Jrp,L,Lz_m,alpha)
    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrmLzp = _Frot(E_Jrm,L,Lz_p,alpha)

    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_LpLzp = _Frot(E_Lp,L_p,Lz_p,alpha)
    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_LmLzm = _Frot(E_Lm,L_m,Lz_m,alpha)

    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_LpLzm = _Frot(E_Lp,L_p,Lz_m,alpha)
    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_LmLzp =_Frot(E_Lm,L_m,Lz_p,alpha)




    DJrL_F_JrL = (dJrL_JrpLp*Frot_JrpLp+dJrL_JrmLm*Frot_JrmLm-dJrL_JrpLm*Frot_JrpLm-dJrL_JrmLp*Frot_JrmLp)/(2.0*eps*_L0)^2

    DJrLz_F_JrLz = (dJrLz_JrpLzp*Frot_JrpLzp+dJrLz_JrmLzm*Frot_JrmLzm-dJrLz_JrpLzm*Frot_JrpLzm-dJrLz_JrmLzp*Frot_JrmLzp)/(2.0*eps*_L0)^2

    DLLz_F_LLz = (dLLz_LpLzp*Frot_LpLzp+dLLz_LmLzm*Frot_LmLzm-dLLz_LpLzm*Frot_LpLzm-dLLz_LmLzp*Frot_LmLzp)/(2.0*eps*_L0)^2


    fluxJr_Jr = DJrF_Jr - 0.5*(DJrJr_F_JrJr + DJrL_F_JrL + DJrLz_F_JrLz)
    fluxL_L = DLF_L - 0.5*(DJrL_F_JrL + DLL_F_LL + DLLz_F_LLz)
    fluxLz_Lz = DLzF_Lz - 0.5*(DJrLz_F_JrLz + DLLz_F_LLz + DLzLz_F_LzLz)

    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz)

end

function dFdtOptiExactPar(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    # Value at point

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot= _Frot(E,L,Lz,alpha)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_Lzm = _Frot(E,L,Lz_m,alpha)


    DJrF_Jr = (dJr_Jrp*Frot_Jrp-dJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    DLF_L = (dL_Lp*Frot_Lp-dL_Lm*Frot_Lm)/(2.0*eps*_L0)
    DLzF_Lz = (dLz_Lzp*Frot_Lzp-dLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    DJrJr_F_JrJr = (dJrJr_Jrp*Frot_Jrp + dJrJr_Jrm*Frot_Jrm - 2.0*dJrJr*Frot)/(eps*_L0)^2
    DLL_F_LL = (dLL_Lp*Frot_Lp + dLL_Lm*Frot_Lm - 2.0*dLL*Frot)/(eps*_L0)^2
    DLzLz_F_LzLz = (dLzLz_Lzp*Frot_Lzp + dLzLz_Lzm*Frot_Lzm - 2.0*dLzLz*Frot)/(eps*_L0)^2





    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrpLp = _Frot(E_Jrp_Lp,L_p,Lz,alpha)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrmLm = _Frot(E_Jrm_Lm,L_m,Lz,alpha)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsOptiExactPar(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrpLm = _Frot(E_Jrp_Lm,L_m,Lz,alpha)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsOptiExactPar(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrmLp = _Frot(E_Jrm_Lp,L_p,Lz,alpha)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrpLzp = _Frot(E_Jrp,L,Lz_p,alpha)
    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrmLzm = _Frot(E_Jrm,L,Lz_m,alpha)

    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrpLzm = _Frot(E_Jrp,L,Lz_m,alpha)
    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_JrmLzp = _Frot(E_Jrm,L,Lz_p,alpha)

    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_LpLzp = _Frot(E_Lp,L_p,Lz_p,alpha)
    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_LmLzm = _Frot(E_Lm,L_m,Lz_m,alpha)

    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_LpLzm = _Frot(E_Lp,L_p,Lz_m,alpha)
    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Frot_LmLzp =_Frot(E_Lm,L_m,Lz_p,alpha)




    DJrL_F_JrL = (dJrL_JrpLp*Frot_JrpLp+dJrL_JrmLm*Frot_JrmLm-dJrL_JrpLm*Frot_JrpLm-dJrL_JrmLp*Frot_JrmLp)/(2.0*eps*_L0)^2

    DJrLz_F_JrLz = (dJrLz_JrpLzp*Frot_JrpLzp+dJrLz_JrmLzm*Frot_JrmLzm-dJrLz_JrpLzm*Frot_JrpLzm-dJrLz_JrmLzp*Frot_JrmLzp)/(2.0*eps*_L0)^2

    DLLz_F_LLz = (dLLz_LpLzp*Frot_LpLzp+dLLz_LmLzm*Frot_LmLzm-dLLz_LpLzm*Frot_LpLzm-dLLz_LmLzp*Frot_LmLzp)/(2.0*eps*_L0)^2


    fluxJr_Jr = DJrF_Jr - 0.5*(DJrJr_F_JrJr + DJrL_F_JrL + DJrLz_F_JrLz)
    fluxL_L = DLF_L - 0.5*(DJrL_F_JrL + DLL_F_LL + DLLz_F_LLz)
    fluxLz_Lz = DLzF_Lz - 0.5*(DJrLz_F_JrLz + DLLz_F_LLz + DLzLz_F_LzLz)

    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz)

end



function dFdt_Par(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    # Value at point

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot= _Frot(E,L,Lz,alpha)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzm = _Frot(E,L,Lz_m,alpha)


    DJrF_Jr = (dJr_Jrp*Frot_Jrp-dJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    DLF_L = (dL_Lp*Frot_Lp-dL_Lm*Frot_Lm)/(2.0*eps*_L0)
    DLzF_Lz = (dLz_Lzp*Frot_Lzp-dLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    DJrJr_F_JrJr = (dJrJr_Jrp*Frot_Jrp + dJrJr_Jrm*Frot_Jrm - 2.0*dJrJr*Frot)/(eps*_L0)^2
    DLL_F_LL = (dLL_Lp*Frot_Lp + dLL_Lm*Frot_Lm - 2.0*dLL*Frot)/(eps*_L0)^2
    DLzLz_F_LzLz = (dLzLz_Lzp*Frot_Lzp + dLzLz_Lzm*Frot_Lzm - 2.0*dLzLz*Frot)/(eps*_L0)^2





    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsPar(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLp = _Frot(E_Jrp_Lp,L_p,Lz,alpha)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsPar(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLm = _Frot(E_Jrm_Lm,L_m,Lz,alpha)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsPar(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLm = _Frot(E_Jrp_Lm,L_m,Lz,alpha)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsPar(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLp = _Frot(E_Jrm_Lp,L_p,Lz,alpha)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsPar(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLzp = _Frot(E_Jrp,L,Lz_p,alpha)
    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsPar(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLzm = _Frot(E_Jrm,L,Lz_m,alpha)

    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsPar(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLzm = _Frot(E_Jrp,L,Lz_m,alpha)
    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsPar(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLzp = _Frot(E_Jrm,L,Lz_p,alpha)

    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsPar(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LpLzp = _Frot(E_Lp,L_p,Lz_p,alpha)
    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsPar(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LmLzm = _Frot(E_Lm,L_m,Lz_m,alpha)

    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsPar(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LpLzm = _Frot(E_Lp,L_p,Lz_m,alpha)
    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsPar(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LmLzp =_Frot(E_Lm,L_m,Lz_p,alpha)




    DJrL_F_JrL = (dJrL_JrpLp*Frot_JrpLp+dJrL_JrmLm*Frot_JrmLm-dJrL_JrpLm*Frot_JrpLm-dJrL_JrmLp*Frot_JrmLp)/(2.0*eps*_L0)^2

    DJrLz_F_JrLz = (dJrLz_JrpLzp*Frot_JrpLzp+dJrLz_JrmLzm*Frot_JrmLzm-dJrLz_JrpLzm*Frot_JrpLzm-dJrLz_JrmLzp*Frot_JrmLzp)/(2.0*eps*_L0)^2

    DLLz_F_LLz = (dLLz_LpLzp*Frot_LpLzp+dLLz_LmLzm*Frot_LmLzm-dLLz_LpLzm*Frot_LpLzm-dLLz_LmLzp*Frot_LmLzp)/(2.0*eps*_L0)^2


    fluxJr_Jr = DJrF_Jr - 0.5*(DJrJr_F_JrJr + DJrL_F_JrL + DJrLz_F_JrLz)
    fluxL_L = DLF_L - 0.5*(DJrL_F_JrL + DLL_F_LL + DLLz_F_LLz)
    fluxLz_Lz = DLzF_Lz - 0.5*(DJrLz_F_JrLz + DLLz_F_LLz + DLzLz_F_LzLz)

    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz)

end

function dFdtOpti_Par(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)

    # Value at point

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot= _Frot(E,L,Lz,alpha)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_Lzm = _Frot(E,L,Lz_m,alpha)


    DJrF_Jr = (dJr_Jrp*Frot_Jrp-dJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    DLF_L = (dL_Lp*Frot_Lp-dL_Lm*Frot_Lm)/(2.0*eps*_L0)
    DLzF_Lz = (dLz_Lzp*Frot_Lzp-dLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    DJrJr_F_JrJr = (dJrJr_Jrp*Frot_Jrp + dJrJr_Jrm*Frot_Jrm - 2.0*dJrJr*Frot)/(eps*_L0)^2
    DLL_F_LL = (dLL_Lp*Frot_Lp + dLL_Lm*Frot_Lm - 2.0*dLL*Frot)/(eps*_L0)^2
    DLzLz_F_LzLz = (dLzLz_Lzp*Frot_Lzp + dLzLz_Lzm*Frot_Lzm - 2.0*dLzLz*Frot)/(eps*_L0)^2





    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsOptiPar(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLp = _Frot(E_Jrp_Lp,L_p,Lz,alpha)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsOptiPar(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLm = _Frot(E_Jrm_Lm,L_m,Lz,alpha)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsOptiPar(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLm = _Frot(E_Jrp_Lm,L_m,Lz,alpha)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsOptiPar(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLp = _Frot(E_Jrm_Lp,L_p,Lz,alpha)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsOptiPar(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLzp = _Frot(E_Jrp,L,Lz_p,alpha)
    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsOptiPar(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLzm = _Frot(E_Jrm,L,Lz_m,alpha)

    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsOptiPar(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrpLzm = _Frot(E_Jrp,L,Lz_m,alpha)
    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsOptiPar(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_JrmLzp = _Frot(E_Jrm,L,Lz_p,alpha)

    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsOptiPar(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LpLzp = _Frot(E_Lp,L_p,Lz_p,alpha)
    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsOptiPar(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LmLzm = _Frot(E_Lm,L_m,Lz_m,alpha)

    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsOptiPar(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LpLzm = _Frot(E_Lp,L_p,Lz_m,alpha)
    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsOptiPar(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)
    Frot_LmLzp =_Frot(E_Lm,L_m,Lz_p,alpha)




    DJrL_F_JrL = (dJrL_JrpLp*Frot_JrpLp+dJrL_JrmLm*Frot_JrmLm-dJrL_JrpLm*Frot_JrpLm-dJrL_JrmLp*Frot_JrmLp)/(2.0*eps*_L0)^2

    DJrLz_F_JrLz = (dJrLz_JrpLzp*Frot_JrpLzp+dJrLz_JrmLzm*Frot_JrmLzm-dJrLz_JrpLzm*Frot_JrpLzm-dJrLz_JrmLzp*Frot_JrmLzp)/(2.0*eps*_L0)^2

    DLLz_F_LLz = (dLLz_LpLzp*Frot_LpLzp+dLLz_LmLzm*Frot_LmLzm-dLLz_LpLzm*Frot_LpLzm-dLLz_LmLzp*Frot_LmLzp)/(2.0*eps*_L0)^2


    fluxJr_Jr = DJrF_Jr - 0.5*(DJrJr_F_JrJr + DJrL_F_JrL + DJrLz_F_JrLz)
    fluxL_L = DLF_L - 0.5*(DJrL_F_JrL + DLL_F_LL + DLLz_F_LLz)
    fluxLz_Lz = DLzF_Lz - 0.5*(DJrLz_F_JrLz + DLLz_F_LLz + DLzLz_F_LzLz)

    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz)

end


##################################################################
# dF/dt 3D for g(x)=sign(x)
# Consider dirac(x) = 0 for x != 0
##################################################################



function dFdtOptiExactSign(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0

    epsL = min(eps,0.1*abs(1.0-abs(Lz/L)),0.1*abs(1.0-abs(Lz/L)))

    L_p = L + epsL*_L0
    L_m = L - epsL*_L0

    #println((eps,0.5*abs(L_m-Lz),0.5*abs(L_m+Lz)))



    Lz_p = Lz + epsL*_L0
    Lz_m = Lz - epsL*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)



    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExact(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot= _F(E,L)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrp = _F(E_Jrp,L)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrm = _F(E_Jrm,L)



    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lp = _F(E_Lp,L_p)

    #println((Lz,L_m,Lz/L_m))

    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lm = _F(E_Lm,L_m)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)







    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsOptiExact(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLp = _F(E_Jrp_Lp,L_p)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsOptiExact(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLm = _F(E_Jrm_Lm,L_m)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsOptiExact(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLm = _F(E_Jrp_Lm,L_m)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsOptiExact(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLp = _F(E_Jrm_Lp,L_p)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)



    # Jr component

    dJr_Jr = (dJr_Jrp-dJr_Jrm)/(2.0*_L0*eps)
    dJrJr_JrJr = (dJrJr_Jrp+dJrJr_Jrm-2.0*dJrJr)/(eps*_L0)^2
    dJrL_JrL = (dJrL_JrpLp+dJrL_JrmLm-dJrL_JrpLm-dJrL_JrmLp)/(2.0*_L0*eps)^2
    dJrLz_JrLz = (dJrLz_JrpLzp+dJrLz_JrmLzm-dJrLz_JrpLzm-dJrLz_JrmLzp)/(2.0*_L0*eps)^2

    dJrJr_Jr = (dJrJr_Jrp-dJrJr_Jrm)/(2.0*_L0*eps)
    dJrL_L = (dJrL_Lp-dJrL_Lm)/(2.0*_L0*eps)
    dJrLz_Lz = (dJrLz_Lzp-dJrLz_Lzm)/(2.0*_L0*eps)
    Ftot_Jr = (Ftot_Jrp-Ftot_Jrm)/(2.0*_L0*eps)

    dJrL_Jr = (dJrL_Jrp-dJrL_Jrm)/(2.0*_L0*eps)
    dJrLz_Jr = (dJrLz_Jrp-dJrLz_Jrm)/(2.0*_L0*eps)
    Ftot_L = (Ftot_Lp-Ftot_Lm)/(2.0*_L0*eps)

    Ftot_JrJr = (Ftot_Jrp+Ftot_Jrm-2.0*Ftot)/(_L0*eps)^2
    Ftot_JrL = (Ftot_JrpLp+Ftot_JrmLm-Ftot_JrpLm-Ftot_JrmLp)/(2.0*_L0*eps)^2

    fluxJr_Jr = ((dJr_Jr-0.5*(dJrJr_JrJr+dJrL_JrL+dJrLz_JrLz))*Ftot
            +(dJr-0.5*(dJrJr_Jr+dJrL_L+dJrLz_Lz))*Ftot_Jr
            -0.5*(dJrJr_Jr*Ftot_Jr+dJrL_Jr*Ftot_L)
            -0.5*(dJrJr*Ftot_JrJr+dJrL*Ftot_JrL))

    # L component

    dL_L = (dL_Lp-dL_Lm)/(2.0*_L0*eps)
    dLL_LL = (dLL_Lp+dLL_Lm-2.0*dLL)/(_L0*eps)^2
    dLLz_LLz = (dLLz_LpLzp+dLLz_LmLzm-dLLz_LpLzm-dLLz_LmLzp)/(2.0*_L0*eps)^2

    dLL_L = (dLL_Lp-dLL_Lm)/(2.0*_L0*eps)
    dLLz_Lz = (dLLz_Lzp-dLLz_Lzm)/(2.0*_L0*eps)

    dLLz_L = (dLLz_Lp-dLLz_Lm)/(2.0*_L0*eps)

    Ftot_LL = (Ftot_Lp+Ftot_Lm-2.0*Ftot)/(_L0*eps)^2


    fluxL_L = ((dL_L-0.5*(dJrL_JrL+dLL_LL+dLLz_LLz))*Ftot
            +(dL-0.5*(dJrL_Jr+dLL_L+dLLz_Lz))*Ftot_L
            -0.5*(dJrL_L*Ftot_Jr+dLL_L*Ftot_L)
            -0.5*(dJrL*Ftot_JrL+dLL*Ftot_LL))

    # Lz component

    dLz_Lz = (dLz_Lzp-dLz_Lzm)/(2.0*_L0*eps)
    dLzLz_LzLz = (dLzLz_Lzp+dLzLz_Lzm-2.0*dLzLz)/(eps*_L0)^2


    fluxLz_Lz = ((dLz_Lz-0.5*(dJrLz_JrLz+dLLz_LLz+dLzLz_LzLz))*Ftot
            -0.5*(dJrLz_Lz*Ftot_Jr+dLLz_Lz*Ftot_L))



    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz), fluxJr_Jr, fluxL_L, fluxLz_Lz
end


function dfluxLdL(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-2), m_test::Float64=m_field)


    if (L > eps)

        Jr_p = Jr + eps*_L0
        Jr_m = Jr - eps*_L0
        L_p = L + eps*_L0
        L_m = L - eps*_L0

        epsLz = min(eps,0.5*(Lz+L_m),0.5*(L_m-Lz))

        Lz_m = Lz - epsLz*_L0
        Lz_p = Lz + epsLz*_L0


        E = E_from_Jr_L(Jr,L,nbu)
        E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
        E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
        E_Lp = E_from_Jr_L(Jr,L_p,nbu)
        E_Lm = E_from_Jr_L(Jr,L_m,nbu)

        E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
        E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
        E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
        E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)



        dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExact(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot= _F(E,L)

        dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


        # Partial derivatives

        dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Jrp = _F(E_Jrp,L)
        dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Jrm = _F(E_Jrm,L)

        dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lp = _F(E_Lp,L_p)
        dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_Lm = _F(E_Lm,L_m)


        # Mixed derivatives

        dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsOptiExact(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_JrpLp = _F(E_Jrp_Lp,L_p)
        dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsOptiExact(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_JrmLm = _F(E_Jrm_Lm,L_m)

        dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsOptiExact(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_JrpLm = _F(E_Jrp_Lm,L_m)
        dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsOptiExact(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Ftot_JrmLp = _F(E_Jrm_Lp,L_p)


        dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)


        # Jr component

        dJrL_JrL = (dJrL_JrpLp+dJrL_JrmLm-dJrL_JrpLm-dJrL_JrmLp)/(2.0*_L0*eps)^2
        dJrL_L = (dJrL_Lp-dJrL_Lm)/(2.0*_L0*eps)
        dJrL_Jr = (dJrL_Jrp-dJrL_Jrm)/(2.0*_L0*eps)

        Ftot_Jr = (Ftot_Jrp-Ftot_Jrm)/(2.0*_L0*eps)
        Ftot_L = (Ftot_Lp-Ftot_Lm)/(2.0*_L0*eps)
        Ftot_JrL = (Ftot_JrpLp+Ftot_JrmLm-Ftot_JrpLm-Ftot_JrmLp)/(2.0*_L0*eps)^2
        Ftot_LL = (Ftot_Lp+Ftot_Lm-2.0*Ftot)/(_L0*eps)^2

        # L component

        dL_L = (dL_Lp-dL_Lm)/(2.0*_L0*eps)
        dLL_LL = (dLL_Lp+dLL_Lm-2.0*dLL)/(_L0*eps)^2
        dLLz_LLz = (dLLz_LpLzp+dLLz_LmLzm-dLLz_LpLzm-dLLz_LmLzp)/(2.0*_L0*eps*2.0*_L0*epsLz)

        dLL_L = (dLL_Lp-dLL_Lm)/(2.0*_L0*eps)

        dLLz_L = (dLLz_Lp-dLLz_Lm)/(2.0*_L0*eps)
        dLLz_Lz = (dLLz_Lzp-dLLz_Lzm)/(2.0*_L0*epsLz)




        fluxL_L = ((dL_L-0.5*(dJrL_JrL+dLL_LL+dLLz_LLz))*Ftot
                +(dL-0.5*(dJrL_Jr+dLL_L+dLLz_Lz))*Ftot_L
                -0.5*(dJrL_L*Ftot_Jr+dLL_L*Ftot_L)
                -0.5*(dJrL*Ftot_JrL+dLL*Ftot_LL))



        return fluxL_L
    else # L < eps ; forward finite diff in L

        Jr_p = Jr + eps*_L0
        Jr_m = Jr - eps*_L0
        L_p = L + eps*_L0
        L_pp = L + 2.0*eps*_L0
        L_ppp = L + 3.0*eps*_L0

        epsLz = min(eps,0.5*(Lz+L),0.5*(L-Lz))

        Lz_m = Lz - epsLz*_L0
        Lz_p = Lz + epsLz*_L0



        E = E_from_Jr_L(Jr,L,nbu)
        E_Jrp_L = E_from_Jr_L(Jr_p,L,nbu)
        E_Jrm_L = E_from_Jr_L(Jr_m,L,nbu)
        E_Jr_Lp = E_from_Jr_L(Jr,L_p,nbu)
        E_Jr_Lpp = E_from_Jr_L(Jr,L_pp,nbu)
        E_Jr_Lppp = E_from_Jr_L(Jr,L_ppp,nbu)

        E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
        E_Jrp_Lpp = E_from_Jr_L(Jr_p,L_pp,nbu)
        E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
        E_Jrm_Lpp = E_from_Jr_L(Jr_m,L_pp,nbu)


        Ftot = _F(E,L)
        Ftot_Lp = _F(E_Jr_Lp,L_p)
        Ftot_Lpp = _F(E_Jr_Lpp,L_pp)
        Ftot_Lppp = _F(E_Jr_Lppp,L_ppp)

        Ftot_Jrp_L = _F(E_Jrp_L,L)
        Ftot_Jrp_Lp = _F(E_Jrp_Lp,L_p)
        Ftot_Jrp_Lpp = _F(E_Jrp_Lpp,L_pp)

        Ftot_Jrm_L = _F(E_Jrm_L,L)
        Ftot_Jrm_Lp = _F(E_Jrm_Lp,L_p)
        Ftot_Jrm_Lpp = _F(E_Jrm_Lpp,L_pp)

        dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExact(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

        dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


        dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Lpp, dL_Lpp, dLz_Lpp, dJrJr_Lpp, dLL_Lpp, dLzLz_Lpp, dJrL_Lpp, dJrLz_Lpp, dLLz_Lpp = orbitAverageActionCoeffsOptiExact(Jr,L_pp,Lz/L_pp,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Lppp, dL_Lppp, dLz_Lppp, dJrJr_Lppp, dLL_Lppp, dLzLz_Lppp, dJrL_Lppp, dJrLz_Lppp, dLLz_Lppp = orbitAverageActionCoeffsOptiExact(Jr,L_pp,Lz/L_pp,m_field,alpha,nbAvr,nbK,nbu,m_test)

        dJr_Jrp_L, dL_Jrp_L, dLz_Jrp_L, dJrJr_Jrp_L, dLL_Jrp_L, dLzLz_Jrp_L, dJrL_Jrp_L, dJrLz_Jrp_L, dLLz_Jrp_L = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Jrp_Lp, dL_Jrp_Lp, dLz_Jrp_Lp, dJrJr_Jrp_Lp, dLL_Jrp_Lp, dLzLz_Jrp_Lp, dJrL_Jrp_Lp, dJrLz_Jrp_Lp, dLLz_Jrp_Lp = orbitAverageActionCoeffsOptiExact(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Jrp_Lpp, dL_Jrp_Lpp, dLz_Jrp_Lpp, dJrJr_Jrp_Lpp, dLL_Jrp_Lpp, dLzLz_Jrp_Lpp, dJrL_Jrp_Lpp, dJrLz_Jrp_Lpp, dLLz_Jrp_Lpp = orbitAverageActionCoeffsOptiExact(Jr_p,L_pp,Lz/L_pp,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Jrm_L, dL_Jrm_L, dLz_Jrm_L, dJrJr_Jrm_L, dLL_Jrm_L, dLzLz_Jrm_L, dJrL_Jrm_L, dJrLz_Jrm_L, dLLz_Jrm_L = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Jrm_Lp, dL_Jrm_Lp, dLz_Jrm_Lp, dJrJr_Jrm_Lp, dLL_Jrm_Lp, dLzLz_Jrm_Lp, dJrL_Jrm_Lp, dJrLz_Jrm_Lp, dLLz_Jrm_Lp = orbitAverageActionCoeffsOptiExact(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Jrm_Lpp, dL_Jrm_Lpp, dLz_Jrm_Lpp, dJrJr_Jrm_Lpp, dLL_Jrm_Lpp, dLzLz_Jrm_Lpp, dJrL_Jrm_Lpp, dJrLz_Jrm_Lpp, dLLz_Jrm_Lpp = orbitAverageActionCoeffsOptiExact(Jr_m,L_pp,Lz/L_pp,m_field,alpha,nbAvr,nbK,nbu,m_test)

        dJr_L_Lzp, dL_L_Lzp, dLz_L_Lzp, dJrJr_L_Lzp, dLL_L_Lzp, dLzLz_L_Lzp, dJrL_L_Lzp, dJrLz_L_Lzp, dLLz_L_Lzp = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Lp_Lzp, dL_Lp_Lzp, dLz_Lp_Lzp, dJrJr_Lp_Lzp, dLL_Lp_Lzp, dLzLz_Lp_Lzp, dJrL_Lp_Lzp, dJrLz_Lp_Lzp, dLLz_Lp_Lzp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Lpp_Lzp, dL_Lpp_Lzp, dLz_Lpp_Lzp, dJrJr_Lpp_Lzp, dLL_Lpp_Lzp, dLzLz_Lpp_Lzp, dJrL_Lpp_Lzp, dJrLz_Lpp_Lzp, dLLz_Lpp_Lzp = orbitAverageActionCoeffsOptiExact(Jr,L_pp,Lz_p/L_pp,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_L_Lzm, dL_L_Lzm, dLz_L_Lzm, dJrJr_L_Lzm, dLL_L_Lzm, dLzLz_L_Lzm, dJrL_L_Lzm, dJrLz_L_Lzm, dLLz_L_Lzm = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Lp_Lzm, dL_Lp_Lzm, dLz_Lp_Lzm, dJrJr_Lp_Lzm, dLL_Lp_Lzm, dLzLz_Lp_Lzm, dJrL_Lp_Lzm, dJrLz_Lp_Lzm, dLLz_Lp_Lzm = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
        dJr_Lpp_Lzm, dL_Lpp_Lzm, dLz_Lpp_Lzm, dJrJr_Lpp_Lzm, dLL_Lpp_Lzm, dLzLz_Lpp_Lzm, dJrL_Lpp_Lzm, dJrLz_Lpp_Lzm, dLLz_Lpp_Lzm = orbitAverageActionCoeffsOptiExact(Jr,L_pp,Lz_m/L_pp,m_field,alpha,nbAvr,nbK,nbu,m_test)


        dL_L = (-1.5*dL+2.0*dL_Lp-0.5*dL_Lpp)/(_L0*eps)
        dJrL_JrL = (-1.5*dJrL_Jrp_L+2.0*dJrL_Jrp_Lp-0.5*dJrL_Jrp_Lpp+1.5*dJrL_Jrm_L-2.0*dJrL_Jrm_Lp+0.5*dJrL_Jrm_Lpp)/(2.0*(_L0*eps)^2)
        dLL_LL = (2.0*dLL-5.0*dLL_Lp+4.0*dLL_Lpp-dLL_Lppp)/(_L0*eps)^2
        dLLz_LLz = (-1.5*dLLz_L_Lzp+2.0*dLLz_Lp_Lzp-0.5*dLLz_Lpp_Lzp+1.5*dLLz_L_Lzm-2.0*dLLz_Lp_Lzm+0.5*dLLz_Lpp_Lzm)/(2.0*(_L0*eps)*(_L0*epsLz))
        dJrL_Jr = (dJrL_Jrp_L - dJrL_Jrm_L)/(2.0*_L0*eps)
        dLL_L = (-1.5*dLL+2.0*dLL_Lp-0.5*dLL_Lpp)/(_L0*eps)
        dLLz_Lz = (dLLz_Lzp-dLLz_Lzm)/(2.0*_L0*epsLz)
        dJrL_L = (-1.5*dJrL+2.0*dJrL_Lp-0.5*dJrL_Lpp)/(_L0*eps)

        Ftot_L = (-1.5*Ftot+2.0*Ftot_Lp-0.5*Ftot_Lpp)/(_L0*eps)
        Ftot_LL = (2.0*Ftot-5.0*Ftot_Lp+4.0*Ftot_Lpp-Ftot_Lppp)/(_L0*eps)^2
        Ftot_Jr = (Ftot_Jrp_L - dJrL_Jrm_L)/(2.0*_L0*eps)
        Ftot_JrL = (-1.5*Ftot_Jrp_L+2.0*Ftot_Jrp_Lp-0.5*Ftot_Jrp_Lpp+1.5*Ftot_Jrm_L-2.0*Ftot_Jrm_Lp+0.5*Ftot_Jrm_Lpp)/(2.0*(_L0*eps)^2)


        fluxL_L = ((dL_L-0.5*(dJrL_JrL+dLL_LL+dLLz_LLz))*Ftot
                +(dL-0.5*(dJrL_Jr+dLL_L+dLLz_Lz))*Ftot_L
                -0.5*(dJrL_L*Ftot_Jr+dLL_L*Ftot_L)
                -0.5*(dJrL*Ftot_JrL+dLL*Ftot_LL))



        return fluxL_L

    end
end


function dFdtOptiExactSignPar(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)



    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot= _F(E,L)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrp = _F(E_Jrp,L)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrm = _F(E_Jrm,L)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lp = _F(E_Lp,L_p)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lm = _F(E_Lm,L_m)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)







    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLp = _F(E_Jrp_Lp,L_p)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLm = _F(E_Jrm_Lm,L_m)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsOptiExactPar(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLm = _F(E_Jrp_Lm,L_m)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsOptiExactPar(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLp = _F(E_Jrm_Lp,L_p)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)



    # Jr component

    dJr_Jr = (dJr_Jrp-dJr_Jrm)/(2.0*_L0*eps)
    dJrJr_JrJr = (dJrJr_Jrp+dJrJr_Jrm-2.0*dJrJr)/(eps*_L0)^2
    dJrL_JrL = (dJrL_JrpLp+dJrL_JrmLm-dJrL_JrpLm-dJrL_JrmLp)/(2.0*_L0*eps)^2
    dJrLz_JrLz = (dJrLz_JrpLzp+dJrLz_JrmLzm-dJrLz_JrpLzm-dJrLz_JrmLzp)/(2.0*_L0*eps)^2

    dJrJr_Jr = (dJrJr_Jrp-dJrJr_Jrm)/(2.0*_L0*eps)
    dJrL_L = (dJrL_Lp-dJrL_Lm)/(2.0*_L0*eps)
    dJrLz_Lz = (dJrLz_Lzp-dJrLz_Lzm)/(2.0*_L0*eps)
    Ftot_Jr = (Ftot_Jrp-Ftot_Jrm)/(2.0*_L0*eps)

    dJrL_Jr = (dJrL_Jrp-dJrL_Jrm)/(2.0*_L0*eps)
    dJrLz_Jr = (dJrLz_Jrp-dJrLz_Jrm)/(2.0*_L0*eps)
    Ftot_L = (Ftot_Lp-Ftot_Lm)/(2.0*_L0*eps)

    Ftot_JrJr = (Ftot_Jrp+Ftot_Jrm-2.0*Ftot)/(_L0*eps)^2
    Ftot_JrL = (Ftot_JrpLp+Ftot_JrmLm-Ftot_JrpLm-Ftot_JrmLp)/(2.0*_L0*eps)^2

    fluxJr_Jr = ((dJr_Jr-0.5*(dJrJr_JrJr+dJrL_JrL+dJrLz_JrLz))*Ftot
            +(dJr-0.5*(dJrJr_Jr+dJrL_L+dJrLz_Lz))*Ftot_Jr
            -0.5*(dJrJr_Jr*Ftot_Jr+dJrL_Jr*Ftot_L)
            -0.5*(dJrJr*Ftot_JrJr+dJrL*Ftot_JrL))

    # L component

    dL_L = (dL_Lp-dL_Lm)/(2.0*_L0*eps)
    dLL_LL = (dLL_Lp+dLL_Lm-2.0*dLL)/(_L0*eps)^2
    dLLz_LLz = (dLLz_LpLzp+dLLz_LmLzm-dLLz_LpLzm-dLLz_LmLzp)/(2.0*_L0*eps)^2

    dLL_L = (dLL_Lp-dLL_Lm)/(2.0*_L0*eps)
    dLLz_Lz = (dLLz_Lzp-dLLz_Lzm)/(2.0*_L0*eps)

    dLLz_L = (dLLz_Lp-dLLz_Lm)/(2.0*_L0*eps)

    Ftot_LL = (Ftot_Lp+Ftot_Lm-2.0*Ftot)/(_L0*eps)^2


    fluxL_L = ((dL_L-0.5*(dJrL_JrL+dLL_LL+dLLz_LLz))*Ftot
            +(dL-0.5*(dJrL_Jr+dLL_L+dLLz_Lz))*Ftot_L
            -0.5*(dJrL_L*Ftot_Jr+dLL_L*Ftot_L)
            -0.5*(dJrL*Ftot_JrL+dLL*Ftot_LL))

    # Lz component

    dLz_Lz = (dLz_Lzp-dLz_Lzm)/(2.0*_L0*eps)
    dLzLz_LzLz = (dLzLz_Lzp+dLzLz_Lzm-2.0*dLzLz)/(eps*_L0)^2


    fluxLz_Lz = ((dLz_Lz-0.5*(dJrLz_JrLz+dLLz_LLz+dLzLz_LzLz))*Ftot
            -0.5*(dJrLz_Lz*Ftot_Jr+dLLz_Lz*Ftot_L))


    # println("Lz     = ",Lz)
    # println("FJr_Jr = ",fluxJr_Jr)
    # println("FL_L   = ",fluxL_L)
    # println("DLzLz = ",dLzLz)
    # println("-----------------------")

    #dLzLz_Lz = (dLzLz_Lzp-dLzLz_Lzm)/(2.0*eps*_L0)

    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz)#, fluxJr_Jr,fluxL_L,fluxLz_Lz, dLz_Lz, dJrLz_JrLz, dLLz_LLz, dLzLz_LzLz, dJrLz_Lz, dLLz_Lz, dLzLz_Lz
end

function dFdtOptiExactSignParDecomp(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)



    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot= _F(E,L)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrp = _F(E_Jrp,L)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrm = _F(E_Jrm,L)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lp = _F(E_Lp,L_p)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lm = _F(E_Lm,L_m)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)







    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLp = _F(E_Jrp_Lp,L_p)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLm = _F(E_Jrm_Lm,L_m)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsOptiExactPar(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLm = _F(E_Jrp_Lm,L_m)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsOptiExactPar(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLp = _F(E_Jrm_Lp,L_p)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)



    # Jr component

    dJr_Jr = (dJr_Jrp-dJr_Jrm)/(2.0*_L0*eps)
    dJrJr_JrJr = (dJrJr_Jrp+dJrJr_Jrm-2.0*dJrJr)/(eps*_L0)^2
    dJrL_JrL = (dJrL_JrpLp+dJrL_JrmLm-dJrL_JrpLm-dJrL_JrmLp)/(2.0*_L0*eps)^2
    dJrLz_JrLz = (dJrLz_JrpLzp+dJrLz_JrmLzm-dJrLz_JrpLzm-dJrLz_JrmLzp)/(2.0*_L0*eps)^2

    dJrJr_Jr = (dJrJr_Jrp-dJrJr_Jrm)/(2.0*_L0*eps)
    dJrL_L = (dJrL_Lp-dJrL_Lm)/(2.0*_L0*eps)
    dJrLz_Lz = (dJrLz_Lzp-dJrLz_Lzm)/(2.0*_L0*eps)
    Ftot_Jr = (Ftot_Jrp-Ftot_Jrm)/(2.0*_L0*eps)

    dJrL_Jr = (dJrL_Jrp-dJrL_Jrm)/(2.0*_L0*eps)
    dJrLz_Jr = (dJrLz_Jrp-dJrLz_Jrm)/(2.0*_L0*eps)
    Ftot_L = (Ftot_Lp-Ftot_Lm)/(2.0*_L0*eps)

    Ftot_JrJr = (Ftot_Jrp+Ftot_Jrm-2.0*Ftot)/(_L0*eps)^2
    Ftot_JrL = (Ftot_JrpLp+Ftot_JrmLm-Ftot_JrpLm-Ftot_JrmLp)/(2.0*_L0*eps)^2

    fluxJr_Jr = ((dJr_Jr-0.5*(dJrJr_JrJr+dJrL_JrL+dJrLz_JrLz))*Ftot
            +(dJr-0.5*(dJrJr_Jr+dJrL_L+dJrLz_Lz))*Ftot_Jr
            -0.5*(dJrJr_Jr*Ftot_Jr+dJrL_Jr*Ftot_L)
            -0.5*(dJrJr*Ftot_JrJr+dJrL*Ftot_JrL))

    # L component

    dL_L = (dL_Lp-dL_Lm)/(2.0*_L0*eps)
    dLL_LL = (dLL_Lp+dLL_Lm-2.0*dLL)/(_L0*eps)^2
    dLLz_LLz = (dLLz_LpLzp+dLLz_LmLzm-dLLz_LpLzm-dLLz_LmLzp)/(2.0*_L0*eps)^2

    dLL_L = (dLL_Lp-dLL_Lm)/(2.0*_L0*eps)
    dLLz_Lz = (dLLz_Lzp-dLLz_Lzm)/(2.0*_L0*eps)

    dLLz_L = (dLLz_Lp-dLLz_Lm)/(2.0*_L0*eps)

    Ftot_LL = (Ftot_Lp+Ftot_Lm-2.0*Ftot)/(_L0*eps)^2


    fluxL_L = ((dL_L-0.5*(dJrL_JrL+dLL_LL+dLLz_LLz))*Ftot
            +(dL-0.5*(dJrL_Jr+dLL_L+dLLz_Lz))*Ftot_L
            -0.5*(dJrL_L*Ftot_Jr+dLL_L*Ftot_L)
            -0.5*(dJrL*Ftot_JrL+dLL*Ftot_LL))

    # Lz component

    dLz_Lz = (dLz_Lzp-dLz_Lzm)/(2.0*_L0*eps)
    dLzLz_LzLz = (dLzLz_Lzp+dLzLz_Lzm-2.0*dLzLz)/(eps*_L0)^2


    fluxLz_Lz = ((dLz_Lz-0.5*(dJrLz_JrLz+dLLz_LLz+dLzLz_LzLz))*Ftot
            -0.5*(dJrLz_Lz*Ftot_Jr+dLLz_Lz*Ftot_L))


    # println("Lz     = ",Lz)
    # println("FJr_Jr = ",fluxJr_Jr)
    # println("FL_L   = ",fluxL_L)
    # println("DLzLz = ",dLzLz)
    # println("-----------------------")

    #dLzLz_Lz = (dLzLz_Lzp-dLzLz_Lzm)/(2.0*eps*_L0)

    return fluxJr_Jr,fluxL_L,fluxLz_Lz#, dLz_Lz, dJrLz_JrLz, dLLz_LLz, dLzLz_LzLz, dJrLz_Lz, dLLz_Lz, dLzLz_Lz
end




function FluxLzOptiExactSignParDecomp(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)

    E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
    E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
    E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
    E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)



    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot= _F(E,L)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrp = _F(E_Jrp,L)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrm = _F(E_Jrm,L)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lp = _F(E_Lp,L_p)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lm = _F(E_Lm,L_m)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)







    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLp = _F(E_Jrp_Lp,L_p)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLm = _F(E_Jrm_Lm,L_m)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsOptiExactPar(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLm = _F(E_Jrp_Lm,L_m)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsOptiExactPar(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLp = _F(E_Jrm_Lp,L_p)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)


    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)



    # Jr component

    dJr_Jr = (dJr_Jrp-dJr_Jrm)/(2.0*_L0*eps)
    dJrJr_JrJr = (dJrJr_Jrp+dJrJr_Jrm-2.0*dJrJr)/(eps*_L0)^2
    dJrL_JrL = (dJrL_JrpLp+dJrL_JrmLm-dJrL_JrpLm-dJrL_JrmLp)/(2.0*_L0*eps)^2
    dJrLz_JrLz = (dJrLz_JrpLzp+dJrLz_JrmLzm-dJrLz_JrpLzm-dJrLz_JrmLzp)/(2.0*_L0*eps)^2

    dJrJr_Jr = (dJrJr_Jrp-dJrJr_Jrm)/(2.0*_L0*eps)
    dJrL_L = (dJrL_Lp-dJrL_Lm)/(2.0*_L0*eps)
    dJrLz_Lz = (dJrLz_Lzp-dJrLz_Lzm)/(2.0*_L0*eps)
    Ftot_Jr = (Ftot_Jrp-Ftot_Jrm)/(2.0*_L0*eps)

    dJrL_Jr = (dJrL_Jrp-dJrL_Jrm)/(2.0*_L0*eps)
    dJrLz_Jr = (dJrLz_Jrp-dJrLz_Jrm)/(2.0*_L0*eps)
    Ftot_L = (Ftot_Lp-Ftot_Lm)/(2.0*_L0*eps)

    Ftot_JrJr = (Ftot_Jrp+Ftot_Jrm-2.0*Ftot)/(_L0*eps)^2
    Ftot_JrL = (Ftot_JrpLp+Ftot_JrmLm-Ftot_JrpLm-Ftot_JrmLp)/(2.0*_L0*eps)^2

    fluxJr_Jr = ((dJr_Jr-0.5*(dJrJr_JrJr+dJrL_JrL+dJrLz_JrLz))*Ftot
            +(dJr-0.5*(dJrJr_Jr+dJrL_L+dJrLz_Lz))*Ftot_Jr
            -0.5*(dJrJr_Jr*Ftot_Jr+dJrL_Jr*Ftot_L)
            -0.5*(dJrJr*Ftot_JrJr+dJrL*Ftot_JrL))

    # L component

    dL_L = (dL_Lp-dL_Lm)/(2.0*_L0*eps)
    dLL_LL = (dLL_Lp+dLL_Lm-2.0*dLL)/(_L0*eps)^2
    dLLz_LLz = (dLLz_LpLzp+dLLz_LmLzm-dLLz_LpLzm-dLLz_LmLzp)/(2.0*_L0*eps)^2

    dLL_L = (dLL_Lp-dLL_Lm)/(2.0*_L0*eps)
    dLLz_Lz = (dLLz_Lzp-dLLz_Lzm)/(2.0*_L0*eps)

    dLLz_L = (dLLz_Lp-dLLz_Lm)/(2.0*_L0*eps)

    Ftot_LL = (Ftot_Lp+Ftot_Lm-2.0*Ftot)/(_L0*eps)^2


    fluxL_L = ((dL_L-0.5*(dJrL_JrL+dLL_LL+dLLz_LLz))*Ftot
            +(dL-0.5*(dJrL_Jr+dLL_L+dLLz_Lz))*Ftot_L
            -0.5*(dJrL_L*Ftot_Jr+dLL_L*Ftot_L)
            -0.5*(dJrL*Ftot_JrL+dLL*Ftot_LL))

    # Lz component

    dLz_Lz = (dLz_Lzp-dLz_Lzm)/(2.0*_L0*eps)
    dLzLz_LzLz = (dLzLz_Lzp+dLzLz_Lzm-2.0*dLzLz)/(eps*_L0)^2


    fluxLz_Lz = ((dLz_Lz-0.5*(dJrLz_JrLz+dLLz_LLz+dLzLz_LzLz))*Ftot
            -0.5*(dJrLz_Lz*Ftot_Jr+dLLz_Lz*Ftot_L))


    # println("Lz     = ",Lz)
    # println("FJr_Jr = ",fluxJr_Jr)
    # println("FL_L   = ",fluxL_L)
    # println("DLzLz = ",dLzLz)
    # println("-----------------------")

    #dLzLz_Lz = (dLzLz_Lzp-dLzLz_Lzm)/(2.0*eps*_L0)

    return dLz_Lz*Ftot, dJrLz_JrLz*Ftot,dLLz_LLz*Ftot,dLzLz_LzLz*Ftot,dJrLz_Lz*Ftot_Jr,dLLz_Lz*Ftot_L#, dLz_Lz, dJrLz_JrLz, dLLz_LLz, dLzLz_LzLz, dJrLz_Lz, dLLz_Lz, dLzLz_Lz
end








##################################################################
# Testing flux components
##################################################################




function FluxComponentTestJr(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    # Jr component


end













##################################################################
##################################################################

# Intergrate div(F) over Lz
function integrate_dFdt_Lz(Jr::Float64, L::Float64, m_field::Float64, nbLz::Int64=12,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)


    sum = 0.0

    for iLz=1:nbLz
        Lz = -L + 2.0*L*(iLz-0.5)/nbLz
        dfdt = dFdt(Jr,L,Lz,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,eps,m_test)

        sum += dfdt
    end

    sum *= 2.0*L/nbLz

    return sum

end

function integrate_dFdtOpti_Lz(Jr::Float64, L::Float64, m_field::Float64, nbLz::Int64=12,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)


    sum = 0.0

    for iLz=1:nbLz
        Lz = -L + 2.0*L*(iLz-0.5)/nbLz
        dfdt = dFdtOpti(Jr,L,Lz,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,eps,m_test)

        sum += dfdt
    end

    sum *= 2.0*L/nbLz

    return sum

end

function integrate_dFdtOpti_Lz(Jr::Float64, L::Float64, m_field::Float64, nbLz::Int64=12,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)


    sum = 0.0

    for iLz=1:nbLz
        Lz = -L + 2.0*L*(iLz-0.5)/nbLz
        dfdt = dFdtOpti(Jr,L,Lz,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,eps,m_test)

        sum += dfdt
    end

    sum *= 2.0*L/nbLz

    return sum

end

function integrate_dFdt_Lz_Par(Jr::Float64, L::Float64, m_field::Float64, nbLz::Int64=12,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)


    sum = 0.0

    for iLz=1:nbLz
        Lz = -L + 2.0*L*(iLz-0.5)/nbLz
        dfdt = dFdt_Par(Jr,L,Lz,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,eps,m_test)

        sum += dfdt
    end

    sum *= 2.0*L/nbLz

    return sum

end

function integrate_dFdtOpti_Lz_Par(Jr::Float64, L::Float64, m_field::Float64, nbLz::Int64=12,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)


    sum = 0.0

    for iLz=1:nbLz
        Lz = -L + 2.0*L*(iLz-0.5)/nbLz
        dfdt = dFdtOpti_Par(Jr,L,Lz,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,eps,m_test)

        sum += dfdt
    end

    sum *= 2.0*L/nbLz

    return sum

end

##################################################################
# dF/dt 3D for g(x)=sign(x)
# Lz-integration
# Explicitely takes into account integration of dirac(Lz)
##################################################################

function integrate_dFdtOptiExactSign_Lz(Jr::Float64, L::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbLz::Int64=12, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5))


    sumReg = 0.0

    # Integrate regular part

    for iLz=1:nbLz
        Lz = -L + 2.0*L*(iLz-0.5)/nbLz

        ##### Compute dfdt, regular part
        epsGrad = min(0.5*Jr,0.5*L,0.5*(Lz+L),0.5*(L-Lz),eps)
        dfdt, _ = dFdtOptiExactSign(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,epsGrad,m_test)

        sumReg += dfdt
    end

    sumReg *= 2.0*L/nbLz

    # println("sumReg = ",sumReg)


    # Integrate the singular (Dirac) part
    # Evaluation is at Lz=0

    Lz = 0.0

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)


    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExact(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot= _F(E,L)



    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExact(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrp = _F(E_Jrp,L)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExact(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrm = _F(E_Jrm,L)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExact(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lp = _F(E_Lp,L_p)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExact(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lm = _F(E_Lm,L_m)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExact(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)




    # Jr component

    dJrLz_Jr = (dJrLz_Jrp-dJrLz_Jrm)/(2.0*eps*_L0)
    Ftot_Jr = (Ftot_Jrp-Ftot_Jrm)/(2.0*eps*_L0)

    sumSingJr = -0.5*alpha*Ftot*dJrLz_Jr-0.5*alpha*dJrLz*Ftot_Jr


    # L component

    dLLz_L = (dLLz_Lp-dLLz_Lm)/(2.0*eps*_L0)
    Ftot_L = (Ftot_Lp-Ftot_Lm)/(2.0*eps*_L0)

    sumSingL = -0.5*alpha*Ftot*dLLz_L-0.5*alpha*dLLz*Ftot_L

    # Lz component

    dJrLz_Jr = (dJrLz_Jrp-dJrLz_Jrm)/(2.0*eps*_L0)
    dLLz_L = (dLLz_Lp-dLLz_Lm)/(2.0*eps*_L0)
    dLzLz_Lz = (dLzLz_Lzp-dLzLz_Lzm)/(2.0*eps*_L0)

    sumSingLz = (alpha*Ftot*(dLz-0.5*(dJrLz_Jr+dLLz_L+dLzLz_Lz))
                -0.5*alpha*Ftot*dLzLz_Lz
                -0.5*alpha*(Ftot_Jr*dJrLz+Ftot_L*dLLz-Ftot*dLzLz_Lz))



    return sumReg-sumSingJr-sumSingL-sumSingLz

end


function integrate_dFdtOptiExactSign_Lz_Par(Jr::Float64, L::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field,
            nbK::Int64=nbK_default, nbu::Int64=nbu0, eps::Float64=10^(-5), nbLz::Int64=12)


    sumReg = 0.0

    # Integrate regular part

    for iLz=1:nbLz
        Lz = -L + 2.0*L*(iLz-0.5)/nbLz

        ##### Compute dfdt, regular part
        dfdt = dFdtOptiExactSignPar(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,eps,m_test)

        sumReg += dfdt
    end

    sumReg *= 2.0*L/nbLz

    # println("sumReg = ",sumReg)


    # Integrate the singular (Dirac) part
    # Evaluation is at Lz=0

    Lz = 0.0

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    E = E_from_Jr_L(Jr,L,nbu)
    E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
    E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
    E_Lp = E_from_Jr_L(Jr,L_p,nbu)
    E_Lm = E_from_Jr_L(Jr,L_m,nbu)


    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot= _F(E,L)



    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsOptiExactPar(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrp = _F(E_Jrp,L)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsOptiExactPar(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Jrm = _F(E_Jrm,L)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsOptiExactPar(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lp = _F(E_Lp,L_p)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsOptiExactPar(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
    Ftot_Lm = _F(E_Lm,L_m)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbK,nbu,m_test)




    # Jr component

    dJrLz_Jr = (dJrLz_Jrp-dJrLz_Jrm)/(2.0*eps*_L0)
    Ftot_Jr = (Ftot_Jrp-Ftot_Jrm)/(2.0*eps*_L0)

    sumSingJr = -0.5*alpha*Ftot*dJrLz_Jr-0.5*alpha*dJrLz*Ftot_Jr


    # L component

    dLLz_L = (dLLz_Lp-dLLz_Lm)/(2.0*eps*_L0)
    Ftot_L = (Ftot_Lp-Ftot_Lm)/(2.0*eps*_L0)

    sumSingL = -0.5*alpha*Ftot*dLLz_L-0.5*alpha*dLLz*Ftot_L

    # Lz component

    dJrLz_Jr = (dJrLz_Jrp-dJrLz_Jrm)/(2.0*eps*_L0)
    dLLz_L = (dLLz_Lp-dLLz_Lm)/(2.0*eps*_L0)
    dLzLz_Lz = (dLzLz_Lzp-dLzLz_Lzm)/(2.0*eps*_L0)

    sumSingLz = (alpha*Ftot*(dLz-0.5*(dJrLz_Jr+dLLz_L+dLzLz_Lz))
                -0.5*alpha*Ftot*dLzLz_Lz
                -0.5*alpha*(Ftot_Jr*dJrLz+Ftot_L*dLLz-Ftot*dLzLz_Lz))



    return sumReg-sumSingJr-sumSingL-sumSingLz

end

# ##################################################################
# # dFdt_2D : integrate over cos I
# ##################################################################
#
# function dFdt_Rot2D(Jr::Float64, L::Float64, m_field::Float64, alpha::Float64=alphaRot,
#             nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
#             m_test::Float64=m_field, eps::Float64=10^(-5))
#
#     Jr_p = Jr + eps*_L0
#     Jr_m = Jr - eps*_L0
#     L_p = L + eps*_L0
#     L_m = L - eps*_L0
#
#     E = E_from_Jr_L(Jr,L,nbu)
#     E_Jrp = E_from_Jr_L(Jr_p,L,nbu)
#     E_Jrm = E_from_Jr_L(Jr_m,L,nbu)
#     E_Lp = E_from_Jr_L(Jr,L_p,nbu)
#     E_Lm = E_from_Jr_L(Jr,L_m,nbu)
#
#
#     dJr_1, dL_1, dI_1, dJrJr_1, dLL_1, dII_1, dJrL_1, dJrI_1, dLI_1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     Fred = 2.0*L*_F(E,L)
#
#     # Jr
#     dJr_Jrp_1, dL_Jrp_1, dI_Jrp_1, dJrJr_Jrp_1, dLL_Jrp_1, dII_Jrp_1, dJrL_Jrp_1, dJrI_Jrp_1, dLI_Jrp_1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr_p,L,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     dJr_Jrm_1, dL_Jrm_1, dI_Jrm_1, dJrJr_Jrm_1, dLL_Jrm_1, dII_Jrm_1, dJrL_Jrm_1, dJrI_Jrm_1, dLI_Jrm_1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr_m,L,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     Fred_Jrp = 2.0*L*_F(E_Jrp,L)
#     Fred_Jrm = 2.0*L*_F(E_Jrm,L)
#
#     dJr = (dJr_Jrp_1*Fred_Jrp-dJr_Jrm_1*Fred_Jrm)/(2.0*_L0*eps)
#
#     # L
#     dJr_Lp_1, dL_Lp_1, dI_Lp_1, dJrJr_Lp_1, dLL_Lp_1, dII_Lp_1, dJrL_Lp_1, dJrI_Lp_1, dLI_Lp_1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     dJr_Lm_1, dL_Lm_1, dI_Lm_1, dJrJr_Lm_1, dLL_Lm_1, dII_Lm_1, dJrL_Lm_1, dJrI_Lm_1, dLI_Lm_1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_m,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     Fred_Lp = 2.0*L_p*_F(E_Lp,L_p)
#     Fred_Lm = 2.0*L_m*_F(E_Lm,L_m)
#
#     dL = (dL_Lp_1*Fred_Lp-dL_Lm_1*Fred_Lm)/(2.0*_L0*eps)
#
#     # JrJr
#     dJrJr = (dJrJr_Jrp_1*Fred_Jrp + dJrJr_Jrm_1*Fred_Jrm - 2.0*dJrJr_1*Fred)/(_L0*eps)^2
#
#     # LL
#     dLL = (dLL_Lp_1*Fred_Lp + dLL_Lm_1*Fred_Lm - 2.0*dLL_1*Fred)/(_L0*eps)^2
#
#     # JrL
#     E_Jrp_Lp = E_from_Jr_L(Jr_p,L_p,nbu)
#     E_Jrp_Lm = E_from_Jr_L(Jr_p,L_m,nbu)
#     E_Jrm_Lp = E_from_Jr_L(Jr_m,L_p,nbu)
#     E_Jrm_Lm = E_from_Jr_L(Jr_m,L_m,nbu)
#
#     _, _, _, _, _, _, dJrL_JrpLp, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_p,L_p,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     _, _, _, _, _, _, dJrL_JrmLm, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_m,L_m,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     _, _, _, _, _, _, dJrL_JrpLm, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_p,L_m,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     _, _, _, _, _, _, dJrL_JrmLp, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr_m,L_p,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#
#     Fred_JrpLp = 2.0*L_p*_F(E_Jrp_Lp,L_p)
#     Fred_JrmLm = 2.0*L_m*_F(E_Jrm_Lm,L_m)
#     Fred_JrpLm = 2.0*L_m*_F(E_Jrp_Lm,L_m)
#     Fred_JrmLp = 2.0*L_p*_F(E_Jrm_Lp,L_p)
#
#     dJrL = (dJrL_JrpLp*Fred_JrpLp + dJrL_JrmLm*Fred_JrmLm - dJrL_JrpLm*Fred_JrpLm - dJrL_JrmLp*Fred_JrmLp)/(2.0*_L0*eps)^2
#
#
#     # I
#     dJr_m1, dL_m1, dI_m1, dJrJr_m1, dLL_m1, dII_m1, dJrL_m1, dJrI_m1, dLI_m1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     Frot_1 = _Frot(E,L,L,alpha)
#     Frot_m1 = _Frot(E,L,-L,alpha)
#
#     dI = dI_1*Frot_1 - dI_m1*Frot_m1
#
#     # JrI
#     dJr_Jrp_m1, dL_Jrp_m1, dI_Jrp_m1, dJrJr_Jrp_m1, dLL_Jrp_m1, dII_Jrp_m1, dJrL_Jrp_m1, dJrI_Jrp_m1, dLI_Jrp_m1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr_p,L,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     dJr_Jrm_m1, dL_Jrm_m1, dI_Jrm_m1, dJrJr_Jrm_m1, dLL_Jrm_m1, dII_Jrm_m1, dJrL_Jrm_m1, dJrI_Jrm_m1, dLI_Jrm_m1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr_m,L,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#
#     Frot_Jrp_1 = _Frot(E_Jrp,L,L,alpha)
#     Frot_Jrm_1 = _Frot(E_Jrm,L,L,alpha)
#     Frot_Jrp_m1 = _Frot(E_Jrp,L,-L,alpha)
#     Frot_Jrm_m1 = _Frot(E_Jrm,L,-L,alpha)
#
#     dJrI = (dJrI_Jrp_1*Frot_Jrp_1-dJrI_Jrm_1*Frot_Jrm_1)/(2.0*_L0*eps) - (dJrI_Jrp_m1*Frot_Jrp_m1-dJrI_Jrm_m1*Frot_Jrm_m1)/(2.0*_L0*eps)
#
#     # LI
#     dJr_Lp_m1, dL_Lp_m1, dI_Lp_m1, dJrJr_Lp_m1, dLL_Lp_m1, dII_Lp_m1, dJrL_Lp_m1, dJrI_Lp_m1, dLI_Lp_m1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_p,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     dJr_Lm_m1, dL_Lm_m1, dI_Lm_m1, dJrJr_Lm_m1, dLL_Lm_m1, dII_Lm_m1, dJrL_Lm_m1, dJrI_Lm_m1, dLI_Lm_m1 = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L_m,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#
#     Frot_Lp_1 = _Frot(E_Lp,L_p,L_p,alpha)
#     Frot_Lm_1 = _Frot(E_Lm,L_m,L_m,alpha)
#     Frot_Lp_m1 = _Frot(E_Lp,L_p,-L_p,alpha)
#     Frot_Lm_m1 = _Frot(E_Lm,L_m,-L_m,alpha)
#
#     dLI = (dLI_Lp_1*Frot_Lp_1-dLI_Lm_1*Frot_Lm_1)/(2.0*_L0*eps) - (dLI_Lp_m1*Frot_Lp_m1-dLI_Lm_m1*Frot_Lm_m1)/(2.0*_L0*eps)
#
#     # II
#     cosI_1_m = 1.0 - eps
#     cosI_1_mm = 1.0 - 2.0*eps
#
#     cosI_m1_p = -1.0 + eps
#     cosI_m1_pp = -1.0 + 2.0*eps
#
#     _, _, _, _, _, dII_1, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     _, _, _, _, _, dII_m1, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,-1.0,m_field,alpha,nbAvr,nbK,nbu,m_test)
#
#     _, _, _, _, _, dII_1_m, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_1_m,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     _, _, _, _, _, dII_1_mm, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_1_mm,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     _, _, _, _, _, dII_m1_p, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_m1_p,m_field,alpha,nbAvr,nbK,nbu,m_test)
#     _, _, _, _, _, dII_m1_pp, _, _, _ = orbitAverageActionCoeffs_cosI_OptiExact(Jr,L,cosI_m1_pp,m_field,alpha,nbAvr,nbK,nbu,m_test)
#
#     Frot_1 = _Frot(E,L,L,alpha)
#     Frot_1_m = _Frot(E,L,L*cosI_1_m,alpha)
#     Frot_1_mm = _Frot(E,L,L*cosI_1_mm,alpha)
#     Frot_m1 = _Frot(E,L,-L,alpha)
#     Frot_m1_p = _Frot(E,L,L*cosI_m1_p,alpha)
#     Frot_m1_pp = _Frot(E,L,L*cosI_m1_pp,alpha)
#
#     # f'(1) = [1.5*f(1)-2.0*f(1-eps)+0.5*f(1-2*eps)/eps
#     dcosI_1 = (1.5*dII_1*Frot_1 - 2.0*dII_1_m*Frot_1_m + 0.5*dII_1_mm*Frot_1_mm)/(eps)
#
#     # f'(-1) = [-1.5*f(-1)+2.0*f(-1+eps)-0.5*f(-1+2*eps)/eps
#     dcosI_m1 = (-1.5*dII_m1*Frot_m1 + 2.0*dII_m1_p*Frot_m1_p - 0.5*dII_m1_pp*Frot_m1_pp)/(eps)
#
#     dII = dcosI_1-dcosI_m1
#
#     println("dJr:",dJr-0.5*dJrJr-0.5*dJrL)
#     println("dL:",dL-0.5*dLL-0.5*dJrL)
#
#     #leftover = Leftover(Jr,L,m_field,alpha,nbAvr,nbK,nbu,eps,m_test)
#
#
#     return -dJr-dL+0.5*dJrJr+0.5*dLL+dJrL - (dI - dJrI - dLI - 0.5*dII)
# end

##################################################################
# Compute dFdt(cosI)
# F(cosI) = (2pi)^3 int dJr dL L F(Jr,L,LcosI)
##################################################################

function integrateJrLPow_dFdtPar(cosI::Float64, m_field::Float64,
            alpha::Float64=alphaRot, pow::Int64=3,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbAvrTh::Int64=nbAvrTh_default=20, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-2),
            Jrmax::Float64=3.0*_L0, nbJr::Int64=30,
            Lmax::Float64=2.0*_L0, nbL::Int64=20)

    sum = Threads.Atomic{Float64}(0.0)
    xmax = (Jrmax/_L0)^(1/pow)
    ymax = (Lmax/_L0)^(1/pow)

    nbJrL = nbJr*nbL

    tab_xy = zeros(Float64,2,nbJrL)

    iGrid = 1
    for ix=1:nbJr, iy=1:nbL
        x = xmax * (ix-0.5)/nbJr
        y = ymax * (iy-0.5)/nbL
        tab_xy[1,iGrid], tab_xy[2,iGrid] = x, y
        iGrid += 1
    end



    Threads.@threads for iGrid=1:nbJrL
        x = tab_xy[1,iGrid]
        Jr = _L0 * x^pow
        y = tab_xy[2,iGrid]
        L = _L0 * y^pow
        Lz = L*cosI

        # This is so that finite differences are well-defined, i.e. L,Jr>0 and -L<Lz<L as arguments
        epsEff = min(eps,0.01*Jr/_L0,0.01*L/_L0,0.01*abs(L-Lz)/_L0,0.01*abs(-L-Lz)/_L0)


        dfdt = dFdt(Jr,L,L*cosI,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,epsEff,m_test)


        Threads.atomic_add!(sum,x^(pow-1)*y^(pow-1)*L*dfdt)

    end

    sum[] *=  (2*pi)^3*xmax/nbJr * pow*_L0 * ymax/nbL * pow*_L0

    return sum[]

end

function integrateJrLOptiExactSignPow_dFdtPar(cosI::Float64, m_field::Float64,
            alpha::Float64=alphaRot, pow::Int64=3,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-2),
            Jrmax::Float64=3.0*_L0, nbJr::Int64=30,
            Lmax::Float64=2.0*_L0, nbL::Int64=20,
            Lmin::Float64=2.0*10^(-2))

    sum = Threads.Atomic{Float64}(0.0)
    xmax = (Jrmax/_L0)^(1/pow)

    ymax = (Lmax/_L0)^(1/pow)
    ymin = (Lmin/_L0)^(1/pow)

    nbJrL = nbJr*nbL

    tab_xy = zeros(Float64,2,nbJrL)

    iGrid = 1
    for ix=1:nbJr, iy=1:nbL
        x = xmax * (ix-0.5)/nbJr
        y = ymin + (ymax-ymin) * (iy-0.5)/nbL
        tab_xy[1,iGrid], tab_xy[2,iGrid] = x, y
        iGrid += 1
    end



    Threads.@threads for iGrid=1:nbJrL
        x = tab_xy[1,iGrid]
        Jr = _L0 * x^pow
        y = tab_xy[2,iGrid]
        L = _L0 * y^pow
        Lz = L*cosI

        # This is so that finite differences are well-defined, i.e. L,Jr>0 and -L<Lz<L as arguments
        epsEff = min(eps,0.01*Jr/_L0,0.01*L/_L0,0.01*abs(L-Lz)/_L0,0.01*abs(-L-Lz)/_L0)


        dfdt, _ = dFdtOptiExactSign(Jr,L,L*cosI,m_field,alpha,nbAvr,nbK,nbu,epsEff,m_test)


        Threads.atomic_add!(sum,x^(pow-1)*y^(pow-1)*L*dfdt)

    end

    sum[] *=  (2*pi)^3*xmax/nbJr * pow*_L0 * (ymax-ymin)/nbL * pow*_L0

    return sum[], Lmin

end
