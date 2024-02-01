
# Computes the 3D flux in (Jr,L,Lz) space
function flux3D(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default, nbvartheta::Int64=nbvartheta_default, nbphi::Int64=nbphi_default, 
            nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

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

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffs(Jr,L,Lz/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)


    Frot = _Frot(E,L,Lz,alpha)


    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffs(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffs(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffs(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffs(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffs(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffs(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
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


# Computes the 2D flux in (Jr,L) space
function flux2D_JrL(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default,
            nbvartheta::Int64=nbvartheta_default, nbphi::Int64=nbphi_default,
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

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffs(Jr,L,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)


    Frot = _Frot_cosI(E,L,cosI,alpha)


    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffs(Jr_p,L,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Jrp = _Frot_cosI(E_Jrp,L,cosI,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffs(Jr_m,L,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Jrm = _Frot_cosI(E_Jrm,L,cosI,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffs(Jr,L_p,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Lp = _Frot_cosI(E_Lp,L_p,cosI,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffs(Jr,L_m,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
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

# Computes the 2D-diffusion rate dF/dt in (Jr,L) space
function dFdt2D_JrL(Jr::Float64, L::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbCosI::Int64=50, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default,
            nbvartheta::Int64=nbvartheta_default, nbphi::Int64=nbphi_default,
            nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    sumJr = 0.0
    sumL = 0.0

    for i=1:nbCosI
        cosI = -1.0 + 2.0/nbCosI*(i-0.5)

        fJr_p, _ = flux2D_JrL(Jr+eps,L,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,eps,m_test)
        fJr_m, _ = flux2D_JrL(Jr-eps,L,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,eps,m_test)

        dJr = (fJr_p-fJr_m)/(2.0*eps)

        _, fL_p  = flux2D_JrL(Jr,L+eps,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,eps,m_test)
        _, fL_m  = flux2D_JrL(Jr,L-eps,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,eps,m_test)

        dL = (fL_p-fL_m)/(2.0*eps)

        sumJr += dJr
        sumL  += dL
    end

    sumJr *= 2.0/nbCosI
    sumL *= 2.0/nbCosI

    return -(sumJr+sumL)
end

# Computes the 2D-diffusion rate dF/dt in (Jr,cos I) space
function dFdt2D_JrcosI(Jr::Float64, cosI::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbL::Int64=50, Lmax::Float64=3.0, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default,
            nbvartheta::Int64=nbvartheta_default, nbphi::Int64=nbphi_default,
            nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    sumJr = 0.0
    sumcosI = 0.0

    for i=1:nbL
        L = Lmax/nbL*(i-0.5)

        

        fJr_p, _ = flux2D_JrL(Jr+eps,L,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,eps,m_test)
        fJr_m, _ = flux2D_JrL(Jr-eps,L,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,eps,m_test)

        dJr = (fJr_p-fJr_m)/(2.0*eps)

        fcosI_p  = FluxCosI(Jr,L,cosI+eps,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
        fcosI_m  = FluxCosI(Jr,L,cosI-eps,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)


        dcosI = (fcosI_p-fcosI_m)/(2.0*eps)

        sumJr += dJr
        sumcosI  += dcosI
    end


    sumJr *= Lmax/nbL
    sumcosI *= Lmax/nbL

    return -(sumJr+sumcosI)
end

# Computes the 2D-diffusion rate dF/dt in (L,cos I) space
function dFdt2D_LcosI(L::Float64, cosI::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbJr::Int64=50, Jrmax::Float64=3.0, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default,
            nbvartheta::Int64=nbvartheta_default, nbphi::Int64=nbphi_default,
            nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

    sumL = 0.0
    sumcosI = 0.0

    for i=1:nbJr
        Jr = Jrmax/nbJr*(i-0.5)

        _ , fL_p = flux2D_JrL(Jr,L+eps,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,eps,m_test)
        _ , fL_m = flux2D_JrL(Jr,L-eps,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,eps,m_test)

        dL = (fL_p-fL_m)/(2.0*eps)

        fcosI_p  = FluxCosI(Jr,L,cosI+eps,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
        fcosI_m  = FluxCosI(Jr,L,cosI-eps,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)


        dcosI = (fcosI_p-fcosI_m)/(2.0*eps)

        sumL += dL
        sumcosI  += dcosI
    end

    sumL *= Jrmax/nbJr
    sumcosI *= Jrmax/nbJr

    return -(sumL+sumcosI)
end


# Drift (orbit-averaged) in cos I space
function orbitAverageDriftCosI(sp::Float64, sa::Float64, cosI::Float64,
    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
    nbw::Int64=nbw_default,
    nbvartheta::Int64=nbvartheta_default, nbphi::Int64=nbphi_default,
    m_test::Float64=m_field)

drift = 0.0

halfperiod = 0.0

E, L = E_L_from_sp_sa(sp,sa)


sma, ecc = sma_ecc_from_sp_sa(sp,sa)


for iu=1:nbAvr
uloc = -1+2*(iu-0.5)/nbAvr
sloc = s_from_u_sma_ecc(uloc,sma,ecc)
rloc = r_from_s(sloc)
jac_loc = Theta(uloc,sp,sa)

vr = sqrt(2.0*abs(E - psiEff(rloc,L)))
vt = L/rloc

halfperiod += jac_loc

driftloc = localDriftCosIAngleAverage(rloc,vr,vt,cosI,m_field,alpha,nbw,nbvartheta,nbphi,m_test)

drift += jac_loc*driftloc

end
drift /= halfperiod


return drift
end



# Computes the 1D-flux in cos I space.
function FluxCosI(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
    alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
    nbw::Int64=nbw_default,
    nbvartheta::Int64=nbvartheta_default, nbphi::Int64=nbphi_default,
    nbu::Int64=nbu0, m_test::Float64=m_field)

E = E_from_Jr_L(Jr,L,nbu)
if (Jr > 0.0)
sp, sa = sp_sa_from_E_L(E,L)
else
sc = _sc(E/_E0)
sp, sa = sc, sc
end
sma, ecc = sma_ecc_from_sp_sa(sp,sa)


drift = orbitAverageDriftCosI(sp,sa,cosI,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,m_test)

Frot = _Frot_cosI(E,L,cosI,alpha)

return drift*Frot
end



###########################################################
###########################################################

# Computes the 3D-diffusion rate dF/dt in (Jr,L,Lz) space
function dFdt3D(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
            nbw::Int64=nbw_default, nbvartheta::Int64=nbvartheta_default, nbphi::Int64=nbphi_default,
            nbu::Int64=nbu0, eps::Float64=10^(-5), m_test::Float64=m_field)

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

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffs(Jr,L,Lz/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot= _Frot(E,L,Lz,alpha)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffs(Jr_p,L,Lz/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Jrp = _Frot(E_Jrp,L,Lz,alpha)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffs(Jr_m,L,Lz/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Jrm = _Frot(E_Jrm,L,Lz,alpha)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffs(Jr,L_p,Lz/L_p,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Lp = _Frot(E_Lp,L_p,Lz,alpha)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffs(Jr,L_m,Lz/L_m,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Lm = _Frot(E_Lm,L_m,Lz,alpha)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffs(Jr,L,Lz_p/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Lzp = _Frot(E,L,Lz_p,alpha)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffs(Jr,L,Lz_m/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_Lzm = _Frot(E,L,Lz_m,alpha)


    DJrF_Jr = (dJr_Jrp*Frot_Jrp-dJr_Jrm*Frot_Jrm)/(2.0*eps*_L0)
    DLF_L = (dL_Lp*Frot_Lp-dL_Lm*Frot_Lm)/(2.0*eps*_L0)
    DLzF_Lz = (dLz_Lzp*Frot_Lzp-dLz_Lzm*Frot_Lzm)/(2.0*eps*_L0)

    DJrJr_F_JrJr = (dJrJr_Jrp*Frot_Jrp + dJrJr_Jrm*Frot_Jrm - 2.0*dJrJr*Frot)/(eps*_L0)^2
    DLL_F_LL = (dLL_Lp*Frot_Lp + dLL_Lm*Frot_Lm - 2.0*dLL*Frot)/(eps*_L0)^2
    DLzLz_F_LzLz = (dLzLz_Lzp*Frot_Lzp + dLzLz_Lzm*Frot_Lzm - 2.0*dLzLz*Frot)/(eps*_L0)^2





    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffs(Jr_p,L_p,Lz/L_p,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_JrpLp = _Frot(E_Jrp_Lp,L_p,Lz,alpha)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffs(Jr_m,L_m,Lz/L_m,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_JrmLm = _Frot(E_Jrm_Lm,L_m,Lz,alpha)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffs(Jr_p,L_m,Lz/L_m,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_JrpLm = _Frot(E_Jrp_Lm,L_m,Lz,alpha)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffs(Jr_m,L_p,Lz/L_p,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_JrmLp = _Frot(E_Jrm_Lp,L_p,Lz,alpha)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffs(Jr_p,L,Lz_p/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_JrpLzp = _Frot(E_Jrp,L,Lz_p,alpha)
    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffs(Jr_m,L,Lz_m/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_JrmLzm = _Frot(E_Jrm,L,Lz_m,alpha)

    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffs(Jr_p,L,Lz_m/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_JrpLzm = _Frot(E_Jrp,L,Lz_m,alpha)
    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffs(Jr_m,L,Lz_p/L,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_JrmLzp = _Frot(E_Jrm,L,Lz_p,alpha)

    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffs(Jr,L_p,Lz_p/L_p,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_LpLzp = _Frot(E_Lp,L_p,Lz_p,alpha)
    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffs(Jr,L_m,Lz_m/L_m,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_LmLzm = _Frot(E_Lm,L_m,Lz_m,alpha)

    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffs(Jr,L_p,Lz_m/L_p,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_LpLzm = _Frot(E_Lp,L_p,Lz_m,alpha)
    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffs(Jr,L_m,Lz_p/L_m,m_field,alpha,nbAvr,nbw,nbvartheta,nbphi,nbu,m_test)
    Frot_LmLzp =_Frot(E_Lm,L_m,Lz_p,alpha)




    DJrL_F_JrL = (dJrL_JrpLp*Frot_JrpLp+dJrL_JrmLm*Frot_JrmLm-dJrL_JrpLm*Frot_JrpLm-dJrL_JrmLp*Frot_JrmLp)/(2.0*eps*_L0)^2

    DJrLz_F_JrLz = (dJrLz_JrpLzp*Frot_JrpLzp+dJrLz_JrmLzm*Frot_JrmLzm-dJrLz_JrpLzm*Frot_JrpLzm-dJrLz_JrmLzp*Frot_JrmLzp)/(2.0*eps*_L0)^2

    DLLz_F_LLz = (dLLz_LpLzp*Frot_LpLzp+dLLz_LmLzm*Frot_LmLzm-dLLz_LpLzm*Frot_LpLzm-dLLz_LmLzp*Frot_LmLzp)/(2.0*eps*_L0)^2


    fluxJr_Jr = DJrF_Jr - 0.5*(DJrJr_F_JrJr + DJrL_F_JrL + DJrLz_F_JrLz)
    fluxL_L = DLF_L - 0.5*(DJrL_F_JrL + DLL_F_LL + DLLz_F_LLz)
    fluxLz_Lz = DLzF_Lz - 0.5*(DJrLz_F_JrLz + DLLz_F_LLz + DLzLz_F_LzLz)

    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz)

end
