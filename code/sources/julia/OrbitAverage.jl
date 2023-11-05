##################################################
# Local diffusion coefficients
##################################################

# USED
function localOrbitChangeAngleAverageExact_nbKseparate(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbw::Int64=nbw_default,
                        nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                        m_test::Float64=m_field)

    v = sqrt(vr^2+vt^2)
    vr_v = vr/v
    vt_v = vt/v

    sinI = sqrt(abs(1.0 - cosI^2))

    dvPar, dvPar2, dvPerp2, sinSqdvPerSq = localVelChange3DAngleAverageExact_nbKseparate(r,vr,vt,cosI,m_field,alpha,nbw,nbvarphi,nbphi,m_test)

    dE   = 0.5*dvPar2 + 0.5*dvPerp2 + v*dvPar
    dE2  = v^2* dvPar2
    dL   = r*vt_v*dvPar + 0.25*(r/vt)*dvPerp2
    dL2  = r^2*vt_v^2*dvPar2 + 0.5*r^2*vr_v^2*dvPerp2
    dEdL = r*vt*dvPar2

    # Lz coefficients

    dLz = r*vt_v*cosI*dvPar
    dLz2 = r^2*cosI^2*(vt_v^2*dvPar2 + 0.5*vr_v^2*dvPerp2) + 0.5*r^2*sinSqdvPerSq*sinI^2
    dEdLz = r*vt*cosI*dvPar2
    dLdLz = r^2*cosI*(vt_v^2*dvPar2 + 0.5*vr_v^2*dvPerp2)

    return dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz
end


##################################################
# Orbit-averaged (E,L,Lz)-diffusion coefficients
##################################################



# USED
function orbitAverageEnergyCoeffsOptiExact_spsa_nbKseparate(sp::Float64, sa::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                    nbw::Int64=nbw_default,
                    nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                    m_test::Float64=m_field)

    avrDE = 0.0
    avrDL = 0.0
    avrDEE = 0.0
    avrDEL = 0.0
    avrDLL = 0.0
    avrDLz = 0.0
    avrDLzLz = 0.0
    avrDLLz = 0.0
    avrDELz = 0.0

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


        dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChangeAngleAverageExact_nbKseparate(rloc,vr,vt,cosI,m_field,alpha,nbw,nbvarphi,nbphi,m_test)

        avrDE += jac_loc*dE
        avrDL += jac_loc*dL
        avrDEE += jac_loc*dE2
        avrDEL += jac_loc*dEdL
        avrDLL += jac_loc*dL2

        avrDLz += jac_loc*dLz
        avrDLzLz += jac_loc*dLz2
        avrDLLz += jac_loc*dLdLz
        avrDELz += jac_loc*dEdLz





    end
    avrDE /= halfperiod
    avrDL /= halfperiod
    avrDEE /= halfperiod
    avrDEL /= halfperiod
    avrDLL /= halfperiod
    avrDLz /= halfperiod
    avrDLzLz /= halfperiod
    avrDLLz /= halfperiod
    avrDELz /= halfperiod

    return avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz
end



# USED
function orbitAverageDriftCosIOptiExact_spsa_nbKseparate(sp::Float64, sa::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                    nbw::Int64=nbw_default,
                    nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
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




        driftloc = localDriftCosIAngleAverageExact_nbKseparate(rloc,vr,vt,cosI,m_field,alpha,nbw,nbvarphi,nbphi,m_test)

        drift += jac_loc*driftloc






    end
    drift /= halfperiod


    return drift
end


##################################################
# Orbit-averaged (Jr,L,Lz)-diffusion coefficients
##################################################


# USED
function orbitAverageActionCoeffsOptiExact_nbKseparate(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                                nbw::Int64=nbw_default,
                                nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                                nbu::Int64=nbu0, m_test::Float64=m_field)

    E = E_from_Jr_L(Jr,L,nbu)
    if (Jr > 0.0)
        sp, sa = sp_sa_from_E_L(E,L)
    else
        sc = _sc(E/_E0)
        sp, sa = sc, sc
    end
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOptiExact_spsa_nbKseparate(sp,sa,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,m_test)




    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end



# USED
function FluxCosIOptiExact_nbKseparate(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                                nbw::Int64=nbw_default,
                                nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                                nbu::Int64=nbu0, m_test::Float64=m_field)

    E = E_from_Jr_L(Jr,L,nbu)
    #sp, sa = sp_sa_from_E_L(E,L)
    if (Jr > 0.0)
        sp, sa = sp_sa_from_E_L(E,L)
    else
        sc = _sc(E/_E0)
        #println(sc)
        sp, sa = sc, sc
    end
    #println((sp, sa))
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    #avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOptiExact(E,L,cosI,m_field,alpha,nbAvr,nbK,m_test)
    drift = orbitAverageDriftCosIOptiExact_spsa_nbKseparate(sp,sa,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,m_test)

    Frot = _Frot_cosI(E,L,cosI,alpha)

    return drift*Frot
end



#####################################################################
# (Jr,L,cos I) space
#####################################################################

# ok
function orbitAverageActionCoeffs_cosI_OptiExact_nbKseparate(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                                nbw::Int64=nbw_default,
                                nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                                 nbu::Int64=nbu0, m_test::Float64=m_field)


    avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz = orbitAverageActionCoeffsOptiExact_nbKseparate(Jr,L,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,nbu,m_test)

    avrDI = -cosI/L*avrDL + 1.0/L*avrDLz + cosI/L^2*avrDLL- 1.0/L^2*avrDLLz
    avrDJrI = -cosI/L*avrDJrL + 1.0/L*avrDJrLz
    avrDLI = -cosI/L*avrDLL + 1.0/L*avrDLLz
    avrDII = cosI^2/L^2*avrDLL - 2.0*cosI/L^2*avrDLLz + 1.0/L^2*avrDLzLz

    return avrDJr, avrDL, avrDI, avrDJrJr, avrDLL, avrDII, avrDJrL, avrDJrI, avrDLI
end
