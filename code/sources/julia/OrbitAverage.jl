using PolynomialRoots

#
# Redo everything (at least check functions) with tex_RR_Anisotropic function
# (See OrbitalParameters.jl, etc)
#

##################################################
# Local diffusion coefficients
##################################################

function localOrbitChange(r::Float64, theta::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbK::Int64=nbK_default, m_test::Float64=m_field)

    v = sqrt(vr^2+vt^2)
    vr_v = vr/v
    vt_v = vt/v

    sinI = sqrt(abs(1.0 - cosI^2))

    dvPar, dvPar2, dvPerp2 = localVelChange3D(r,theta,vr,vt,cosI,m_field,alpha,nbK,m_test)

    dE   = 0.5*dvPar2 + 0.5*dvPerp2 + v*dvPar
    dE2  = v^2* dvPar2
    dL   = r*vt_v*dvPar + 0.25*(r/vt)*dvPerp2
    dL2  = r^2*vt_v^2*dvPar2 + 0.5*r^2*vr_v^2*dvPerp2
    dEdL = r*vt*dvPar2

    # Lz coefficients

    dLz = r*vt_v*cosI*dvPar
    dLz2 = r^2*cosI^2*(vt_v^2*dvPar2 + 0.5*vr_v^2*dvPerp2) + 0.5*r^2*sin(theta)^2*sinI^2*dvPerp2
    dEdLz = r*vt*cosI*dvPar2
    dLdLz = r^2*cosI*(vt_v^2*dvPar2 + 0.5*vr_v^2*dvPerp2)

    return dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz
end

# Orbit-averaged over theta
# No rot
function localOrbitChangeNoRot(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64,
                        nbK::Int64=nbK_default, m_test::Float64=m_field)

    v = sqrt(vr^2+vt^2)
    vr_v = vr/v
    vt_v = vt/v

    sinI = sqrt(abs(1.0 - cosI^2))

    dvPar, dvPar2, dvPerp2 = localVelChange3D(r,0.0,vr,vt,cosI,m_field,0.0,nbK,m_test)

    dE   = 0.5*dvPar2 + 0.5*dvPerp2 + v*dvPar
    dE2  = v^2* dvPar2
    dL   = r*vt_v*dvPar + 0.25*(r/vt)*dvPerp2
    dL2  = r^2*vt_v^2*dvPar2 + 0.5*r^2*vr_v^2*dvPerp2
    dEdL = r*vt*dvPar2

    # Lz coefficients

    dLz = r*vt_v*cosI*dvPar
    dLz2 = r^2*cosI^2*(vt_v^2*dvPar2 + 0.5*vr_v^2*dvPerp2) + 0.25*r^2*sinI^2*dvPerp2
    dEdLz = r*vt*cosI*dvPar2
    dLdLz = r^2*cosI*(vt_v^2*dvPar2 + 0.5*vr_v^2*dvPerp2)

    return dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz
end

# Orbit-averaged over theta
# With rotation
function localOrbitChangeAngleAverage(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbAvrTh::Int64=nbAvrTh_default, nbK::Int64=nbK_default,
                        m_test::Float64=m_field)

    v = sqrt(vr^2+vt^2)
    vr_v = vr/v
    vt_v = vt/v

    sinI = sqrt(abs(1.0 - cosI^2))

    dvPar, dvPar2, dvPerp2, sinSqdvPerSq = localVelChange3DAngleAverage(r,vr,vt,cosI,m_field,alpha,nbAvrTh,nbK,m_test)

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

function localOrbitChangeAngleAverage_nbKseparate(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbAvrTh::Int64=nbAvrTh_default, nbw::Int64=nbw_default,
                        nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                        m_test::Float64=m_field)

    v = sqrt(vr^2+vt^2)
    vr_v = vr/v
    vt_v = vt/v

    sinI = sqrt(abs(1.0 - cosI^2))

    dvPar, dvPar2, dvPerp2, sinSqdvPerSq = localVelChange3DAngleAverage_nbKseparate(r,vr,vt,cosI,m_field,alpha,nbAvrTh,nbw,nbvarphi,nbphi,m_test)

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


function localOrbitChangeAngleAverageExact(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbK::Int64=nbK_default,
                        m_test::Float64=m_field)

    v = sqrt(vr^2+vt^2)
    vr_v = vr/v
    vt_v = vt/v

    sinI = sqrt(abs(1.0 - cosI^2))

    dvPar, dvPar2, dvPerp2, sinSqdvPerSq = localVelChange3DAngleAverageExact(r,vr,vt,cosI,m_field,alpha,nbK,m_test)

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

function radius_s_bounds(E::Float64, L::Float64)
    if (E >= 0.0)
        return "Unbounded orbit"
    elseif (L > _Lc(E))
        return "Not a possible orbit"
    else
        tE = E/_E0
        tL = L/_L0
        if (L != 0.0)
            rts = sort(real(roots([1,-(tE-tL^2/2),-1,tE])))
            return Float64(rts[2]), Float64(rts[3])
        else
            return 1.0, 1/tE
        end
    end
end


function orbitAverageEnergyCoeffs(E::Float64, L::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                    nbAvrTh::Int64=nbAvrTh_default, nbK::Int64=nbK_default, m_test::Float64=m_field)

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

    sp, sa = radius_s_bounds(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        sloc = s_from_u_sma_ecc(uloc,sma,ecc)
        rloc = r_from_s(sloc)
        jac_loc = Theta(uloc,sp,sa)

        vr = sqrt(2.0*(E - psiEff(rloc,L)))
        vt = L/rloc

        halfperiod += jac_loc


        for ith=1:nbAvrTh
            thloc = 2.0*pi*(ith-0.5)/nbAvrTh



            # Positive vr contribution

            dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChange(rloc,thloc,vr,vt,cosI,m_field,alpha,nbK,m_test)

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


    end
    avrDE /= halfperiod*nbAvrTh
    avrDL /= halfperiod*nbAvrTh
    avrDEE /= halfperiod*nbAvrTh
    avrDEL /= halfperiod*nbAvrTh
    avrDLL /= halfperiod*nbAvrTh
    avrDLz /= halfperiod*nbAvrTh
    avrDLzLz /= halfperiod*nbAvrTh
    avrDLLz /= halfperiod*nbAvrTh
    avrDELz /= halfperiod*nbAvrTh

    return avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz
end

# Angle-averaged first in a clever way
function orbitAverageEnergyCoeffsOpti(E::Float64, L::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                    nbAvrTh::Int64=nbAvrTh_default, nbK::Int64=nbK_default, m_test::Float64=m_field)

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

    sp, sa = radius_s_bounds(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        sloc = s_from_u_sma_ecc(uloc,sma,ecc)
        rloc = r_from_s(sloc)
        jac_loc = Theta(uloc,sp,sa)

        vr = sqrt(2.0*(E - psiEff(rloc,L)))
        vt = L/rloc

        halfperiod += jac_loc




        dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChangeAngleAverage(rloc,vr,vt,cosI,m_field,alpha,nbAvrTh,nbK,m_test)

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


function orbitAverageEnergyCoeffsOptiExact(E::Float64, L::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                    nbK::Int64=nbK_default, m_test::Float64=m_field)

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

    sp, sa = radius_s_bounds(E,L)
    #println("ooooo :",(sp,sa))
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        sloc = s_from_u_sma_ecc(uloc,sma,ecc)
        rloc = r_from_s(sloc)
        jac_loc = Theta(uloc,sp,sa)

        vr = sqrt(2.0*(E - psiEff(rloc,L)))
        vt = L/rloc

        halfperiod += jac_loc




        dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChangeAngleAverageExact(rloc,vr,vt,cosI,m_field,alpha,nbK,m_test)

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



function orbitAverageEnergyCoeffsOptiExact_spsa(sp::Float64, sa::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                    nbK::Int64=nbK_default, m_test::Float64=m_field)

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




        dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChangeAngleAverageExact(rloc,vr,vt,cosI,m_field,alpha,nbK,m_test)

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


function orbitAverageEnergyCoeffsOpti_spsa_nbKseparate(sp::Float64, sa::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
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




        dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChangeAngleAverage_nbKseparate(rloc,vr,vt,cosI,m_field,alpha,nbAvrTh,nbw,nbvarphi,nbphi,m_test)

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




# No rotation: only integral over "radius"
function orbitAverageEnergyCoeffsNoRot(E::Float64, L::Float64, cosI::Float64,
                    m_field::Float64, nbAvr::Int64=nbAvr_default,
                    nbK::Int64=nbK_default, m_test::Float64=m_field)

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

    sp, sa = radius_s_bounds(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        sloc = s_from_u_sma_ecc(uloc,sma,ecc)
        rloc = r_from_s(sloc)
        jac_loc = Theta(uloc,sp,sa)

        vr = sqrt(2.0*(E - psiEff(rloc,L)))
        vt = L/rloc

        halfperiod += jac_loc

        # Positive vr contribution

        dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChangeNoRot(rloc,vr,vt,cosI,m_field,nbK,m_test)

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


function orbitAverageEnergyCoeffsNoRotPar(E::Float64, L::Float64, cosI::Float64,
                    m_field::Float64, nbAvr::Int64=nbAvr_default,
                    nbK::Int64=nbK_default, m_test::Float64=m_field)

    avrDE = Threads.Atomic{Float64}(0.0)
    avrDL = Threads.Atomic{Float64}(0.0)
    avrDEE = Threads.Atomic{Float64}(0.0)
    avrDEL = Threads.Atomic{Float64}(0.0)
    avrDLL = Threads.Atomic{Float64}(0.0)
    avrDLz = Threads.Atomic{Float64}(0.0)
    avrDLzLz = Threads.Atomic{Float64}(0.0)
    avrDLLz = Threads.Atomic{Float64}(0.0)
    avrDELz = Threads.Atomic{Float64}(0.0)

    halfperiod = Threads.Atomic{Float64}(0.0)

    sp, sa = radius_s_bounds(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    Threads.@threads for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        sloc = s_from_u_sma_ecc(uloc,sma,ecc)
        rloc = r_from_s(sloc)
        jac_loc = Theta(uloc,sp,sa)

        vr = sqrt(2.0*(E - psiEff(rloc,L)))
        vt = L/rloc



        # Positive vr contribution

        dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChangeNoRot(rloc,vr,vt,cosI,m_field,nbK,m_test)


        Threads.atomic_add!(avrDE,jac_loc*dE)
        Threads.atomic_add!(avrDL,jac_loc*dL)
        Threads.atomic_add!(avrDEE,jac_loc*dE2)
        Threads.atomic_add!(avrDEL,jac_loc*dEdL)
        Threads.atomic_add!(avrDLL,jac_loc*dL2)
        Threads.atomic_add!(avrDLz,jac_loc*dLz)
        Threads.atomic_add!(avrDLzLz,jac_loc*dLz2)
        Threads.atomic_add!(avrDLLz,jac_loc*dLdLz)
        Threads.atomic_add!(avrDELz,jac_loc*dEdLz)

        Threads.atomic_add!(halfperiod,jac_loc)





    end
    avrDE[] /= halfperiod[]
    avrDL[] /= halfperiod[]
    avrDEE[] /= halfperiod[]
    avrDEL[] /= halfperiod[]
    avrDLL[] /= halfperiod[]
    avrDLz[] /= halfperiod[]
    avrDLzLz[] /= halfperiod[]
    avrDLLz[] /= halfperiod[]
    avrDELz[] /= halfperiod[]


    return avrDE[], avrDL[], avrDLz[], avrDEE[], avrDLL[], avrDLzLz[], avrDEL[], avrDELz[], avrDLLz[]
end



function orbitAverageEnergyCoeffsPar(E::Float64, L::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                    nbAvrTh::Int64=nbAvrTh_default, nbK::Int64=nbK_default, m_test::Float64=m_field)


    avrDE = Threads.Atomic{Float64}(0.0)
    avrDL = Threads.Atomic{Float64}(0.0)
    avrDEE = Threads.Atomic{Float64}(0.0)
    avrDEL = Threads.Atomic{Float64}(0.0)
    avrDLL = Threads.Atomic{Float64}(0.0)
    avrDLz = Threads.Atomic{Float64}(0.0)
    avrDLzLz = Threads.Atomic{Float64}(0.0)
    avrDLLz = Threads.Atomic{Float64}(0.0)
    avrDELz = Threads.Atomic{Float64}(0.0)

    halfperiod = Threads.Atomic{Float64}(0.0)

    sp, sa = radius_s_bounds(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)



    Threads.@threads for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        sloc = s_from_u_sma_ecc(uloc,sma,ecc)
        rloc = r_from_s(sloc)
        jac_loc = Theta(uloc,sp,sa)

        vr = sqrt(2.0*(E - psiEff(rloc,L)))
        vt = L/rloc

        avrDE_th = 0.0
        avrDL_th = 0.0
        avrDEE_th = 0.0
        avrDEL_th = 0.0
        avrDLL_th = 0.0
        avrDLz_th = 0.0
        avrDLzLz_th = 0.0
        avrDLLz_th = 0.0
        avrDELz_th = 0.0

        Threads.atomic_add!(halfperiod,jac_loc)

        for ith=1:nbAvrTh
            thloc = 2.0*pi*(ith-0.5)/nbAvrTh



            # Positive vr contribution

            dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChange(rloc,thloc,vr,vt,cosI,m_field,alpha,nbK,m_test)

            avrDE_th += jac_loc*dE
            avrDL_th += jac_loc*dL
            avrDEE_th += jac_loc*dE2
            avrDEL_th += jac_loc*dEdL
            avrDLL_th += jac_loc*dL2

            avrDLz_th += jac_loc*dLz
            avrDLzLz_th += jac_loc*dLz2
            avrDLLz_th += jac_loc*dLdLz
            avrDELz_th += jac_loc*dEdLz

        end

        Threads.atomic_add!(avrDE,avrDE_th)
        Threads.atomic_add!(avrDL,avrDL_th)
        Threads.atomic_add!(avrDEE,avrDEE_th)
        Threads.atomic_add!(avrDEL,avrDEL_th)
        Threads.atomic_add!(avrDLL,avrDLL_th)
        Threads.atomic_add!(avrDLz,avrDLz_th)
        Threads.atomic_add!(avrDLzLz,avrDLzLz_th)
        Threads.atomic_add!(avrDLLz,avrDLLz_th)
        Threads.atomic_add!(avrDELz,avrDELz_th)




    end
    avrDE[] /= halfperiod[]*nbAvrTh
    avrDL[] /= halfperiod[]*nbAvrTh
    avrDEE[] /= halfperiod[]*nbAvrTh
    avrDEL[] /= halfperiod[]*nbAvrTh
    avrDLL[] /= halfperiod[]*nbAvrTh
    avrDLz[] /= halfperiod[]*nbAvrTh
    avrDLzLz[] /= halfperiod[]*nbAvrTh
    avrDLLz[] /= halfperiod[]*nbAvrTh
    avrDELz[] /= halfperiod[]*nbAvrTh

    return avrDE[], avrDL[], avrDLz[], avrDEE[], avrDLL[], avrDLzLz[], avrDEL[], avrDELz[], avrDLLz[]
end

function orbitAverageEnergyCoeffsOptiPar(E::Float64, L::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                    nbAvrTh::Int64=nbAvrTh_default, nbK::Int64=nbK_default, m_test::Float64=m_field)

    avrDE = Threads.Atomic{Float64}(0.0)
    avrDL = Threads.Atomic{Float64}(0.0)
    avrDEE = Threads.Atomic{Float64}(0.0)
    avrDEL = Threads.Atomic{Float64}(0.0)
    avrDLL = Threads.Atomic{Float64}(0.0)
    avrDLz = Threads.Atomic{Float64}(0.0)
    avrDLzLz = Threads.Atomic{Float64}(0.0)
    avrDLLz = Threads.Atomic{Float64}(0.0)
    avrDELz = Threads.Atomic{Float64}(0.0)

    halfperiod = Threads.Atomic{Float64}(0.0)

    sp, sa = radius_s_bounds(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)



    Threads.@threads for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        sloc = s_from_u_sma_ecc(uloc,sma,ecc)
        rloc = r_from_s(sloc)
        jac_loc = Theta(uloc,sp,sa)

        vr = sqrt(2.0*(E - psiEff(rloc,L)))
        vt = L/rloc

        avrDE_th = 0.0
        avrDL_th = 0.0
        avrDEE_th = 0.0
        avrDEL_th = 0.0
        avrDLL_th = 0.0
        avrDLz_th = 0.0
        avrDLzLz_th = 0.0
        avrDLLz_th = 0.0
        avrDELz_th = 0.0

        Threads.atomic_add!(halfperiod,jac_loc)





        # Positive vr contribution

        dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChangeAngleAverage(rloc,vr,vt,cosI,m_field,alpha,nbAvrTh,nbK,m_test)


        Threads.atomic_add!(avrDE,jac_loc*dE)
        Threads.atomic_add!(avrDL,jac_loc*dL)
        Threads.atomic_add!(avrDEE,jac_loc*dE2)
        Threads.atomic_add!(avrDEL,jac_loc*dEdL)
        Threads.atomic_add!(avrDLL,jac_loc*dL2)
        Threads.atomic_add!(avrDLz,jac_loc*dLz)
        Threads.atomic_add!(avrDLzLz,jac_loc*dLz2)
        Threads.atomic_add!(avrDLLz,jac_loc*dLdLz)
        Threads.atomic_add!(avrDELz,jac_loc*dEdLz)




    end
    avrDE[] /= halfperiod[]
    avrDL[] /= halfperiod[]
    avrDEE[] /= halfperiod[]
    avrDEL[] /= halfperiod[]
    avrDLL[] /= halfperiod[]
    avrDLz[] /= halfperiod[]
    avrDLzLz[] /= halfperiod[]
    avrDLLz[] /= halfperiod[]
    avrDELz[] /= halfperiod[]

    return avrDE[], avrDL[], avrDLz[], avrDEE[], avrDLL[], avrDLzLz[], avrDEL[], avrDELz[], avrDLLz[]
end

function orbitAverageEnergyCoeffsOptiExactPar(E::Float64, L::Float64, cosI::Float64,
                    m_field::Float64, alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                    nbK::Int64=nbK_default, m_test::Float64=m_field)

    avrDE = Threads.Atomic{Float64}(0.0)
    avrDL = Threads.Atomic{Float64}(0.0)
    avrDEE = Threads.Atomic{Float64}(0.0)
    avrDEL = Threads.Atomic{Float64}(0.0)
    avrDLL = Threads.Atomic{Float64}(0.0)
    avrDLz = Threads.Atomic{Float64}(0.0)
    avrDLzLz = Threads.Atomic{Float64}(0.0)
    avrDLLz = Threads.Atomic{Float64}(0.0)
    avrDELz = Threads.Atomic{Float64}(0.0)

    halfperiod = Threads.Atomic{Float64}(0.0)

    sp, sa = radius_s_bounds(E,L)
    #println("(sp,sa) = ",(sp,sa))
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)



    Threads.@threads for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        sloc = s_from_u_sma_ecc(uloc,sma,ecc)
        rloc = r_from_s(sloc)
        jac_loc = Theta(uloc,sp,sa)

        vr = sqrt(2.0*abs(E - psiEff(rloc,L)))
        vt = L/rloc

        avrDE_th = 0.0
        avrDL_th = 0.0
        avrDEE_th = 0.0
        avrDEL_th = 0.0
        avrDLL_th = 0.0
        avrDLz_th = 0.0
        avrDLzLz_th = 0.0
        avrDLLz_th = 0.0
        avrDELz_th = 0.0

        Threads.atomic_add!(halfperiod,jac_loc)





        # Positive vr contribution

        dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz = localOrbitChangeAngleAverageExact(rloc,vr,vt,cosI,m_field,alpha,nbK,m_test)


        Threads.atomic_add!(avrDE,jac_loc*dE)
        Threads.atomic_add!(avrDL,jac_loc*dL)
        Threads.atomic_add!(avrDEE,jac_loc*dE2)
        Threads.atomic_add!(avrDEL,jac_loc*dEdL)
        Threads.atomic_add!(avrDLL,jac_loc*dL2)
        Threads.atomic_add!(avrDLz,jac_loc*dLz)
        Threads.atomic_add!(avrDLzLz,jac_loc*dLz2)
        Threads.atomic_add!(avrDLLz,jac_loc*dLdLz)
        Threads.atomic_add!(avrDELz,jac_loc*dEdLz)




    end
    avrDE[] /= halfperiod[]
    avrDL[] /= halfperiod[]
    avrDEE[] /= halfperiod[]
    avrDEL[] /= halfperiod[]
    avrDLL[] /= halfperiod[]
    avrDLz[] /= halfperiod[]
    avrDLzLz[] /= halfperiod[]
    avrDLLz[] /= halfperiod[]
    avrDELz[] /= halfperiod[]

    return avrDE[], avrDL[], avrDLz[], avrDEE[], avrDLL[], avrDLzLz[], avrDEL[], avrDELz[], avrDLLz[]
end


##################################################
# Orbit-averaged (Jr,L,Lz)-diffusion coefficients
##################################################



function orbitAverageActionCoeffs(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    E = E_from_Jr_L(Jr,L,nbu)
    sp, sa = sp_sa_from_E_L(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffs(E,L,cosI,m_field,alpha,nbAvr,nbAvrTh,nbK,m_test)


    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end

function orbitAverageActionCoeffsOpti(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    E = E_from_Jr_L(Jr,L,nbu)
    sp, sa = sp_sa_from_E_L(E,L)

    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOpti(E,L,cosI,m_field,alpha,nbAvr,nbAvrTh,nbK,m_test)


    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end

function orbitAverageActionCoeffsOptiExact(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

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
    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOptiExact_spsa(sp,sa,cosI,m_field,alpha,nbAvr,nbK,m_test)




    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end

# USED
function orbitAverageActionCoeffsOptiExact_nbKseparate(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
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
    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOptiExact_spsa_nbKseparate(sp,sa,cosI,m_field,alpha,nbAvr,nbw,nbvarphi,nbphi,m_test)




    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end


function orbitAverageActionCoeffsOpti_nbKseparate(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
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
    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOpti_spsa_nbKseparate(sp,sa,cosI,m_field,alpha,nbAvr,nbAvrTh,nbw,nbvarphi,nbphi,m_test)




    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

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



function orbitAverageActionCoeffsOptiExactOld(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

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


    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOptiExact(E,L,cosI,m_field,alpha,nbAvr,nbK,m_test)




    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end




function orbitAverageActionCoeffsPar(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    E = E_from_Jr_L(Jr,L,nbu)
    sp, sa = sp_sa_from_E_L(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsPar(E,L,cosI,m_field,alpha,nbAvr,nbAvrTh,nbK,m_test)


    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end


function orbitAverageActionCoeffsOptiPar(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    E = E_from_Jr_L(Jr,L,nbu)
    #println((Jr,E,L,cosI))
    sp, sa = sp_sa_from_E_L(E,L)
    #println((sp,sa))
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)
    #println((sma,ecc))


    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOptiPar(E,L,cosI,m_field,alpha,nbAvr,nbAvrTh,nbK,m_test)


    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end

function orbitAverageActionCoeffsOptiExactPar(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    # E = E_from_Jr_L(Jr,L,nbu)
    # #println((Jr,E,L,cosI))
    # sp, sa = sp_sa_from_E_L(E,L)
    # #println((sp,sa))
    # sma, ecc = sma_ecc_from_sp_sa(sp,sa)
    # #println((sma,ecc))
    #
    #
    # avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOptiExactPar(E,L,cosI,m_field,alpha,nbAvr,nbK,m_test)

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
    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsOptiExactPar_spsa(sp,sa,cosI,m_field,alpha,nbAvr,nbK,m_test)





    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end

function orbitAverageActionCoeffsNoRot(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                nbAvr::Int64=nbAvr_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    E = E_from_Jr_L(Jr,L,nbu)
    sp, sa = sp_sa_from_E_L(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsNoRot(E,L,cosI,m_field,nbAvr,nbK,m_test)


    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end

function orbitAverageActionCoeffsNoRotPar(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                nbAvr::Int64=nbAvr_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)

    E = E_from_Jr_L(Jr,L,nbu)
    sp, sa = sp_sa_from_E_L(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    avrDE, avrDL, avrDLz, avrDEE, avrDLL, avrDLzLz, avrDEL, avrDELz, avrDLLz = orbitAverageEnergyCoeffsNoRotPar(E,L,cosI,m_field,nbAvr,nbK,m_test)


    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL
    avrDJrLz = dJrdE*avrDELz+dJrdL*avrDLLz

    return avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz
end

#####################################################################
# (Jr,L,cos I) space
#####################################################################


function orbitAverageActionCoeffs_cosI_OptiExact(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)


    avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz = orbitAverageActionCoeffsOptiExact(Jr,L,cosI,m_field,alpha,nbAvr,nbK,nbu,m_test)

    avrDI = -cosI/L*avrDL + 1.0/L*avrDLz + cosI/L^2*avrDLL- 1.0/L^2*avrDLLz
    avrDJrI = -cosI/L*avrDJrL + 1.0/L*avrDJrLz
    avrDLI = -cosI/L*avrDLL + 1.0/L*avrDLLz
    avrDII = cosI^2/L^2*avrDLL - 2.0*cosI/L^2*avrDLLz + 1.0/L^2*avrDLzLz

    return avrDJr, avrDL, avrDI, avrDJrJr, avrDLL, avrDII, avrDJrL, avrDJrI, avrDLI
end

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

function orbitAverageActionCoeffs_cosI_Opti_nbKseparate(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
                                nbw::Int64=nbw_default,
                                nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                                 nbu::Int64=nbu0, m_test::Float64=m_field)


    avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz = orbitAverageActionCoeffsOpti_nbKseparate(Jr,L,cosI,m_field,alpha,nbAvr,nbAvrTh,nbw,nbvarphi,nbphi,nbu,m_test)

    avrDI = -cosI/L*avrDL + 1.0/L*avrDLz + cosI/L^2*avrDLL- 1.0/L^2*avrDLLz
    avrDJrI = -cosI/L*avrDJrL + 1.0/L*avrDJrLz
    avrDLI = -cosI/L*avrDLL + 1.0/L*avrDLLz
    avrDII = cosI^2/L^2*avrDLL - 2.0*cosI/L^2*avrDLLz + 1.0/L^2*avrDLzLz

    return avrDJr, avrDL, avrDI, avrDJrJr, avrDLL, avrDII, avrDJrL, avrDJrI, avrDLI
end

function orbitAverageActionCoeffs_cosI_OptiExactOld(Jr::Float64, L::Float64, cosI::Float64, m_field::Float64,
                                alpha::Float64=alphaRot, nbAvr::Int64=nbAvr_default,
                                nbK::Int64=nbK_default, nbu::Int64=nbu0, m_test::Float64=m_field)


    avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz = orbitAverageActionCoeffsOptiExactOld(Jr,L,cosI,m_field,alpha,nbAvr,nbK,nbu,m_test)

    avrDI = -cosI/L*avrDL + 1.0/L*avrDLz + cosI/L^2*avrDLL- 1.0/L^2*avrDLLz
    avrDJrI = -cosI/L*avrDJrL + 1.0/L*avrDJrLz
    avrDLI = -cosI/L*avrDLL + 1.0/L*avrDLLz
    avrDII = cosI^2/L^2*avrDLL - 2.0*cosI/L^2*avrDLLz + 1.0/L^2*avrDLzLz

    return avrDJr, avrDL, avrDI, avrDJrJr, avrDLL, avrDII, avrDJrL, avrDJrI, avrDLI
end
