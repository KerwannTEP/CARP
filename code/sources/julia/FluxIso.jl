include("Main.jl")

# Compute analytical isotropic 3D-flux
# Compare Lz-integrated 3D-flux to 2D-flux
# Check relation F_Lz = Lz/L F_L


#############################################################
# Local velocity deflections
# And vr,vt gradients
#############################################################


function localVelChangeIsoGrad(r::Float64, vr::Float64, vt::Float64,
                        m_field::Float64,
                        nbK::Int64=nbK_default, m_test::Float64=m_field)


    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + v^2/2.0
    L = r*vt

    dvPar = 0.0
    dvParSq = 0.0
    dvPerSq = 0.0
    dvPar_vr = 0.0
    dvParSq_vr = 0.0
    dvPerSq_vr = 0.0
    dvPar_vt = 0.0
    dvParSq_vt = 0.0
    dvPerSq_vt = 0.0

    dvPar_r = 0.0
    dvParSq_r = 0.0
    dvPerSq_r = 0.0

    for ivarphi=1:nbK
        varphi = pi*(ivarphi-0.5)/nbK
        sinvarphi, cosvarphi = sincos(varphi)

        wmax = v*cosvarphi + sqrt(vSq*cosvarphi^2 - 2.0*E)

        dvPar_phi = 0.0
        dvParSq_phi = 0.0
        dvPerSq_phi = 0.0
        dvPar_vr_phi = 0.0
        dvParSq_vr_phi = 0.0
        dvPerSq_vr_phi = 0.0
        dvPar_vt_phi = 0.0
        dvParSq_vt_phi = 0.0
        dvPerSq_vt_phi = 0.0

        dvPar_r_phi = 0.0
        dvParSq_r_phi = 0.0
        dvPerSq_r_phi = 0.0

        for iphi=1:nbK
            phi = 2.0*pi*(iphi-0.5)/nbK
            sinphi, cosphi = sincos(phi)

            dvPar_w = 0.0
            dvParSq_w = 0.0
            dvPerSq_w = 0.0
            dvPar_vr_w = 0.0
            dvParSq_vr_w = 0.0
            dvPerSq_vr_w = 0.0
            dvPar_vt_w = 0.0
            dvParSq_vt_w = 0.0
            dvPerSq_vt_w = 0.0

            dvPar_r_w = 0.0
            dvParSq_r_w = 0.0
            dvPerSq_r_w = 0.0

            for iw=1:nbK
                w = wmax*(iw-0.5)/nbK

                w1 = w*cosvarphi
                w2 = w*sinvarphi*cosphi
                w3 = w*sinvarphi*sinphi

                Ep = E + w^2/2.0 - v*w1

                v1p = v-w1
                v2p = -w2
                v3p = -w3

                Lp = r * sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )

                # if (Lp < 10^(-3))
                #
                #     println(Lp)
                # end


                Ftot, dFtotdE, dFtotdL = _FdF(Ep,Lp)

                dEpdvr = vr - vr_v*w*cosvarphi
                dEpdvt = vt - vt_v*w*cosvarphi
                dLpdvr = r^2*(vt^2/v^3*w*sinvarphi*sinphi+vr*vt/v^3*w*cosvarphi)*(vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp
                dLpdvt = r^2*(1.0-vr*vt/v^3*w*sinvarphi*sinphi-vr^2/v^3*w*cosvarphi)*(vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp

                dEpdr = dpsidr(r)
                dLpdr = sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )

                dFdvr = dEpdvr*dFtotdE + dLpdvr*dFtotdL
                dFdvt = dEpdvt*dFtotdE + dLpdvt*dFtotdL

                dFdr = dEpdr*dFtotdE + dLpdr*dFtotdL

                # if (abs((vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp) >= 1.0)
                #
                #     println((vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp)
                # end

                dvPar_w += Ftot
                dvParSq_w += w*Ftot
                dvPerSq_w += w*Ftot
                dvPar_vr_w += dFdvr
                dvParSq_vr_w += w*dFdvr
                dvPerSq_vr_w += w*dFdvr
                dvPar_vt_w += dFdvt
                dvParSq_vt_w += w*dFdvt
                dvPerSq_vt_w += w*dFdvt

                dvPar_r_w += dFdr
                dvParSq_r_w += w*dFdr
                dvPerSq_r_w += w*dFdr

            end

            dvPar_w *= wmax
            dvParSq_w *= wmax
            dvPerSq_w *= wmax
            dvPar_vr_w *= wmax
            dvParSq_vr_w *= wmax
            dvPerSq_vr_w *= wmax
            dvPar_vt_w *= wmax
            dvParSq_vt_w *= wmax
            dvPerSq_vt_w *= wmax

            dvPar_r_w *= wmax
            dvParSq_r_w *= wmax
            dvPerSq_r_w *= wmax

            dvPar_phi += dvPar_w
            dvParSq_phi += dvParSq_w
            dvPerSq_phi += dvPerSq_w
            dvPar_vr_phi += dvPar_vr_w
            dvParSq_vr_phi += dvParSq_vr_w
            dvPerSq_vr_phi += dvPerSq_vr_w
            dvPar_vt_phi += dvPar_vt_w
            dvParSq_vt_phi += dvParSq_vt_w
            dvPerSq_vt_phi += dvPerSq_vt_w

            dvPar_r_phi += dvPar_r_w
            dvParSq_r_phi += dvParSq_r_w
            dvPerSq_r_phi += dvPerSq_r_w

        end

        dvPar += 2.0*sinvarphi*cosvarphi*dvPar_phi
        dvParSq += sinvarphi^3*dvParSq_phi
        dvPerSq += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_phi
        dvPar_vr += 2.0*sinvarphi*cosvarphi*dvPar_vr_phi
        dvParSq_vr += sinvarphi^3*dvParSq_vr_phi
        dvPerSq_vr += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_vr_phi
        dvPar_vt += 2.0*sinvarphi*cosvarphi*dvPar_vt_phi
        dvParSq_vt += sinvarphi^3*dvParSq_vt_phi
        dvPerSq_vt += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_vt_phi

        dvPar_r += 2.0*sinvarphi*cosvarphi*dvPar_r_phi
        dvParSq_r += sinvarphi^3*dvParSq_r_phi
        dvPerSq_r += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_r_phi

    end
    # A = 2 pi G^2 ln Lambda

    pref = 2.0*pi^2/nbK^3
    pref *= 2.0*pi*_G^2*logCoulomb

    dvPar *= -pref * (m_field + m_test)
    dvParSq *= 2.0*pref * m_field
    dvPerSq *= 2.0*pref * m_field
    dvPar_vr *= -pref * (m_field + m_test)
    dvParSq_vr *= 2.0*pref * m_field
    dvPerSq_vr *= 2.0*pref * m_field
    dvPar_vt *= -pref * (m_field + m_test)
    dvParSq_vt *= 2.0*pref * m_field
    dvPerSq_vt *= 2.0*pref * m_field

    dvPar_r *= -pref * (m_field + m_test)
    dvParSq_r *= 2.0*pref * m_field
    dvPerSq_r *= 2.0*pref * m_field

    #println(dvPerSq_vt)

    return dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt, dvPar_r, dvParSq_r, dvPerSq_r
end

#############################################################
# Local angle-averaged energy diffusion coeffivients
# E,L,Lz gradients
#############################################################



function dE_Iso(r::Float64, E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, m_test::Float64=m_field)

    vr = vr_from_E_L(r,E,L)
    vt = vt_from_E_L(r,E,L)
    v2 = vSq_from_E_L(r,E,L)
    v = sqrt(v2)
    vr_v = vr/v
    vt_v = vt/v

    dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt, dvPar_r, dvParSq_r, dvPerSq_r = localVelChangeIsoGrad(r,vr,vt,m_field,nbK,m_test)

    dvParSq_E = dvParSq_vr/vr
    dvPerSq_E = dvPerSq_vr/vr
    dvPar_E = dvPar_vr/vr

    dvParSq_L = (1.0/r*dvParSq_vt-L/(r^2*vr)*dvParSq_vr)
    dvPerSq_L = (1.0/r*dvPerSq_vt-L/(r^2*vr)*dvPerSq_vr)
    dvPar_L = (1.0/r*dvPar_vt-L/(r^2*vr)*dvPar_vr)




    # dE
    dE_E = 0.5*dvParSq_E + 0.5*dvPerSq_E + 1.0/v*dvPar + v*dvPar_E
    dE_L = 0.5*dvParSq_L + 0.5*dvPerSq_L + v*dvPar_L

    return dE_E, dE_L


end

function dELLz_Iso(r::Float64, E::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, m_test::Float64=m_field)

    vr = vr_from_E_L(r,E,L)
    vt = vt_from_E_L(r,E,L)
    v2 = vSq_from_E_L(r,E,L)
    v = sqrt(v2)
    vr_v = vr/v
    vt_v = vt/v

    # println("vr = ",vr)

    dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt, dvPar_r, dvParSq_r, dvPerSq_r = localVelChangeIsoGrad(r,vr,vt,m_field,nbK,m_test)

    # println("test = ",(dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt))

    dvParSq_E = dvParSq_vr/vr
    dvPerSq_E = dvPerSq_vr/vr
    dvPar_E = dvPar_vr/vr

    dvParSq_L = (1.0/r*dvParSq_vt-L/(r^2*vr)*dvParSq_vr)
    dvPerSq_L = (1.0/r*dvPerSq_vt-L/(r^2*vr)*dvPerSq_vr)
    dvPar_L = (1.0/r*dvPar_vt-L/(r^2*vr)*dvPar_vr)

    dvParSq_r = dvParSq_r - dpsiEffdr(r,L)/vr*dvParSq_vr - L/r^2*dvParSq_vt
    dvPerSq_r = dvPerSq_r - dpsiEffdr(r,L)/vr*dvPerSq_vr - L/r^2*dvPerSq_vt
    dvPar_r = dvPar_r - dpsiEffdr(r,L)/vr*dvPar_vr - L/r^2*dvPar_vt

    # println("dvParSq_vr = ",dvParSq_vr)
    # println("dvParSq_vr/vr = ",dvParSq_vr/vr)


    # dE
    dE_r = 0.5*dvParSq_r + 0.5*dvPerSq_r - dpsidr(r)/v*dvPar + v*dvPar_r
    dE_E = 0.5*dvParSq_E + 0.5*dvPerSq_E + 1.0/v*dvPar + v*dvPar_E
    dE_L = 0.5*dvParSq_L + 0.5*dvPerSq_L + v*dvPar_L

    # dEE
    dEE_r = -2.0*dpsidr(r)*dvParSq + v^2*dvParSq_r
    dEE_E = 2.0*dvParSq + v^2*dvParSq_E
    dEE_L = v2*dvParSq_L

    # dL
    dL_r = L*dpsidr(r)/v^3*dvPar + L/v*dvPar_r + 0.5*r/L*dvPerSq + 0.25*r^2/L*dvPerSq_r
    dL_E = -L/v^3*dvPar + L/v*dvPar_E + 0.25*r^2/L*dvPerSq_E
    dL_L = 1.0/v*dvPar + L/v*dvPar_L - 0.25*r^2/L^2*dvPerSq + 0.25*r^2/L*dvPerSq_L

    # dLL
    dLL_r = 2.0*L^2*dpsidr(r)/v^4*dvParSq + L^2/v^2*dvParSq_r + 0.5*(2.0*r*vr^2/v^2-2.0*r^2/v^2*dpsiEffdr(r,L)+2.0*r^2*vr^2*dpsidr(r)/v^4)*dvPerSq + 0.5*r^2*vr^2/v^2*dvPerSq_r
    dLL_E = -2.0*L^2/v^4*dvParSq + L^2/v^2*dvParSq_E + L^2/v^4*dvPerSq + 0.5*r^2*vr^2/v^2*dvPerSq_E
    dLL_L = 2.0*L/v^2*dvParSq + L^2/v^2*dvParSq_L - L/v^2*dvPerSq + 0.5*r^2*vr^2/v^2*dvPerSq_L

    # dEL
    dEL_r = L*dvParSq_r
    dEL_E = L*dvParSq_E
    dEL_L = dvParSq + L*dvParSq_L

    # dLz
    dLz_r = Lz*dpsidr(r)/v^3*dvPar + Lz/v*dvPar_r
    dLz_E = -Lz/v^3*dvPar + Lz/v*dvPar_E
    dLz_L = Lz/v*dvPar_L
    dLz_Lz = 1.0/v*dvPar

    # dLzLz
    dLzLz_r = (Lz^2*dpsidr(r)/v^4*dvParSq + Lz^2/v^2*dvParSq_r + 0.5*Lz^2/L^2*(2.0*r*vr^2/v^2-2.0*r^2/v^2*dpsiEffdr(r,L)+2.0*r^2*vr^2*dpsidr(r)/v^4)*dvPerSq
            + Lz^2/L^2*0.5*r^2*vr^2/v^2*dvPerSq_r + 0.5*r*(1.0-Lz^2/L^2)*dvPerSq + 0.25*r^2*(1.0-Lz^2/L^2)*dvPerSq_r)
    dLzLz_E = (-2.0*Lz^2/v^4*dvParSq + Lz^2/v^2*dvParSq_E + Lz^2/v^4*dvPerSq+0.5*Lz^2/L^2*r^2*vr^2/v^2*dvPerSq_E
            + 0.25*r^2*(1.0-Lz^2/L^2)*dvPerSq_E)
    dLzLz_L = (Lz^2/v^2*dvParSq_L - r^2*Lz^2/L^3*dvPerSq + 0.5*Lz^2/L^2*r^2*vr^2/v^2*dvPerSq_L
            + 0.5*r^2*Lz^2/L^3*dvPerSq + 0.25*r^2*(1.0-Lz^2/L^2)*dvPerSq_L)
    # dLzLz_L = (Lz^2/v^2*dvParSq_L - Lz^2*r^2*vr^2/(v^2*L^3)*dvPerSq + Lz^2/L^2*r^2*vr^2/(2.0*v^2)*dvPerSq_L
    #         + 0.5*r^2*Lz^2/L^3**dvPerSq + 0.25*r^2*(1.0-Lz^2/L^2)*dvPerSq_L)
    dLzLz_Lz = 2.0*Lz/v^2*dvParSq + Lz/L^2*r^2*vr^2/v^2*dvPerSq - 0.5*r^2*Lz/L^2*dvPerSq

    # dELz
    dELz_r = Lz*dvParSq_r
    dELz_E = Lz*dvParSq_E
    dELz_L = Lz*dvParSq_L
    dELz_Lz = dvParSq

    # dLLz
    #dLLz_r = 2.0*Lz*L*dpsidr(r)/v^4*dvParSq + Lz*L/v^2*dvParSq_r + Lz/L*0.5*(2.0*r+2.0*L^2*dpsidr(r)/v^4)*dvPerSq + Lz/L*0.5*r^2*vr^2/v^2*dvPerSq_r
    dLLz_r = 2.0*Lz*L*dpsidr(r)/v^4*dvParSq + Lz*L/v^2*dvParSq_r + Lz/L*0.5*(2.0*r*vr^2/v^2-2.0*r^2/v^2*dpsiEffdr(r,L)+2.0*r^2*vr^2*dpsidr(r)/v^4)*dvPerSq + Lz/L*0.5*r^2*vr^2/v^2*dvPerSq_r
    dLLz_E = -2.0*Lz*L/v^4*dvParSq + Lz*L/v^2*dvParSq_E + Lz*L/v^4*dvPerSq + 0.5*Lz/L*r^2*vr^2/v^2*dvPerSq_E
    dLLz_L = Lz/v^2*dvParSq + Lz*L/v^2*dvParSq_L - 0.5*r^2*Lz/L^2*dvPerSq - 0.5*Lz/v^2*dvPerSq + 0.5*Lz/L*r^2*vr^2/v^2*dvPerSq_L
    dLLz_Lz = L/v^2*dvParSq + 0.5*r^2/L*vr^2/v^2*dvPerSq




    return dE_E, dE_L, dL_E, dL_L, dLz_E, dLz_L, dLz_Lz, dEE_E, dEE_L, dLL_E, dLL_L, dLzLz_E, dLzLz_L, dLzLz_Lz, dEL_E, dEL_L, dELz_E, dELz_L, dELz_Lz, dLLz_E, dLLz_L, dLLz_Lz, dE_r, dL_r, dLz_r,  dEE_r, dLL_r, dLzLz_r, dEL_r, dELz_r, dLLz_r
end

#############################################################
# Orbit-average
#############################################################

# compute E,L,Lz-derivatives of DEE, etc...
function DELLz_Iso_OrbitAvg(E::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field)

    cosI = Lz/L



    dLLz_E = 0.0
    dLLz_L = 0.0
    dLLz_Lz = 0.0

    halfperiod = 0.0
    halfperiod_E = 0.0
    halfperiod_L = 0.0

    sp, sa = radius_s_bounds(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)

    dE_1_E = 0.0
    dE_1_L = 0.0
    dE_2 = 0.0
    dL_1_E = 0.0
    dL_1_L = 0.0
    dL_2 = 0.0
    dLz_1_E = 0.0
    dLz_1_L = 0.0
    dLz_1_Lz = 0.0
    dLz_2 = 0.0

    dEE_1_E = 0.0
    dEE_1_L = 0.0
    dEE_2 = 0.0
    dLL_1_E = 0.0
    dLL_1_L = 0.0
    dLL_2 = 0.0
    dLzLz_1_E = 0.0
    dLzLz_1_L = 0.0
    dLzLz_1_Lz = 0.0
    dLzLz_2 = 0.0

    dEL_1_E = 0.0
    dEL_1_L = 0.0
    dEL_2 = 0.0
    dELz_1_E = 0.0
    dELz_1_L = 0.0
    dELz_1_Lz = 0.0
    dELz_2 = 0.0

    dLLz_1_E = 0.0
    dLLz_1_L = 0.0
    dLLz_1_Lz = 0.0
    dLLz_2 = 0.0


    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        sloc = s_from_u_sma_ecc(uloc,sma,ecc)
        rloc = r_from_s(sloc)
        jac_loc = Theta(uloc,sp,sa)

        djacdE, djacdL, dsdL, dsdE = djac_and_ds(uloc,sp,sa)

        vr = sqrt(2.0*(E - psiEff(rloc,L)))
        vt = L/rloc



        drudE = _b^2*sloc/rloc * dsdE
        drudL = _b^2*sloc/rloc * dsdL



        dE, dL, dLz, dEE, dLL, dLzLz, dEL, dELz, dLLz = localOrbitChangeNoRot(rloc,vr,vt,cosI,m_field,nbK,m_test)

        dE_E, dE_L, dL_E, dL_L, dLz_E, dLz_L, dLz_Lz, dEE_E, dEE_L, dLL_E, dLL_L, dLzLz_E, dLzLz_L, dLzLz_Lz, dEL_E, dEL_L, dELz_E, dELz_L, dELz_Lz, dLLz_E, dLLz_L, dLLz_Lz, dE_r, dL_r, dLz_r,  dEE_r, dLL_r, dLzLz_r, dEL_r, dELz_r, dLLz_r = dELLz_Iso(rloc,E,L,Lz,m_field,nbK,m_test)

        # T/2

        halfperiod += jac_loc
        halfperiod_E += djacdE
        halfperiod_L += djacdL

        # dE
        dE_1_E += djacdE*dE + jac_loc*dE_E + jac_loc*drudE*dE_r
        dE_1_L += djacdL*dE + jac_loc*dE_L+ jac_loc*drudL*dE_r
        dE_2 += jac_loc*dE

        # dL
        dL_1_E += djacdE*dL + jac_loc*dL_E + jac_loc*drudE*dL_r
        dL_1_L += djacdL*dL + jac_loc*dL_L + jac_loc*drudL*dL_r
        dL_2 += jac_loc*dL

        # dLz
        dLz_1_E += djacdE*dLz + jac_loc*dLz_E + jac_loc*drudE*dLz_r
        dLz_1_L += djacdL*dLz + jac_loc*dLz_L + jac_loc*drudL*dLz_r
        dLz_1_Lz += jac_loc*dLz_Lz
        dLz_2 += jac_loc*dLz


        # dEE
        dEE_1_E += djacdE*dEE + jac_loc*dEE_E + jac_loc*drudE*dEE_r
        dEE_1_L += djacdL*dEE + jac_loc*dEE_L + jac_loc*drudL*dEE_r
        dEE_2 += jac_loc*dEE

        # dLL
        dLL_1_E += djacdE*dLL + jac_loc*dLL_E + jac_loc*drudE*dLL_r
        dLL_1_L += djacdL*dLL + jac_loc*dLL_L + jac_loc*drudL*dLL_r
        dLL_2 += jac_loc*dLL

        # dLzLz
        dLzLz_1_E += djacdE*dLzLz + jac_loc*dLzLz_E + jac_loc*drudE*dLzLz_r
        dLzLz_1_L += djacdL*dLzLz + jac_loc*dLzLz_L + jac_loc*drudL*dLzLz_r
        dLzLz_1_Lz += jac_loc*dLzLz_Lz
        dLzLz_2 += jac_loc*dLzLz

        # dEL
        dEL_1_E += djacdE*dEL + jac_loc*dEL_E + jac_loc*drudE*dEL_r
        dEL_1_L += djacdL*dEL + jac_loc*dEL_L + jac_loc*drudL*dEL_r
        dEL_2 += jac_loc*dEL

        # dELz
        dELz_1_E += djacdE*dELz + jac_loc*dELz_E + jac_loc*drudE*dELz_r
        dELz_1_L += djacdL*dELz + jac_loc*dELz_L + jac_loc*drudL*dELz_r
        dELz_1_Lz += jac_loc*dELz_Lz
        dELz_2 += jac_loc*dELz

        # dLLz
        dLLz_1_E += djacdE*dLLz + jac_loc*dLLz_E + jac_loc*drudE*dLLz_r
        dLLz_1_L += djacdL*dLLz + jac_loc*dLLz_L + jac_loc*drudL*dLLz_r
        dLLz_1_Lz += jac_loc*dLLz_Lz
        dLLz_2 += jac_loc*dLLz

    end

    dE_E = (halfperiod*dE_1_E - halfperiod_E*dE_2)/halfperiod^2
    dE_L = (halfperiod*dE_1_L - halfperiod_L*dE_2)/halfperiod^2
    dL_E = (halfperiod*dL_1_E - halfperiod_E*dL_2)/halfperiod^2
    dL_L = (halfperiod*dL_1_L - halfperiod_L*dL_2)/halfperiod^2
    dLz_E = (halfperiod*dLz_1_E - halfperiod_E*dLz_2)/halfperiod^2
    dLz_L = (halfperiod*dLz_1_L - halfperiod_L*dLz_2)/halfperiod^2
    dLz_Lz = dLz_1_Lz/halfperiod

    dEE_E = (halfperiod*dEE_1_E - halfperiod_E*dEE_2)/halfperiod^2
    dEE_L = (halfperiod*dEE_1_L - halfperiod_L*dEE_2)/halfperiod^2
    dLL_E = (halfperiod*dLL_1_E - halfperiod_E*dLL_2)/halfperiod^2
    dLL_L = (halfperiod*dLL_1_L - halfperiod_L*dLL_2)/halfperiod^2
    dLzLz_E = (halfperiod*dLzLz_1_E - halfperiod_E*dLzLz_2)/halfperiod^2
    dLzLz_L = (halfperiod*dLzLz_1_L - halfperiod_L*dLzLz_2)/halfperiod^2
    dLzLz_Lz = dLzLz_1_Lz/halfperiod

    dEL_E = (halfperiod*dEL_1_E - halfperiod_E*dEL_2)/halfperiod^2
    dEL_L = (halfperiod*dEL_1_L - halfperiod_L*dEL_2)/halfperiod^2
    dELz_E = (halfperiod*dELz_1_E - halfperiod_E*dELz_2)/halfperiod^2
    dELz_L = (halfperiod*dELz_1_L - halfperiod_L*dELz_2)/halfperiod^2
    dELz_Lz = dELz_1_Lz/halfperiod
    dLLz_E = (halfperiod*dLLz_1_E - halfperiod_E*dLLz_2)/halfperiod^2
    dLLz_L = (halfperiod*dLLz_1_L - halfperiod_L*dLLz_2)/halfperiod^2
    dLLz_Lz = dLLz_1_Lz/halfperiod

    return dE_E, dE_L, dL_E, dL_L, dLz_E, dLz_L, dLz_Lz, dEE_E, dEE_L, dLL_E, dLL_L, dLzLz_E, dLzLz_L, dLzLz_Lz, dEL_E, dEL_L, dELz_E, dELz_L, dELz_Lz, dLLz_E, dLLz_L, dLLz_Lz

end




# compute E,L,Lz-derivatives of DEE, etc...
# Second order diffusion coefficients gradients
function DJrLLz_Iso_OrbitAvg(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field)

    E = E_from_Jr_L(Jr,L,nbu)
    cosI = Lz/L

    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L_Wrap(E,L,nbu)
    dE_E, dE_L, dL_E, dL_L, dLz_E, dLz_L, dLz_Lz, dEE_E, dEE_L, dLL_E, dLL_L, dLzLz_E, dLzLz_L, dLzLz_Lz, dEL_E, dEL_L, dELz_E, dELz_L, dELz_Lz, dLLz_E, dLLz_L, dLLz_Lz = DELLz_Iso_OrbitAvg(E,L,Lz,m_field,nbK,nbAvr,m_test)

    dE, dL, dLz, dEE, dLL, dLzLz, dEL, dELz, dLLz = orbitAverageEnergyCoeffsNoRot(E,L,cosI,m_field,nbAvr,nbK,m_test)

    dEdJr = 1.0/dJrdE
    dEdL = -dJrdL/dJrdE


    # JrJr
    dJrJr_E = (2.0*d2JrdE2*dJrdE*dEE + dJrdE^2*dEE_E + 2.0*d2JrdE2*dJrdL*dEL
            + 2.0*dJrdE*d2JrdEL*dEL + 2.0*dJrdE*dJrdL*dEL_E + 2.0*d2JrdEL*dJrdL*dLL
            + dJrdL^2*dLL_E)

    dJrJr_L = (2.0*d2JrdEL*dJrdE*dEE + dJrdE^2*dEE_L + 2.0*d2JrdEL*dJrdL*dEL
            + 2.0*dJrdE*d2JrdL2*dEL + 2.0*dJrdE*dJrdL*dEL_L + 2.0*d2JrdL2*dJrdL*dLL
            + dJrdL^2*dLL_L)

    # JrL
    dJrL_E = d2JrdE2*dEL + dJrdE*dEL_E + d2JrdEL*dLL + dJrdL*dLL_E
    dJrL_L = d2JrdEL*dEL + dJrdE*dEL_L + d2JrdL2*dLL + dJrdL*dLL_L

    # JrLz
    dJrLz_E = d2JrdE2*dELz + dJrdE*dELz_E + d2JrdEL*dLLz + dJrdL*dLLz_E
    dJrLz_L = d2JrdEL*dELz + dJrdE*dELz_L + d2JrdL2*dLLz + dJrdL*dLLz_L
    dJrLz_Lz = dJrdE*dELz_Lz + dJrdL*dLLz_Lz



    # Obtain Jr, L derivatives

    dJrJr_Jr = dEdJr*dJrJr_E
    dJrJr_L = dEdL*dJrJr_E + dJrJr_L
    dLL_Jr = dEdJr*dLL_E
    dLL_L = dEdL*dLL_E + dLL_L
    dLzLz_Jr = dEdJr*dLzLz_E
    dLzLz_L = dEdL*dLzLz_E + dLzLz_L

    dJrL_Jr = dEdJr*dJrL_E
    dJrL_L = dEdL*dJrL_E + dJrL_L
    dJrLz_Jr = dEdJr*dJrLz_E
    dJrLz_L = dEdL*dJrLz_E + dJrLz_L
    dLLz_Jr = dEdJr*dLLz_E
    dLLz_L = dEdL*dLLz_E + dLLz_L

    return dJrJr_Jr, dJrJr_L, dLL_Jr, dLL_L, dLzLz_Jr, dLzLz_L, dLzLz_Lz, dJrL_Jr, dJrL_L, dJrLz_Jr, dJrLz_L, dJrLz_Lz, dLLz_Jr, dLLz_L, dLLz_Lz
end

#############################################################
# 3D flux without rotstion
#############################################################

function flux_Iso(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field)

    cosI = Lz/L
    E = E_from_Jr_L(Jr,L,nbu)

    dJrJr_Jr, dJrJr_L, dLL_Jr, dLL_L, dLzLz_Jr, dLzLz_L, dLzLz_Lz, dJrL_Jr, dJrL_L, dJrLz_Jr, dJrLz_L, dJrLz_Lz, dLLz_Jr, dLLz_L, dLLz_Lz = DJrLLz_Iso_OrbitAvg(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test)
    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsNoRot(Jr,L,cosI,m_field,nbAvr,nbK,nbu,m_test)

    Ftot, dFdE, dFdL = _FdF(E,L)

    dJrdE, dJrdL, _ = grad_Jr_E_L_Wrap(E,L,nbu)
    dEdJr = 1.0/dJrdE
    dEdL = -dJrdL/dJrdE

    dFdJr = dEdJr*dFdE
    dFdL = dEdL*dFdE + dFdL

    # Jr
    fluxJr = (dJr-0.5*dJrJr_Jr-0.5*dJrL_L-0.5*dJrLz_Lz)*Ftot - 0.5*dJrJr*dFdJr - 0.5*dJrL*dFdL

    # L
    fluxL = (dL-0.5*dJrL_Jr-0.5*dLL_L-0.5*dLLz_Lz)*Ftot - 0.5*dJrL*dFdJr - 0.5*dLL*dFdL

    # Lz
    fluxLz = (dLz-0.5*dJrLz_Jr-0.5*dLLz_L-0.5*dLzLz_Lz)*Ftot - 0.5*dJrLz*dFdJr - 0.5*dLLz*dFdL

    return fluxJr, fluxL, fluxLz
end



function dFdt_NoRot(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    # Jr
    fluxJr_p, _, _ = flux_Iso(Jr_p,L,Lz,m_field,nbK,nbAvr,nbu,m_test)
    fluxJr_m, _, _ = flux_Iso(Jr_m,L,Lz,m_field,nbK,nbAvr,nbu,m_test)

    dJr = (fluxJr_p - fluxJr_m)/(2.0*eps*_L0)

    # L
    _, fluxL_p, _ = flux_Iso(Jr,L_p,Lz,m_field,nbK,nbAvr,nbu,m_test)
    _, fluxL_m, _ = flux_Iso(Jr,L_m,Lz,m_field,nbK,nbAvr,nbu,m_test)

    dL = (fluxL_p - fluxL_m)/(2.0*eps*_L0)

    # Lz
    _, _, fluxLz_p = flux_Iso(Jr,L,Lz_p,m_field,nbK,nbAvr,nbu,m_test)
    _, _, fluxLz_m = flux_Iso(Jr,L,Lz_m,m_field,nbK,nbAvr,nbu,m_test)

    dLz = (fluxLz_p - fluxLz_m)/(2.0*eps*_L0)

    # println("Jr = ",dJr)
    # println("L  = ",dL)
    # println("Lz = ",dLz)

    return -(dJr + dL + dLz)
end


function integrateLz_dFdt_NoRotPar(Jr::Float64, L::Float64, m_field::Float64, nbLz::Int64=12,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    sum = Threads.Atomic{Float64}(0.0)

    Threads.@threads for iLz=1:nbLz
        Lz = -L + 2.0*L/nbLz * (iLz-0.5)
        dFdt = dFdt_NoRot(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test,eps)

        Threads.atomic_add!(sum,dFdt)

    end

    sum[] *=  2.0*L/nbLz

    return sum[]

end

function integrateLz_dFdt_NoRot(Jr::Float64, L::Float64, m_field::Float64, nbLz::Int64=12,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    sum = 0.0

    for iLz=1:nbLz
        Lz = -L + 2.0*L/nbLz * (iLz-0.5)
        dFdt = dFdt_NoRot(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test,eps)

        sum += dFdt

    end

    sum *=  2.0*L/nbLz

    return sum

end


function integrateJr_dFdt_NoRotPar(L::Float64, Lz::Float64, m_field::Float64,
            Jrmax::Float64=1.0*_L0, nbJr::Int64=12,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    sum = Threads.Atomic{Float64}(0.0)

    Threads.@threads for iJr=1:nbJr
        Jr = Jrmax/nbJr * (iJr-0.5)
        dFdt = dFdt_NoRot(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test,eps)

        Threads.atomic_add!(sum,dFdt)

    end

    sum[] *=  Jrmax/nbJr

    return sum[]

end

function integrateJrSq_dFdt_NoRotPar(L::Float64, Lz::Float64, m_field::Float64,
            Jrmax::Float64=1.0*_L0, nbJr::Int64=12,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    sum = Threads.Atomic{Float64}(0.0)
    xmax = sqrt(Jrmax/_L0)

    Threads.@threads for ix=1:nbJr
        x = xmax * (ix-0.5)/nbJr
        Jr = _L0 * x^2
        dFdt = dFdt_NoRot(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test,eps)

        Threads.atomic_add!(sum,x*dFdt)

    end

    sum[] *=  xmax/nbJr * 2.0*_L0

    return sum[]

end

function integrateJrCbrt_dFdt_NoRotPar(L::Float64, Lz::Float64, m_field::Float64,
            Jrmax::Float64=1.0*_L0, nbJr::Int64=12,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    sum = Threads.Atomic{Float64}(0.0)
    xmax = cbrt(Jrmax/_L0)

    Threads.@threads for iJr=1:nbJr
        x = xmax * (ix-0.5)/nbJr
        Jr = _L0 * x^3
        dFdt = dFdt_NoRot(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test,eps)

        Threads.atomic_add!(sum,x^2*dFdt)

    end

    sum[] *=  xmax/nbJr * 3.0*_L0

    return sum[]

end

function integrateJrPow_dFdt_NoRotPar(L::Float64, Lz::Float64, m_field::Float64, pow::Int64=3,
            Jrmax::Float64=3.0*_L0, nbJr::Int64=30,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    sum = Threads.Atomic{Float64}(0.0)
    xmax = (Jrmax/_L0)^(1/pow)

    Threads.@threads for ix=1:nbJr
        x = xmax * (ix-0.5)/nbJr
        Jr = _L0 * x^pow
        dFdt = dFdt_NoRot(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test,eps)

        Threads.atomic_add!(sum,x^(pow-1)*dFdt)

    end

    sum[] *=  xmax/nbJr * pow*_L0

    return sum[]

end

function integrateJrPow_dFdt_NoRot(L::Float64, Lz::Float64, m_field::Float64, pow::Int64=3,
            Jrmax::Float64=3.0*_L0, nbJr::Int64=30,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    sum = 0.0
    xmax = (Jrmax/_L0)^(1/pow)

    for ix=1:nbJr
        x = xmax * (ix-0.5)/nbJr
        Jr = _L0 * x^pow
        dFdt = dFdt_NoRot(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test,eps)

        sum +=x^(pow-1)*dFdt

    end

    sum *=  xmax/nbJr * pow*_L0

    return sum[]

end


# Compute dFdt(cosI)
# F(cosI) = (2pi)^3 int dJr dL L F(Jr,L,LcosI)
function integrateJrLPow_dFdt_NoRotPar(cosI::Float64, m_field::Float64, pow::Int64=3,
            Jrmax::Float64=3.0*_L0, nbJr::Int64=30,
            Lmax::Float64=2.0*_L0, nbL::Int64=20,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

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


        dfdt = dFdt_NoRot(Jr,L,L*cosI,m_field,nbK,nbAvr,nbu,m_test,epsEff)



        Threads.atomic_add!(sum,x^(pow-1)*y^(pow-1)*L*dfdt)

    end

    sum[] *=  (2*pi)^3*xmax/nbJr * pow*_L0 * ymax/nbL * pow*_L0

    return sum[]

end

#############################################################
# 2D flux
#############################################################

function flux2D_Iso(Jr::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field)


    E = E_from_Jr_L(Jr,L,nbu)

    dJrJr_Jr, dJrJr_L, dLL_Jr, dLL_L, _, _, _, dJrL_Jr, dJrL_L, _, _, _, _, _, _ = DJrLLz_Iso_OrbitAvg(Jr,L,L,m_field,nbK,nbAvr,nbu,m_test)
    dJr, dL, dLz, dJrJr, dLL, _, dJrL, _, _ = orbitAverageActionCoeffsNoRot(Jr,L,1.0,m_field,nbAvr,nbK,nbu,m_test)

    Ftot, dFdE, dFdL = _FdF(E,L)

    dJrdE, dJrdL, _ = grad_Jr_E_L_Wrap(E,L,nbu)
    dEdJr = 1.0/dJrdE
    dEdL = -dJrdL/dJrdE

    dFdJr = dEdJr*dFdE*2.0*L
    dFdL = dEdL*dFdE*2.0*L + 2.0*L*dFdL + 2.0*Ftot

    # Jr
    fluxJr = (dJr-0.5*dJrJr_Jr-0.5*dJrL_L)*2.0*L*Ftot - 0.5*dJrJr*dFdJr - 0.5*dJrL*dFdL

    # L
    fluxL = (dL-0.5*dJrL_Jr-0.5*dLL_L)*2.0*L*Ftot - 0.5*dJrL*dFdJr - 0.5*dLL*dFdL


    return fluxJr, fluxL
end

function dFdt_NoRot2D(Jr::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0


    # Jr
    fluxJr_p, _ = flux2D_Iso(Jr_p,L,m_field,nbK,nbAvr,nbu,m_test)
    fluxJr_m, _ = flux2D_Iso(Jr_m,L,m_field,nbK,nbAvr,nbu,m_test)

    dJr = (fluxJr_p - fluxJr_m)/(2.0*eps*_L0)

    # println("dJr:",dJr)

    # L
    _, fluxL_p = flux2D_Iso(Jr,L_p,m_field,nbK,nbAvr,nbu,m_test)
    _, fluxL_m = flux2D_Iso(Jr,L_m,m_field,nbK,nbAvr,nbu,m_test)

    dL = (fluxL_p - fluxL_m)/(2.0*eps*_L0)

    # println("dL:",dL)


    return -(dJr + dL)
end
