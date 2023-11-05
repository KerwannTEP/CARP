include("FluxIso.jl")

#############################################################
# Check local velocity deflection gradients
#############################################################

function localVelChangeIsoNum(r::Float64, vr::Float64, vt::Float64,
                        m_field::Float64,
                        nbK::Int64=nbK_default, m_test::Float64=m_field, eps::Float64=10^(-3))

    vr_p = vr + _v0*eps
    vr_m = vr - _v0*eps
    vt_p = vt + _v0*eps
    vt_m = vt - _v0*eps

    r_p = r + eps*_b
    r_m = r - eps*_b



    dvPar_rp, dvParSq_rp, dvPerSq_rp, _ = localVelChangeIso(r,vr_p,vt,m_field,nbK,m_test)
    dvPar_rm, dvParSq_rm, dvPerSq_rm, _ = localVelChangeIso(r,vr_m,vt,m_field,nbK,m_test)
    dvPar_tp, dvParSq_tp, dvPerSq_tp, _ = localVelChangeIso(r,vr,vt_p,m_field,nbK,m_test)
    dvPar_tm, dvParSq_tm, dvPerSq_tm, _ = localVelChangeIso(r,vr,vt_m,m_field,nbK,m_test)

    dvPar_rrp, dvParSq_rrp, dvPerSq_rrp, _ = localVelChangeIso(r_p,vr,vt,m_field,nbK,m_test)
    dvPar_rrm, dvParSq_rrm, dvPerSq_rrm, _ = localVelChangeIso(r_m,vr,vt,m_field,nbK,m_test)


    # dvPar
    dvPar_vr = (dvPar_rp-dvPar_rm)/(2.0*eps*_v0)
    dvPar_vt = (dvPar_tp-dvPar_tm)/(2.0*eps*_v0)
    dvPar_r = (dvPar_rrp-dvPar_rrm)/(2.0*eps*_b)

    # dvParSq
    dvParSq_vr = (dvParSq_rp-dvParSq_rm)/(2.0*eps*_v0)
    dvParSq_vt = (dvParSq_tp-dvParSq_tm)/(2.0*eps*_v0)
    dvParSq_r = (dvParSq_rrp-dvParSq_rrm)/(2.0*eps*_b)

    # dvPerSq
    dvPerSq_vr = (dvPerSq_rp-dvPerSq_rm)/(2.0*eps*_v0)
    dvPerSq_vt = (dvPerSq_tp-dvPerSq_tm)/(2.0*eps*_v0)
    dvPerSq_r = (dvPerSq_rrp-dvPerSq_rrm)/(2.0*eps*_b)

    return dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt, dvPar_r, dvParSq_r, dvPerSq_r
end


# TODO
# Returns dvPar, dvPar2, dvPerp2 and vr,vt gradients of those
function localVelChangeIso(r::Float64, vr::Float64, vt::Float64,
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

                dFdvr = dEpdvr*dFtotdE + dLpdvr*dFtotdL
                dFdvt = dEpdvt*dFtotdE + dLpdvt*dFtotdL

                # if (abs((vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp) >= 1.0)
                #
                #     println((vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp)
                # end

                vxp = -vt/v*(v-w*cosvarphi)-vr/v*w*sinvarphi*sinphi
                vzp = -vr/v*(v-w*cosvarphi)+vt/v*w*sinvarphi*sinphi

                dvPar_w += Ftot
                dvParSq_w += w*Ftot
                dvPerSq_w += w*Ftot
                dvPar_vr_w += dFdvr
                dvParSq_vr_w += (2.0*sinvarphi*sinphi*vt^2*vxp/v^3 - 3.0*sinvarphi^2*(cosvarphi*(vr/v-vt^2*vzp/v^3)+sinvarphi*sinphi*vt^2*vxp/v^3))*Ftot
                dvPerSq_vr_w += (4.0*cosvarphi*(vr/v-vt^2*vzp/v^3)+2.0*sinvarphi*sinphi*vt^2*vxp/v^3-3.0*(1.0+cosvarphi^2)*(cosvarphi*(vr/v-vt^2*vzp/v^3)+sinvarphi*sinphi*vt^2*vxp/v^3))*Ftot
                dvPar_vt_w += dFdvt
                dvParSq_vt_w += (-2.0*sinvarphi*sinphi*vr^2*vzp/v^3 - 3.0*sinvarphi^2*(cosvarphi*(vt/v-vr^2*vxp/v^3)-sinvarphi*sinphi*vr^2*vzp/v^3))*Ftot
                dvPerSq_vt_w += (4.0*cosvarphi*(vt/v-vr^2*vxp/v^3)-2.0*sinvarphi*sinphi*vr^2*vzp/v^3-3.0*(1.0+cosvarphi^2)*(cosvarphi*(vt/v-vr^2*vxp/v^3)-sinvarphi*sinphi*vr^2*vzp/v^3))*Ftot

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



            dvPar_phi += dvPar_w
            dvParSq_phi += dvParSq_w
            dvPerSq_phi += dvPerSq_w
            dvPar_vr_phi += dvPar_vr_w
            dvParSq_vr_phi += dvParSq_vr_w
            dvPerSq_vr_phi += dvPerSq_vr_w
            dvPar_vt_phi += dvPar_vt_w
            dvParSq_vt_phi += dvParSq_vt_w
            dvPerSq_vt_phi += dvPerSq_vt_w

        end

        dvPar += 2.0*sinvarphi*cosvarphi*dvPar_phi
        dvParSq += sinvarphi^3*dvParSq_phi
        dvPerSq += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_phi
        dvPar_vr += 2.0*sinvarphi*cosvarphi*dvPar_vr_phi
        dvParSq_vr += sinvarphi*dvParSq_vr_phi
        dvPerSq_vr += sinvarphi*dvPerSq_vr_phi
        dvPar_vt += 2.0*sinvarphi*cosvarphi*dvPar_vt_phi
        dvParSq_vt += sinvarphi*dvParSq_vt_phi
        dvPerSq_vt += sinvarphi*dvPerSq_vt_phi

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

    #println(dvPerSq_vt)

    return dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt
end


function localVelChangeIsoCompare(r::Float64, vr::Float64, vt::Float64,
                        m_field::Float64,
                        nbK::Int64=nbK_default, m_test::Float64=m_field, eps::Float64=10^(-5))

    _, _, _, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt = localVelChangeIso(r,vr,vt,m_field,nbK,m_test)
    _, _, _, dvPar_vrGrad, dvParSq_vrGrad, dvPerSq_vrGrad, dvPar_vtGrad, dvParSq_vtGrad, dvPerSq_vtGrad, dvPar_rGrad, dvParSq_rGrad, dvPerSq_rGrad = localVelChangeIsoGrad(r,vr,vt,m_field,nbK,m_test)
    dvPar_vr_num, dvParSq_vr_num, dvPerSq_vr_num, dvPar_vt_num, dvParSq_vt_num, dvPerSq_vt_num, dvPar_r_num, dvParSq_r_num, dvPerSq_r_num = localVelChangeIsoNum(r,vr,vt,m_field,nbK,m_test,eps)

    println("dvPar_vr : ")
    println("Analytical = ",dvPar_vr)
    println("Ana. Grad. = ",dvPar_vrGrad)
    println("Numerical  = ",dvPar_vr_num)
    println("-------------------")
    println("dvPar_vt : ")
    println("Analytical = ",dvPar_vt)
    println("Ana. Grad. = ",dvPar_vtGrad)
    println("Numerical  = ",dvPar_vt_num)
    println("-------------------")
    println("dvPar_r : ")
    println("Ana. Grad. = ",dvPar_rGrad)
    println("Numerical  = ",dvPar_r_num)
    println("-------------------")
    println("dvParSq_vr : ")
    println("Analytical = ",dvParSq_vr)
    println("Ana. Grad. = ",dvParSq_vrGrad)
    println("Numerical  = ",dvParSq_vr_num)
    println("-------------------")
    println("dvParSq_vt : ")
    println("Analytical = ",dvParSq_vt)
    println("Ana. Grad. = ",dvParSq_vtGrad)
    println("Numerical  = ",dvParSq_vt_num)
    println("-------------------")
    println("dvParSq_r : ")
    println("Ana. Grad. = ",dvParSq_rGrad)
    println("Numerical  = ",dvParSq_r_num)
    println("-------------------")
    println("dvPerpSq_vr : ")
    println("Analytical = ",dvPerSq_vr)
    println("Ana. Grad. = ",dvPerSq_vrGrad)
    println("Numerical  = ",dvPerSq_vr_num)
    println("-------------------")
    println("dvPerpSq_vt : ")
    println("Analytical = ",dvPerSq_vt)
    println("Ana. Grad. = ",dvPerSq_vtGrad)
    println("Numerical  = ",dvPerSq_vt_num)
    println("-------------------")
    println("dvPerpSq_r : ")
    println("Ana. Grad. = ",dvPerSq_rGrad)
    println("Numerical  = ",dvPerSq_r_num)
    println("-------------------")

end

#############################################################
# Check local energy diffusion gradients
#############################################################

function LOverv(r::Float64,E::Float64,L::Float64,m_field::Float64)
    v2 = vSq_from_E_L(r,E,L)
    vr = vr_from_E_L(r,E,L)
    vt = vt_from_E_L(r,E,L)
    v = sqrt(v2)
    dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt, dvPar_r, dvParSq_r, dvPerSq_r = localVelChangeIsoGrad(r,vr,vt,m_field)

    return L/v*dvPar
end

function LOverv_r_num(r::Float64,E::Float64,L::Float64,m_field::Float64,eps::Float64=10^(-5))

    r_p = r + eps*_b
    r_m = r - eps*_b

    v2_p = vSq_from_E_L(r_p,E,L)
    v_p = sqrt(v2_p)
    vr_p = vr_from_E_L(r_p,E,L)
    vt_p = vt_from_E_L(r_p,E,L)

    v2_m = vSq_from_E_L(r_m,E,L)
    v_m = sqrt(v2_m)
    vr_m = vr_from_E_L(r_m,E,L)
    vt_m = vt_from_E_L(r_m,E,L)

    dvPar_rrp, _ = localVelChangeIso(r_p,vr_p,vt_p,m_field)
    dvPar_rrm, _ = localVelChangeIso(r_m,vr_m,vt_m,m_field)


    ratio_r = (L/v_p*dvPar_rrp - L/v_m*dvPar_rrm)/(2.0*eps*_b)

    return ratio_r
end

function LOverv_r(r::Float64,E::Float64,L::Float64,m_field::Float64)

    vr = vr_from_E_L(r,E,L)
    vt = vt_from_E_L(r,E,L)
    v2 = vSq_from_E_L(r,E,L)
    v = sqrt(v2)
    vr_v = vr/v
    vt_v = vt/v

    dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt, dvPar_r, dvParSq_r, dvPerSq_r = localVelChangeIsoGrad(r,vr,vt,m_field)
    dvPar_r = dvPar_r - dpsiEffdr(r,L)/vr*dvPar_vr - L/r^2*dvPar_vt


    ratio_r = L*dpsidr(r)/v^3*dvPar + L/v*dvPar_r

    return ratio_r
end

function LOverv_r_comp(r::Float64,E::Float64,L::Float64,m_field::Float64,eps::Float64=10^(-5))

    ratio_r = LOverv_r(r,E,L,m_field)
    ratio_r_num = LOverv_r_num(r,E,L,m_field,eps)

    println("ratio_r")
    println("Analytical = ",ratio_r)
    println("Numerical  = ",ratio_r_num)

end






function localOrbitalChangeEL(r::Float64, E::Float64, L::Float64, Lz::Float64,
                        m_field::Float64,
                        nbK::Int64=nbK_default, m_test::Float64=m_field)

    vr = vr_from_E_L(r,E,L)
    vt = vt_from_E_L(r,E,L)
    v = sqrt(vr^2 + vt^2)
    cosI = Lz/L

    # set theta=pi/4 so that sin^2(theta)=1/2, i.e. the mean value after angle-average
    coeffs = localOrbitalChange(r,0.25*pi,vr,vt,cosI,m_field,0.0,nbK,m_test)

    return coeffs
end

function dELLz_Iso_Num(r::Float64, E::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, m_test::Float64=m_field, eps::Float64=10^(-5))

    E_p = E + eps*abs(_E0)
    E_m = E - eps*abs(_E0)
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    r_p = r + eps*_b
    r_m = r - eps*_b

    # r-grad

    dE_rp, dL_rp, dLz_rp, dEE_rp, dLL_rp, dLzLz_rp, dEL_rp, dELz_rp, dLLz_rp = localOrbitalChangeEL(r_p,E,L,Lz,m_field,nbK,m_test)
    dE_rm, dL_rm, dLz_rm, dEE_rm, dLL_rm, dLzLz_rm, dEL_rm, dELz_rm, dLLz_rm = localOrbitalChangeEL(r_m,E,L,Lz,m_field,nbK,m_test)

    dE_r = (dE_rp - dE_rm)/(2.0*eps*_b)
    dL_r = (dL_rp - dL_rm)/(2.0*eps*_b)
    dLz_r = (dLz_rp - dLz_rm)/(2.0*eps*_b)
    dEE_r = (dEE_rp - dEE_rm)/(2.0*eps*_b)
    dLL_r = (dLL_rp - dLL_rm)/(2.0*eps*_b)
    dLzLz_r = (dLzLz_rp - dLzLz_rm)/(2.0*eps*_b)
    dEL_r = (dEL_rp - dEL_rm)/(2.0*eps*_b)
    dELz_r = (dELz_rp - dELz_rm)/(2.0*eps*_b)
    dLLz_r = (dLLz_rp - dLLz_rm)/(2.0*eps*_b)

    # E-grad

    dE_Ep, dL_Ep, dLz_Ep, dEE_Ep, dLL_Ep, dLzLz_Ep, dEL_Ep, dELz_Ep, dLLz_Ep = localOrbitalChangeEL(r,E_p,L,Lz,m_field,nbK,m_test)
    dE_Em, dL_Em, dLz_Em, dEE_Em, dLL_Em, dLzLz_Em, dEL_Em, dELz_Em, dLLz_Em = localOrbitalChangeEL(r,E_m,L,Lz,m_field,nbK,m_test)

    dE_E = (dE_Ep - dE_Em)/(2.0*eps*abs(_E0))
    dL_E = (dL_Ep - dL_Em)/(2.0*eps*abs(_E0))
    dLz_E = (dLz_Ep - dLz_Em)/(2.0*eps*abs(_E0))
    dEE_E = (dEE_Ep - dEE_Em)/(2.0*eps*abs(_E0))
    dLL_E = (dLL_Ep - dLL_Em)/(2.0*eps*abs(_E0))
    dLzLz_E = (dLzLz_Ep - dLzLz_Em)/(2.0*eps*abs(_E0))
    dEL_E = (dEL_Ep - dEL_Em)/(2.0*eps*abs(_E0))
    dELz_E = (dELz_Ep - dELz_Em)/(2.0*eps*abs(_E0))
    dLLz_E = (dLLz_Ep - dLLz_Em)/(2.0*eps*abs(_E0))

    # L-grad

    dE_Lp, dL_Lp, dLz_Lp, dEE_Lp, dLL_Lp, dLzLz_Lp, dEL_Lp, dELz_Lp, dLLz_Lp = localOrbitalChangeEL(r,E,L_p,Lz,m_field,nbK,m_test)
    dE_Lm, dL_Lm, dLz_Lm, dEE_Lm, dLL_Lm, dLzLz_Lm, dEL_Lm, dELz_Lm, dLLz_Lm = localOrbitalChangeEL(r,E,L_m,Lz,m_field,nbK,m_test)

    dE_L = (dE_Lp - dE_Lm)/(2.0*eps*abs(_L0))
    dL_L = (dL_Lp - dL_Lm)/(2.0*eps*abs(_L0))
    dLz_L = (dLz_Lp - dLz_Lm)/(2.0*eps*abs(_L0))
    dEE_L = (dEE_Lp - dEE_Lm)/(2.0*eps*abs(_L0))
    dLL_L = (dLL_Lp - dLL_Lm)/(2.0*eps*abs(_L0))
    dLzLz_L = (dLzLz_Lp - dLzLz_Lm)/(2.0*eps*abs(_L0))
    dEL_L = (dEL_Lp - dEL_Lm)/(2.0*eps*abs(_L0))
    dELz_L = (dELz_Lp - dELz_Lm)/(2.0*eps*abs(_L0))
    dLLz_L = (dLLz_Lp - dLLz_Lm)/(2.0*eps*abs(_L0))

    # Lz-grad

    _, _, dLz_Lzp, _, _, dLzLz_Lzp,  _, dELz_Lzp, dLLz_Lzp = localOrbitalChangeEL(r,E,L,Lz_p,m_field,nbK,m_test)
    _, _, dLz_Lzm, _, _, dLzLz_Lzm,  _, dELz_Lzm, dLLz_Lzm = localOrbitalChangeEL(r,E,L,Lz_m,m_field,nbK,m_test)

    dLz_Lz = (dLz_Lzp - dLz_Lzm)/(2.0*eps*abs(_L0))
    dLzLz_Lz = (dLzLz_Lzp - dLzLz_Lzm)/(2.0*eps*abs(_L0))
    dELz_Lz = (dELz_Lzp - dELz_Lzm)/(2.0*eps*abs(_L0))
    dLLz_Lz = (dLLz_Lzp - dLLz_Lzm)/(2.0*eps*abs(_L0))

    return dE_E, dE_L, dL_E, dL_L, dLz_E, dLz_L, dLz_Lz, dEE_E, dEE_L, dLL_E, dLL_L, dLzLz_E, dLzLz_L, dLzLz_Lz, dEL_E, dEL_L, dELz_E, dELz_L, dELz_Lz, dLLz_E, dLLz_L, dLLz_Lz, dE_r, dL_r, dLz_r, dEE_r, dLL_r, dLzLz_r, dEL_r, dELz_r, dLLz_r
end



function dELLz_Iso_Compare(r::Float64, E::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, m_test::Float64=m_field, eps::Float64=10^(-5))

    num = dELLz_Iso_Num(r,E,L,Lz,m_field,nbK,m_test,eps)
    ana = dELLz_Iso(r,E,L,Lz,m_field,nbK,m_test)

    println("dE_E")
    println("Analytical = ",ana[1])
    println("Numerical  = ",num[1])
    println("-------------------")
    println("dE_L")
    println("Analytical = ",ana[2])
    println("Numerical  = ",num[2])
    println("-------------------")
    println("dL_E")
    println("Analytical = ",ana[3])
    println("Numerical  = ",num[3])
    println("-------------------")
    println("dL_L")
    println("Analytical = ",ana[4])
    println("Numerical  = ",num[4])
    println("-------------------")
    println("dLz_E")
    println("Analytical = ",ana[5])
    println("Numerical  = ",num[5])
    println("-------------------")
    println("dLz_L")
    println("Analytical = ",ana[6])
    println("Numerical  = ",num[6])
    println("-------------------")
    println("dLz_Lz")
    println("Analytical = ",ana[7])
    println("Numerical  = ",num[7])
    println("-------------------")
    println("dEE_E")
    println("Analytical = ",ana[8])
    println("Numerical  = ",num[8])
    println("-------------------")
    println("dEE_L")
    println("Analytical = ",ana[9])
    println("Numerical  = ",num[9])
    println("-------------------")
    println("dLL_E")
    println("Analytical = ",ana[10])
    println("Numerical  = ",num[10])
    println("-------------------")
    println("dLL_L")
    println("Analytical = ",ana[11])
    println("Numerical  = ",num[11])
    println("-------------------")
    println("dLzLz_E")
    println("Analytical = ",ana[12])
    println("Numerical  = ",num[12])
    println("-------------------")
    println("dLzLz_L")
    println("Analytical = ",ana[13])
    println("Numerical  = ",num[13])
    println("-------------------")
    println("dLzLz_Lz")
    println("Analytical = ",ana[14])
    println("Numerical  = ",num[14])
    println("-------------------")
    println("dEL_E")
    println("Analytical = ",ana[15])
    println("Numerical  = ",num[15])
    println("-------------------")
    println("dEL_L")
    println("Analytical = ",ana[16])
    println("Numerical  = ",num[16])
    println("-------------------")
    println("dELz_E")
    println("Analytical = ",ana[17])
    println("Numerical  = ",num[17])
    println("-------------------")
    println("dELz_L")
    println("Analytical = ",ana[18])
    println("Numerical  = ",num[18])
    println("-------------------")
    println("dELz_Lz")
    println("Analytical = ",ana[19])
    println("Numerical  = ",num[19])
    println("-------------------")
    println("dLLz_E")
    println("Analytical = ",ana[20])
    println("Numerical  = ",num[20])
    println("-------------------")
    println("dLLz_L")
    println("Analytical = ",ana[21])
    println("Numerical  = ",num[21])
    println("-------------------")
    println("dLLz_Lz")
    println("Analytical = ",ana[22])
    println("Numerical  = ",num[22])
    println("-------------------")
    println("dE_r")
    println("Analytical = ",ana[23])
    println("Numerical  = ",num[23])
    println("-------------------")
    println("dL_r")
    println("Analytical = ",ana[24])
    println("Numerical  = ",num[24])
    println("-------------------")
    println("dLz_r")
    println("Analytical = ",ana[25])
    println("Numerical  = ",num[25])
    println("-------------------")
    println("dEE_r")
    println("Analytical = ",ana[26])
    println("Numerical  = ",num[26])
    println("-------------------")
    println("dLL_r")
    println("Analytical = ",ana[27])
    println("Numerical  = ",num[27])
    println("-------------------")
    println("dLzLz_r")
    println("Analytical = ",ana[28])
    println("Numerical  = ",num[28])
    println("-------------------")
    println("dEL_r")
    println("Analytical = ",ana[29])
    println("Numerical  = ",num[29])
    println("-------------------")
    println("dELz_r")
    println("Analytical = ",ana[30])
    println("Numerical  = ",num[30])
    println("-------------------")
    println("dLLz_r")
    println("Analytical = ",ana[31])
    println("Numerical  = ",num[31])

end

#############################################################
# Check orbit-average diffusion gradients Energy
#############################################################

# orbitAverageEnergyCoeffsNoRot(E,L,Lz/L,m_field,nbK,m_test)

function DELLz_Iso_OrbitAvg_Num(E::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbAvr::Int64=nbAvr_default, nbK::Int64=nbK_default, m_test::Float64=m_field, eps::Float64=10^(-5))

    E_p = E + eps*abs(_E0)
    E_m = E - eps*abs(_E0)
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0

    # E-grad

    dE_Ep, dL_Ep, dLz_Ep, _ = orbitAverageEnergyCoeffsNoRot(E_p,L,Lz/L,m_field,nbAvr,nbK,m_test)
    dE_Em, dL_Em, dLz_Em, _ = orbitAverageEnergyCoeffsNoRot(E_m,L,Lz/L,m_field,nbAvr,nbK,m_test)

    dE_E = (dE_Ep - dE_Em)/(2.0*eps*abs(_E0))
    dL_E = (dL_Ep - dL_Em)/(2.0*eps*abs(_E0))
    dLz_E = (dLz_Ep - dLz_Em)/(2.0*eps*abs(_E0))


    # L-grad

    dE_Lp, dL_Lp, dLz_Lp, _ = orbitAverageEnergyCoeffsNoRot(E,L_p,Lz/L_p,m_field,nbAvr,nbK,m_test)
    dE_Lm, dL_Lm, dLz_Lm, _ = orbitAverageEnergyCoeffsNoRot(E,L_m,Lz/L_m,m_field,nbAvr,nbK,m_test)

    dE_L = (dE_Lp - dE_Lm)/(2.0*eps*abs(_L0))
    dL_L = (dL_Lp - dL_Lm)/(2.0*eps*abs(_L0))
    dLz_L = (dLz_Lp - dLz_Lm)/(2.0*eps*abs(_L0))

    # Lz-grad

    _, _, dLz_Lzp, _ = orbitAverageEnergyCoeffsNoRot(E,L,Lz_p/L,m_field,nbAvr,nbK,m_test)
    _, _, dLz_Lzm, _ = orbitAverageEnergyCoeffsNoRot(E,L,Lz_m/L,m_field,nbAvr,nbK,m_test)

    dLz_Lz = (dLz_Lzp - dLz_Lzm)/(2.0*eps*abs(_L0))


    return dE_E, dE_L, dL_E, dL_L, dLz_E, dLz_L, dLz_Lz
end

function DELLz_Iso_OrbitAvg_Compare(E::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field,
            eps::Float64=10^(-5))

    num = DELLz_Iso_OrbitAvg_Num(E,L,Lz,m_field,nbK,nbAvr,m_test,eps)
    ana = DELLz_Iso_OrbitAvg(E,L,Lz,m_field,nbK,nbAvr,m_test)

    println("dE_E")
    println("Analytical = ",ana[1])
    println("Numerical  = ",num[1])
    println("-------------------")
    println("dE_L")
    println("Analytical = ",ana[2])
    println("Numerical  = ",num[2])
    println("-------------------")
    println("dL_E")
    println("Analytical = ",ana[3])
    println("Numerical  = ",num[3])
    println("-------------------")
    println("dL_L")
    println("Analytical = ",ana[4])
    println("Numerical  = ",num[4])
    println("-------------------")
    println("dLz_E")
    println("Analytical = ",ana[5])
    println("Numerical  = ",num[5])
    println("-------------------")
    println("dLz_L")
    println("Analytical = ",ana[6])
    println("Numerical  = ",num[6])
    println("-------------------")
    println("dLz_Lz")
    println("Analytical = ",ana[7])
    println("Numerical  = ",num[7])
    println("-------------------")

end

#############################################################
# Check orbit-average diffusion gradients actions
#############################################################

function DJrLLz_Iso_OrbitAvg_Num(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64=nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

    Jr_p = Jr + eps*_L0
    Jr_m = Jr - eps*_L0
    L_p = L + eps*_L0
    L_m = L - eps*_L0
    Lz_p = Lz + eps*_L0
    Lz_m = Lz - eps*_L0


    # Jr-grad
    _, _, _, dJrJr_Jrp, _ = orbitAverageActionCoeffsNoRot(Jr_p,L,Lz/L,m_field,nbAvr,nbK,nbu,m_test)
    _, _, _, dJrJr_Jrm, _ = orbitAverageActionCoeffsNoRot(Jr_m,L,Lz/L,m_field,nbAvr,nbK,nbu,m_test)

    dJrJr_Jr = (dJrJr_Jrp - dJrJr_Jrm)/(2.0*eps*_L0)


    # L-grad
    _, _, _, dJrJr_Lp, _ = orbitAverageActionCoeffsNoRot(Jr,L_p,Lz/L_p,m_field,nbAvr,nbK,nbu,m_test)
    _, _, _, dJrJr_Lm, _ = orbitAverageActionCoeffsNoRot(Jr,L_m,Lz/L_m,m_field,nbAvr,nbK,nbu,m_test)

    dJrJr_L = (dJrJr_Lp - dJrJr_Lm)/(2.0*eps*_L0)


    # Lz-grad
    _, _, _, dJrJr_Lzp, _ = orbitAverageActionCoeffsNoRot(Jr,L,Lz_p/L,m_field,nbAvr,nbK,nbu,m_test)
    _, _, _, dJrJr_Lzm, _ = orbitAverageActionCoeffsNoRot(Jr,L,Lz_m/L,m_field,nbAvr,nbK,nbu,m_test)




    return dJrJr_Jr, dJrJr_L
end



function DJrLLz_Iso_OrbitAvg_Compare(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64=nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))



    num = DJrLLz_Iso_OrbitAvg_Num(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test,eps)
    ana = DJrLLz_Iso_OrbitAvg(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test)


    println("dJrJr_Jr")
    println("Analytical = ",ana[1])
    println("Numerical  = ",num[1])
    println("-------------------")
    println("dJrJr_L")
    println("Analytical = ",ana[2])
    println("Numerical  = ",num[2])
end

#############################################################
# Check flux
#############################################################
function flux_Iso_numerical(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

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

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsNoRot(Jr,L,Lz/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot= _F(E,L)




    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsNoRot(Jr_p,L,Lz/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Jrp = _F(E_Jrp,L)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsNoRot(Jr_m,L,Lz/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Jrm = _F(E_Jrm,L)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsNoRot(Jr,L_p,Lz/L_p,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Lp = _F(E_Lp,L_p)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsNoRot(Jr,L_m,Lz/L_m,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Lm = _F(E_Lm,L_m)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsNoRot(Jr,L,Lz_p/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Lzp = _F(E,L)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsNoRot(Jr,L,Lz_m/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Lzm = _F(E,L)

    # Jr-component

    DJrF = dJr*Ftot

    dJrJrF_Jr = (dJrJr_Jrp*Ftot_Jrp-dJrJr_Jrm*Ftot_Jrm)/(2.0*eps*_L0)
    dJrLF_L = (dJrL_Lp*Ftot_Lp-dJrL_Lm*Ftot_Lm)/(2.0*eps*_L0)
    dJrLzF_Lz = (dJrLz_Lzp*Ftot_Lzp-dJrLz_Lzm*Ftot_Lzm)/(2.0*eps*_L0)

    fluxJr = DJrF - 0.5*(dJrJrF_Jr+dJrLF_L+dJrLzF_Lz)

    # L-component

    DLF = dL*Ftot

    dJrLF_Jr = (dJrL_Jrp*Ftot_Jrp-dJrL_Jrm*Ftot_Jrm)/(2.0*eps*_L0)
    dLLF_L = (dLL_Lp*Ftot_Lp-dLL_Lm*Ftot_Lm)/(2.0*eps*_L0)
    dLLzF_Lz = (dLLz_Lzp*Ftot_Lzp-dLLz_Lzm*Ftot_Lzm)/(2.0*eps*_L0)

    fluxL = DLF - 0.5*(dJrLF_Jr+dLLF_L+dLLzF_Lz)

    # L-component

    DLzF = dLz*Ftot

    dJrLzF_Jr = (dJrLz_Jrp*Ftot_Jrp-dJrLz_Jrm*Ftot_Jrm)/(2.0*eps*_L0)
    dLLzF_L = (dLLz_Lp*Ftot_Lp-dLLz_Lm*Ftot_Lm)/(2.0*eps*_L0)
    dLzLzF_Lz = (dLzLz_Lzp*Ftot_Lzp-dLzLz_Lzm*Ftot_Lzm)/(2.0*eps*_L0)

    fluxLz = DLzF - 0.5*(dJrLzF_Jr+dLLzF_L+dLzLzF_Lz)

    return fluxJr, fluxL, fluxLz
end



function dFdt_Iso_numerical(Jr::Float64, L::Float64, Lz::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-5))

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

    dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsNoRot(Jr,L,Lz/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot= _F(E,L)

    # Partial derivatives

    dJr_Jrp, dL_Jrp, dLz_Jrp, dJrJr_Jrp, dLL_Jrp, dLzLz_Jrp, dJrL_Jrp, dJrLz_Jrp, dLLz_Jrp = orbitAverageActionCoeffsNoRot(Jr_p,L,Lz/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Jrp = _F(E_Jrp,L)
    dJr_Jrm, dL_Jrm, dLz_Jrm, dJrJr_Jrm, dLL_Jrm, dLzLz_Jrm, dJrL_Jrm, dJrLz_Jrm, dLLz_Jrm = orbitAverageActionCoeffsNoRot(Jr_m,L,Lz/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Jrm = _F(E_Jrm,L)

    dJr_Lp, dL_Lp, dLz_Lp, dJrJr_Lp, dLL_Lp, dLzLz_Lp, dJrL_Lp, dJrLz_Lp, dLLz_Lp = orbitAverageActionCoeffsNoRot(Jr,L_p,Lz/L_p,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Lp = _F(E_Lp,L_p)
    dJr_Lm, dL_Lm, dLz_Lm, dJrJr_Lm, dLL_Lm, dLzLz_Lm, dJrL_Lm, dJrLz_Lm, dLLz_Lm = orbitAverageActionCoeffsNoRot(Jr,L_m,Lz/L_m,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Lm = _F(E_Lm,L_m)

    dJr_Lzp, dL_Lzp, dLz_Lzp, dJrJr_Lzp, dLL_Lzp, dLzLz_Lzp, dJrL_Lzp, dJrLz_Lzp, dLLz_Lzp = orbitAverageActionCoeffsNoRot(Jr,L,Lz_p/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Lzp = _F(E,L)
    dJr_Lzm, dL_Lzm, dLz_Lzm, dJrJr_Lzm, dLL_Lzm, dLzLz_Lzm, dJrL_Lzm, dJrLz_Lzm, dLLz_Lzm = orbitAverageActionCoeffsNoRot(Jr,L,Lz_m/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_Lzm = _F(E,L)


    DJrF_Jr = (dJr_Jrp*Ftot_Jrp-dJr_Jrm*Ftot_Jrm)/(2.0*eps*_L0)
    DLF_L = (dL_Lp*Ftot_Lp-dL_Lm*Ftot_Lm)/(2.0*eps*_L0)
    DLzF_Lz = (dLz_Lzp*Ftot_Lzp-dLz_Lzm*Ftot_Lzm)/(2.0*eps*_L0)

    DJrJr_F_JrJr = (dJrJr_Jrp*Ftot_Jrp + dJrJr_Jrm*Ftot_Jrm - 2.0*dJrJr*Ftot)/(eps*_L0)^2
    DLL_F_LL = (dLL_Lp*Ftot_Lp + dLL_Lm*Ftot_Lm - 2.0*dLL*Ftot)/(eps*_L0)^2
    DLzLz_F_LzLz = (dLzLz_Lzp*Ftot_Lzp + dLzLz_Lzm*Ftot_Lzm - 2.0*dLzLz*Ftot)/(eps*_L0)^2





    # Mixed derivatives

    dJr_JrpLp, dL_JrpLp, dLz_JrpLp, dJrJr_JrpLp, dLL_JrpLp, dLzLz_JrpLp, dJrL_JrpLp, dJrLz_JrpLp, dLLz_JrpLp = orbitAverageActionCoeffsNoRot(Jr_p,L_p,Lz/L_p,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLp = _F(E_Jrp_Lp,L_p)
    dJr_JrmLm, dL_JrmLm, dLz_JrmLm, dJrJr_JrmLm, dLL_JrmLm, dLzLz_JrmLm, dJrL_JrmLm, dJrLz_JrmLm, dLLz_JrmLm = orbitAverageActionCoeffsNoRot(Jr_m,L_m,Lz/L_m,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLm = _F(E_Jrm_Lm,L_m)

    dJr_JrpLm, dL_JrpLm, dLz_JrpLm, dJrJr_JrpLm, dLL_JrpLm, dLzLz_JrpLm, dJrL_JrpLm, dJrLz_JrpLm, dLLz_JrpLm = orbitAverageActionCoeffsNoRot(Jr_p,L_m,Lz/L_m,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLm = _F(E_Jrp_Lm,L_m)
    dJr_JrmLp, dL_JrmLp, dLz_JrmLp, dJrJr_JrmLp, dLL_JrmLp, dLzLz_JrmLp, dJrL_JrmLp, dJrLz_JrmLp, dLLz_JrmLp = orbitAverageActionCoeffsNoRot(Jr_m,L_p,Lz/L_p,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLp = _F(E_Jrm_Lp,L_p)

    dJr_JrpLzp, dL_JrpLzp, dLz_JrpLzp, dJrJr_JrpLzp, dLL_JrpLzp, dLzLz_JrpLzp, dJrL_JrpLzp, dJrLz_JrpLzp, dLLz_JrpLzp = orbitAverageActionCoeffsNoRot(Jr_p,L,Lz_p/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLzp = _F(E_Jrp,L)
    dJr_JrmLzm, dL_JrmLzm, dLz_JrmLzm, dJrJr_JrmLzm, dLL_JrmLzm, dLzLz_JrmLzm, dJrL_JrmLzm, dJrLz_JrmLzm, dLLz_JrmLzm = orbitAverageActionCoeffsNoRot(Jr_m,L,Lz_m/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLzm = _F(E_Jrm,L)

    dJr_JrpLzm, dL_JrpLzm, dLz_JrpLzm, dJrJr_JrpLzm, dLL_JrpLzm, dLzLz_JrpLzm, dJrL_JrpLzm, dJrLz_JrpLzm, dLLz_JrpLzm = orbitAverageActionCoeffsNoRot(Jr_p,L,Lz_m/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_JrpLzm = _F(E_Jrp,L)
    dJr_JrmLzp, dL_JrmLzp, dLz_JrmLzp, dJrJr_JrmLzp, dLL_JrmLzp, dLzLz_JrmLzp, dJrL_JrmLzp, dJrLz_JrmLzp, dLLz_JrmLzp = orbitAverageActionCoeffsNoRot(Jr_m,L,Lz_p/L,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_JrmLzp = _F(E_Jrm,L)

    dJr_LpLzp, dL_LpLzp, dLz_LpLzp, dJrJr_LpLzp, dLL_LpLzp, dLzLz_LpLzp, dJrL_LpLzp, dJrLz_LpLzp, dLLz_LpLzp = orbitAverageActionCoeffsNoRot(Jr,L_p,Lz_p/L_p,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_LpLzp = _F(E_Lp,L_p)
    dJr_LmLzm, dL_LmLzm, dLz_LmLzm, dJrJr_LmLzm, dLL_LmLzm, dLzLz_LmLzm, dJrL_LmLzm, dJrLz_LmLzm, dLLz_LmLzm = orbitAverageActionCoeffsNoRot(Jr,L_m,Lz_m/L_m,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_LmLzm = _F(E_Lm,L_m)

    dJr_LpLzm, dL_LpLzm, dLz_LpLzm, dJrJr_LpLzm, dLL_LpLzm, dLzLz_LpLzm, dJrL_LpLzm, dJrLz_LpLzm, dLLz_LpLzm = orbitAverageActionCoeffsNoRot(Jr,L_p,Lz_m/L_p,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_LpLzm = _F(E_Lp,L_p)
    dJr_LmLzp, dL_LmLzp, dLz_LmLzp, dJrJr_LmLzp, dLL_LmLzp, dLzLz_LmLzp, dJrL_LmLzp, dJrLz_LmLzp, dLLz_LmLzp = orbitAverageActionCoeffsNoRot(Jr,L_m,Lz_p/L_m,m_field,nbAvr,nbK,nbu,m_test)
    Ftot_LmLzp = _F(E_Lm,L_m)




    DJrL_F_JrL = (dJrL_JrpLp*Ftot_JrpLp+dJrL_JrmLm*Ftot_JrmLm-dJrL_JrpLm*Ftot_JrpLm-dJrL_JrmLp*Ftot_JrmLp)/(2.0*eps*_L0)^2

    DJrLz_F_JrLz = (dJrLz_JrpLzp*Ftot_JrpLzp+dJrLz_JrmLzm*Ftot_JrmLzm-dJrLz_JrpLzm*Ftot_JrpLzm-dJrLz_JrmLzp*Ftot_JrmLzp)/(2.0*eps*_L0)^2

    DLLz_F_LLz = (dLLz_LpLzp*Ftot_LpLzp+dLLz_LmLzm*Ftot_LmLzm-dLLz_LpLzm*Ftot_LpLzm-dLLz_LmLzp*Ftot_LmLzp)/(2.0*eps*_L0)^2


    fluxJr_Jr = DJrF_Jr - 0.5*(DJrJr_F_JrJr + DJrL_F_JrL + DJrLz_F_JrLz)
    fluxL_L = DLF_L - 0.5*(DJrL_F_JrL + DLL_F_LL + DLLz_F_LLz)
    fluxLz_Lz = DLzF_Lz - 0.5*(DJrLz_F_JrLz + DLLz_F_LLz + DLzLz_F_LzLz)

    # println("Jr = ",fluxJr_Jr)
    # println("L  = ",fluxL_L)
    # println("Lz = ",fluxLz_Lz)

    return -(fluxJr_Jr+fluxL_L+fluxLz_Lz)
end

#############################################################
# Check orbit-average relevant quantities: T/2, Theta(u) and their gradients
# And numerator in 2/T * int{Theta(u) dE}
# For the latter one, just check the derivative at some u-point ?
# This derivative is wrong for the moment
#############################################################

# no correct. d/dr gradients in jacobians must be included
function halfperiod_and_grads(E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field)

   halfT = 0.0
   halfT_E = 0.0
   halfT_L = 0.0

   # int Theta(u)*dE
   numerator = 0.0
   numerator_E = 0.0
   numerator_L = 0.0

   sp, sa = radius_s_bounds(E,L)
   sma, ecc = sma_ecc_from_sp_sa(sp,sa)

   for iu=1:nbAvr
       uloc = -1.0+2.0/nbAvr*(iu-0.5)
       sloc = s_from_u_sma_ecc(uloc,sma,ecc)
       rloc = r_from_s(sloc)
       jac_loc = Theta(uloc,sp,sa)

       djacdE, djacdL, _ = djac_and_ds(uloc,sp,sa)

       vr = sqrt(2.0*(E - psiEff(rloc,L)))
       vt = L/rloc

       halfT += jac_loc
       halfT_E += djacdE
       halfT_L += djacdL

       dE, _ = localOrbitChangeNoRot(rloc,vr,vt,1.0,m_field,nbK,m_test)
       dE_E,dE_L, _ = dELLz_Iso(rloc,E,L,L,m_field,nbK,m_test)

       # println("dE = ",dE)
       # println("dE_E, dE_L = ",(dE_E,dE_L))



       numerator += jac_loc*dE
       numerator_E += jac_loc*dE_E + djacdE*dE
       numerator_L += jac_loc*dE_L + djacdL*dE
   end

   halfT *= 2.0/nbAvr
   halfT_E *= 2.0/nbAvr
   halfT_L *= 2.0/nbAvr

   numerator *= 2.0/nbAvr
   numerator_E *= 2.0/nbAvr
   numerator_L *= 2.0/nbAvr

   return halfT, halfT_E , halfT_L, numerator, numerator_E, numerator_L
end

function halfperiod_grads_num(E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field,
            eps::Float64=10^(-5))

    E_p = E + eps*abs(_E0)
    E_m = E - eps*abs(_E0)

    L_p = L + eps*_L0
    L_m = L - eps*_L0

    # E-grad

    halfT_E_p, _, _, numerator_E_p, _ = halfperiod_and_grads(E_p,L,m_field,nbK,nbAvr,m_test)
    halfT_E_m, _, _, numerator_E_m, _ = halfperiod_and_grads(E_m,L,m_field,nbK,nbAvr,m_test)

    halfT_E = (halfT_E_p-halfT_E_m)/(2.0*eps*abs(_E0))
    numerator_E = (numerator_E_p-numerator_E_m)/(2.0*eps*abs(_E0))

    # L-grad

    halfT_L_p, _, _, numerator_L_p, _ = halfperiod_and_grads(E,L_p,m_field,nbK,nbAvr,m_test)
    halfT_L_m, _, _, numerator_L_m, _ = halfperiod_and_grads(E,L_m,m_field,nbK,nbAvr,m_test)

    halfT_L = (halfT_L_p-halfT_L_m)/(2.0*eps*_L0)
    numerator_L = (numerator_L_p - numerator_L_m)/(2.0*eps*_L0)

    return halfT_E, halfT_L, numerator_E, numerator_L
end

function halfperiod_grads_compare(E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field,
            eps::Float64=10^(-5))

    _, halfT_E, halfT_L, _, numerator_E, numerator_L = halfperiod_and_grads(E,L,m_field,nbK,nbAvr,m_test)
    halfT_E_num, halfT_L_num, numerator_E_num, numerator_L_num = halfperiod_grads_num(E,L,m_field,nbK,nbAvr,m_test,eps)

    println("halfT_E")
    println("Analytical = ",halfT_E)
    println("Numerical  = ",halfT_E_num)
    println("-------------------")
    println("halfT_L")
    println("Analytical = ",halfT_L)
    println("Numerical  = ",halfT_L_num)
    println("-------------------")
    println("numerator_E")
    println("Analytical = ",numerator_E)
    println("Numerical  = ",numerator_E_num)
    println("-------------------")
    println("numerator_L")
    println("Analytical = ",numerator_L)
    println("Numerical  = ",numerator_L_num)
    println("-------------------")
end

function integrand_numerator(u::Float64, E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field)

    sp, sa = radius_s_bounds(E,L)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)

    sloc = s_from_u_sma_ecc(u,sma,ecc)
    rloc = r_from_s(sloc)
    jac_loc = Theta(u,sp,sa)

    djacdE, djacdL, _ = djac_and_ds(u,sp,sa)

    vr = sqrt(2.0*(E - psiEff(rloc,L)))
    vt = L/rloc


    dE, _ = localOrbitChangeNoRot(rloc,vr,vt,1.0,m_field,nbK,m_test)
    dE_E,dE_L, _ = dELLz_Iso(rloc,E,L,L,m_field,nbK,m_test)

    numerator = jac_loc*dE


    numerator_E = jac_loc*dE_E + djacdE*dE
    numerator_L = jac_loc*dE_L + djacdL*dE

    return numerator, numerator_E, numerator_L, jac_loc, djacdE, djacdL
end

function integrand_numerator_num(u::Float64, E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field,
            eps::Float64=10^(-5))

    E_p = E + eps*abs(_E0)
    E_m = E - eps*abs(_E0)

    L_p = L + eps*_L0
    L_m = L - eps*_L0

    # E-grad

    num_E_p, _, _, jac_E_p, _ = integrand_numerator(u,E_p,L,m_field,nbK,nbAvr,m_test)
    num_E_m, _, _, jac_E_m, _ = integrand_numerator(u,E_m,L,m_field,nbK,nbAvr,m_test)

    num_E = (num_E_p-num_E_m)/(2.0*eps*abs(_E0))
    jac_E = (jac_E_p-jac_E_m)/(2.0*eps*abs(_E0))

    # L-grad

    num_L_p, _, _, jac_L_p, _ = integrand_numerator(u,E,L_p,m_field,nbK,nbAvr,m_test)
    num_L_m, _, _, jac_L_m, _ = integrand_numerator(u,E,L_m,m_field,nbK,nbAvr,m_test)

    num_L = (num_L_p-num_L_m)/(2.0*eps*_L0)
    jac_L = (jac_L_p-jac_L_m)/(2.0*eps*_L0)

    return num_E, num_L, jac_E, jac_L
end

function integrand_numerator_comp(u::Float64, E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, m_test::Float64=m_field,
            eps::Float64=10^(-5))

    _, num_E, num_L, _, jac_E, jac_L = integrand_numerator(u,E,L,m_field,nbK,nbAvr,m_test)
    num_E_num, num_L_num, jac_E_num, jac_L_num = integrand_numerator_num(u,E,L,m_field,nbK,nbAvr,m_test,eps)

    println("numerator_E")
    println("Analytical = ",num_E)
    println("Numerical  = ",num_E_num)
    println("-------------------")
    println("numerator_L")
    println("Analytical = ",num_L)
    println("Numerical  = ",num_L_num)
    println("-------------------")
    println("jac_E")
    println("Analytical = ",jac_E)
    println("Numerical  = ",jac_E_num)
    println("-------------------")
    println("jac_L")
    println("Analytical = ",jac_L)
    println("Numerical  = ",jac_L_num)
    println("-------------------")

end

#############################################################
# Check dvParSq
#############################################################

function dvParSqChainRule(r::Float64, E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, m_test::Float64=m_field, eps::Float64=10^(-5))

    vr = vr_from_E_L(r,E,L)
    vt = vt_from_E_L(r,E,L)
    v2 = vSq_from_E_L(r,E,L)
    v = sqrt(v2)

    println("(vr,vt) = ",(vr,vt))
    println("v^2_from_EL =",v2)
    println("vr^2+vt^2   =",vr^2+vt^2)

    # dvPar

    dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt = localVelChangeIsoGrad(r,vr,vt,m_field,nbK,m_test)
    println("1/vr d(dvPar)/dvr     = ",1.0/vr*dvPar_vr)

    dvPar_vr_num, dvParSq_vr_num, dvPerSq_vr_num, dvPar_vt_num, dvParSq_vt_num, dvPerSq_vt_num = localVelChangeIsoNum(r,vr,vt,m_field,nbK,m_test,eps)
    println("1/vr d(dvPar)/dvr_num = ",1.0/vr*dvPar_vr_num)

    dvPar_E_num, dvPar_L_num, dvParSq_E_num, dvParSq_L_num, dvPerSq_E_num, dvPerSq_L_num = dvPar_EL_num(r,E,L,m_field,nbK,m_test,eps)
    println("d(dvPar)/dE_num       = ",dvPar_E_num)

    # dvParSq

    dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt = localVelChangeIsoGrad(r,vr,vt,m_field,nbK,m_test)
    println("1/vr d(dvParSq)/dvr     = ",1.0/vr*dvParSq_vr)

    dvPar_vr_num, dvParSq_vr_num, dvPerSq_vr_num, dvPar_vt_num, dvParSq_vt_num, dvPerSq_vt_num = localVelChangeIsoNum(r,vr,vt,m_field,nbK,m_test,eps)
    println("1/vr d(dvParSq)/dvr_num = ",1.0/vr*dvParSq_vr_num)

    dvPar_E_num, dvPar_L_num, dvParSq_E_num, dvParSq_L_num, dvPerSq_E_num, dvPerSq_L_num = dvPar_EL_num(r,E,L,m_field,nbK,m_test,eps)
    println("d(dvParSq)/dE_num       = ",dvParSq_E_num)

    # dvPerSq

    dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt = localVelChangeIsoGrad(r,vr,vt,m_field,nbK,m_test)
    println("1/vr d(dvPerSq)/dvr     = ",1.0/vr*dvPerSq_vr)

    dvPar_vr_num, dvParSq_vr_num, dvPerSq_vr_num, dvPar_vt_num, dvParSq_vt_num, dvPerSq_vt_num = localVelChangeIsoNum(r,vr,vt,m_field,nbK,m_test,eps)
    println("1/vr d(dvPerSq)/dvr_num = ",1.0/vr*dvPerSq_vr_num)

    dvPar_E_num, dvPar_L_num, dvParSq_E_num, dvParSq_L_num, dvPerSq_E_num, dvPerSq_L_num = dvPar_EL_num(r,E,L,m_field,nbK,m_test,eps)
    println("d(dvPerSq)/dE_num       = ",dvPerSq_E_num)


    # v*dvPar

    dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt = localVelChangeIsoGrad(r,vr,vt,m_field,nbK,m_test)
    println("1/vr d(v*dvPar)/dvr     = ",1.0/vr*dvPar_vr + dvPar/v)

    dvPar_vr_num, dvParSq_vr_num, dvPerSq_vr_num, dvPar_vt_num, dvParSq_vt_num, dvPerSq_vt_num = localVelChangeIsoNum(r,vr,vt,m_field,nbK,m_test,eps)
    println("1/vr d(dvPar)/dvr_num   = ",1.0/vr*dvPar_vr_num + dvPar/v)



    #num derivative of v*dvPar


    E_p = E + eps*abs(_E0)
    E_m = E - eps*abs(_E0)

    vr_Ep = vr_from_E_L(r,E_p,L)
    vr_Em = vr_from_E_L(r,E_m,L)

    vt_Ep = vt_from_E_L(r,E_p,L)
    vt_Em = vt_from_E_L(r,E_m,L)


    v_Ep = sqrt(vr_Ep^2 + vt_Ep^2)
    v_Em = sqrt(vr_Em^2 + vt_Em^2)


    dvPar_Ep, _ = localVelChangeIsoGrad(r,vr_Ep,vt_Ep,m_field,nbK,m_test)
    dvPar_Em, _ = localVelChangeIsoGrad(r,vr_Em,vt_Em,m_field,nbK,m_test)



    # println((dvPar_Lp, dvParSq_Lp, dvPerSq_Lp))
    # println((dvPar_Lm, dvParSq_Lm, dvPerSq_Lm))

    # E-grad

    vdvPar_E_num = (v_Ep*dvPar_Ep-v_Em*dvPar_Em)/(2.0*eps*abs(_E0))





    println("d(dvPar)/dE_num         = ",vdvPar_E_num)


    # v

    E_p = E + eps*abs(_E0)
    E_m = E - eps*abs(_E0)

    v_p = sqrt(2.0*(E_p-psi(r)))
    v_m = sqrt(2.0*(E_m-psi(r)))

    dvdE_num = (v_p-v_m)/(2.0*eps*abs(_E0))

    println("dv/dE     = ",1.0/v)
    println("dv/dE_num = ",dvdE_num)



end

#############################################################
# Check E,L,Lz gradients of dv locaux
#############################################################

function dvPar_EL(r::Float64, E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, m_test::Float64=m_field)

    vr = vr_from_E_L(r,E,L)
    vt = vt_from_E_L(r,E,L)
    v2 = vSq_from_E_L(r,E,L)
    v = sqrt(v2)
    vr_v = vr/v
    vt_v = vt/v

    dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt, dvPar_r, dvParSq_r, dvPerSq_r = localVelChangeIsoGrad(r,vr,vt,m_field,nbK,m_test)

    # println((dvPar, dvParSq, dvPerSq, dvPar_vr, dvParSq_vr, dvPerSq_vr, dvPar_vt, dvParSq_vt, dvPerSq_vt))

    dvPar_E = dvPar_vr/vr
    dvPar_L = -L/(r^2*vr)*dvPar_vr + 1.0/r*dvPar_vt

    # println((dvPar_E,dvPar_L))

    dvParSq_E = dvParSq_vr/vr
    dvParSq_L = -L/(r^2*vr)*dvParSq_vr + 1.0/r*dvParSq_vt



    # println((dvParSq_E,dvParSq_L))

    # println(-L/(r^2*vr)*dvParSq_vr)
    # println(1.0/r*dvParSq_vt)

    dvPerSq_E = dvPerSq_vr/vr
    dvPerSq_L = -L/(r^2*vr)*dvPerSq_vr + 1.0/r*dvPerSq_vt

    dvParSq_r = dvParSq_r - dpsiEffdr(r,L)/vr*dvParSq_vr - L/r^2*dvParSq_vt
    dvPerSq_r = dvPerSq_r - dpsiEffdr(r,L)/vr*dvPerSq_vr - L/r^2*dvPerSq_vt
    dvPar_r = dvPar_r - dpsiEffdr(r,L)/vr*dvPar_vr - L/r^2*dvPar_vt


    # println("dE_E = ",1.0/v*dvPar+v*dvPar_E+0.5*dvParSq_E+0.5*dvPerSq_E)
    # println("dE_L = ",v*dvPar_L+0.5*dvParSq_L+0.5*dvPerSq_L)
    # println((dvPerSq_E,dvPerSq_L))

    return dvPar_E, dvPar_L, dvParSq_E, dvParSq_L, dvPerSq_E, dvPerSq_L, dvPar_r, dvParSq_r, dvPerSq_r
end

function dvPar_EL_num(r::Float64, E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, m_test::Float64=m_field, eps::Float64=10^(-5))

    vr = vr_from_E_L(r,E,L)
    vt = vt_from_E_L(r,E,L)
    v2 = vSq_from_E_L(r,E,L)
    v = sqrt(v2)
    vr_v = vr/v
    vt_v = vt/v

    E_p = E + eps*abs(_E0)
    E_m = E - eps*abs(_E0)
    L_p = L + eps*_L0
    L_m = L - eps*_L0

    r_p = r + eps*_b
    r_m = r - eps*_b

    vr_Ep = vr_from_E_L(r,E_p,L)
    vt_Ep = vt_from_E_L(r,E_p,L)
    vr_Em = vr_from_E_L(r,E_m,L)
    vt_Em = vt_from_E_L(r,E_m,L)

    vr_Lp = vr_from_E_L(r,E,L_p)
    vt_Lp = vt_from_E_L(r,E,L_p)
    vr_Lm = vr_from_E_L(r,E,L_m)
    vt_Lm = vt_from_E_L(r,E,L_m)

    vr_rp = vr_from_E_L(r_p,E,L)
    vt_rp = vt_from_E_L(r_p,E,L)
    vr_rm = vr_from_E_L(r_m,E,L)
    vt_rm = vt_from_E_L(r_m,E,L)

    dvPar_Ep, dvParSq_Ep, dvPerSq_Ep, _ = localVelChangeIsoGrad(r,vr_Ep,vt_Ep,m_field,nbK,m_test)
    dvPar_Em, dvParSq_Em, dvPerSq_Em, _ = localVelChangeIsoGrad(r,vr_Em,vt_Em,m_field,nbK,m_test)
    dvPar_Lp, dvParSq_Lp, dvPerSq_Lp, _ = localVelChangeIsoGrad(r,vr_Lp,vt_Lp,m_field,nbK,m_test)
    dvPar_Lm, dvParSq_Lm, dvPerSq_Lm, _ = localVelChangeIsoGrad(r,vr_Lm,vt_Lm,m_field,nbK,m_test)

    dvPar_rp, dvParSq_rp, dvPerSq_rp, _ = localVelChangeIsoGrad(r_p,vr_rp,vt_rp,m_field,nbK,m_test)
    dvPar_rm, dvParSq_rm, dvPerSq_rm, _ = localVelChangeIsoGrad(r_m,vr_rm,vt_rm,m_field,nbK,m_test)

    # println((dvPar_Lp, dvParSq_Lp, dvPerSq_Lp))
    # println((dvPar_Lm, dvParSq_Lm, dvPerSq_Lm))

    # E-grad

    dvPar_E = (dvPar_Ep-dvPar_Em)/(2.0*eps*abs(_E0))
    dvParSq_E = (dvParSq_Ep-dvParSq_Em)/(2.0*eps*abs(_E0))
    dvPerSq_E = (dvPerSq_Ep-dvPerSq_Em)/(2.0*eps*abs(_E0))

    # L-grad
    dvPar_L = (dvPar_Lp-dvPar_Lm)/(2.0*eps*abs(_L0))
    dvParSq_L = (dvParSq_Lp-dvParSq_Lm)/(2.0*eps*abs(_L0))
    dvPerSq_L = (dvPerSq_Lp-dvPerSq_Lm)/(2.0*eps*abs(_L0))

    # r-grad
    dvPar_r = (dvPar_rp-dvPar_rm)/(2.0*eps*_b)
    dvParSq_r = (dvParSq_rp-dvParSq_rm)/(2.0*eps*_b)
    dvPerSq_r = (dvPerSq_rp-dvPerSq_rm)/(2.0*eps*_b)

    # println(dvParSq_L)
    # println(dvPerSq_L)

    return dvPar_E, dvPar_L, dvParSq_E, dvParSq_L, dvPerSq_E, dvPerSq_L, dvPar_r, dvParSq_r, dvPerSq_r
end

function dvPar_Iso_Compare(r::Float64, E::Float64, L::Float64, m_field::Float64,
            nbK::Int64=nbK_default, m_test::Float64=m_field, eps::Float64=10^(-5))

    dvPar_E_num, dvPar_L_num, dvParSq_E_num, dvParSq_L_num, dvPerSq_E_num, dvPerSq_L_num, dvPar_r_num, dvParSq_r_num, dvPerSq_r_num = dvPar_EL_num(r,E,L,m_field,nbK,m_test,eps)
    dvPar_E, dvPar_L, dvParSq_E, dvParSq_L, dvPerSq_E, dvPerSq_L, dvPar_r, dvParSq_r, dvPerSq_r = dvPar_EL(r,E,L,m_field,nbK,m_test)

    println("dvPar_E")
    println("Analytical = ",dvPar_E)
    println("Numerical  = ",dvPar_E_num)
    println("-------------------")
    println("dvPar_L")
    println("Analytical = ",dvPar_L)
    println("Numerical  = ",dvPar_L_num)
    println("-------------------")
    println("dvPar_r")
    println("Analytical = ",dvPar_r)
    println("Numerical  = ",dvPar_r_num)
    println("-------------------")
    println("dvParSq_E")
    println("Analytical = ",dvParSq_E)
    println("Numerical  = ",dvParSq_E_num)
    println("-------------------")
    println("dvParSq_L")
    println("Analytical = ",dvParSq_L)
    println("Numerical  = ",dvParSq_L_num)
    println("-------------------")
    println("dvParSq_r")
    println("Analytical = ",dvParSq_r)
    println("Numerical  = ",dvParSq_r_num)
    println("-------------------")
    println("dvPerSq_E")
    println("Analytical = ",dvPerSq_E)
    println("Numerical  = ",dvPerSq_E_num)
    println("-------------------")
    println("dvPerSq_L")
    println("Analytical = ",dvPerSq_L)
    println("Numerical  = ",dvPerSq_L_num)
    println("-------------------")
    println("dvPerSq_r")
    println("Analytical = ",dvPerSq_r)
    println("Numerical  = ",dvPerSq_r_num)
end


#############################################################
# Test dF/dvr and dF/dvt in integrands for dvPar
#############################################################

function dFdvrTest(r::Float64, vr::Float64, vt::Float64, w::Float64, varphi::Float64,
            phi::Float64, eps::Float64=10^(-5))

    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + v^2/2.0
    L = r*vt

    # analytical

    sinvarphi, cosvarphi = sincos(varphi)
    sinphi, cosphi = sincos(phi)

    w1 = w*cosvarphi
    w2 = w*sinvarphi*cosphi
    w3 = w*sinvarphi*sinphi

    Ep = E + w^2/2.0 - v*w1

    v1p = v-w1
    v2p = -w2
    v3p = -w3

    Lp = r * sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )

    Ftot, dFtotdE, dFtotdL = _FdF(Ep,Lp)

    dEpdvr = vr - vr_v*w*cosvarphi
    dEpdvt = vt - vt_v*w*cosvarphi
    dLpdvr = r^2*(vt^2/v^3*w*sinvarphi*sinphi+vr*vt/v^3*w*cosvarphi)*(vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp
    dLpdvt = r^2*(1.0-vr*vt/v^3*w*sinvarphi*sinphi-vr^2/v^3*w*cosvarphi)*(vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp

    dFdvr = dEpdvr*dFtotdE + dLpdvr*dFtotdL
    dFdvt = dEpdvt*dFtotdE + dLpdvt*dFtotdL





    # num vr

    vr_p = vr + eps*_v0
    vr_m = vr - eps*_v0
    v_p = sqrt(vr_p^2 + vt^2)
    v_m = sqrt(vr_m^2 + vt^2)
    E_p = psi(r) + v_p^2/2.0
    E_m = psi(r) + v_m^2/2.0

    Ep_p = E_p + w^2/2.0 - v_p*w1
    Ep_m = E_m + w^2/2.0 - v_m*w1

    v1p_p = v_p-w1
    v2p_p = -w2
    v3p_p = -w3
    v1p_m = v_m-w1
    v2p_m = -w2
    v3p_m = -w3

    Lp_p = r * sqrt(v2p_p^2 + (vr_p/v_p*v3p_p-vt/v_p*v1p_p)^2 )
    Lp_m = r * sqrt(v2p_m^2 + (vr_m/v_m*v3p_m-vt/v_m*v1p_m)^2 )

    Ftot_p, _ = _FdF(Ep_p,Lp_p)
    Ftot_m, _ = _FdF(Ep_m,Lp_m)

    dFdvr_num = (Ftot_p-Ftot_m)/(2.0*eps*_v0)

    println("dFdvt : ")
    println("Analytical = ",dFdvr)
    println("Numerical  = ",dFdvr_num)
end

#############################################################
# Test dLp/dvr and dLp/dvt in integrands for dvPar
#############################################################

function dLpdvrTest(r::Float64, vr::Float64, vt::Float64, w::Float64, varphi::Float64,
            phi::Float64, eps::Float64=10^(-5))



    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + v^2/2.0
    L = r*vt

    # analytical

    sinvarphi, cosvarphi = sincos(varphi)
    sinphi, cosphi = sincos(phi)

    w1 = w*cosvarphi
    w2 = w*sinvarphi*cosphi
    w3 = w*sinvarphi*sinphi

    Ep = E + w^2/2.0 - v*w1

    v1p = v-w1
    v2p = -w2
    v3p = -w3

    Lp = r * sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )

    Ftot, dFtotdE, dFtotdL = _FdF(Ep,Lp)

    dEpdvr = vr - vr_v*w*cosvarphi
    dEpdvt = vt - vt_v*w*cosvarphi
    dLpdvr = r^2*(vt^2/v^3*w*sinvarphi*sinphi+vr*vt/v^3*w*cosvarphi)*(vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp
    dLpdvt = r^2*(1.0-vr*vt/v^3*w*sinvarphi*sinphi-vr^2/v^3*w*cosvarphi)*(vt+vr_v*w*sinvarphi*sinphi-vt_v*w*cosvarphi)/Lp

    # numerical

    vr_p = vr + eps*_v0
    vr_m = vr - eps*_v0
    v_p = sqrt(vr_p^2 + vt^2)
    v_m = sqrt(vr_m^2 + vt^2)
    E_p = psi(r) + v_p^2/2.0
    E_m = psi(r) + v_m^2/2.0

    Ep_p = E_p + w^2/2.0 - v_p*w1
    Ep_m = E_m + w^2/2.0 - v_m*w1

    v1p_p = v_p-w1
    v2p_p = -w2
    v3p_p = -w3
    v1p_m = v_m-w1
    v2p_m = -w2
    v3p_m = -w3

    Lp_p = r * sqrt(v2p_p^2 + (vr_p/v_p*v3p_p-vt/v_p*v1p_p)^2 )
    Lp_m = r * sqrt(v2p_m^2 + (vr_m/v_m*v3p_m-vt/v_m*v1p_m)^2 )

    dLpdvr_num = (Lp_p-Lp_m)/(2.0*eps*_v0)

    println("dLpdvt : ")
    println("Analytical = ",dLpdvr)
    println("Numerical  = ",dLpdvr_num)
end
