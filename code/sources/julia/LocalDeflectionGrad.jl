


# Compute vr- and vt-gradients of local velocity deflections
# No rotation

function localVelChange3DGrad(r::Float64, vr::Float64, vt::Float64,
                        m_field::Float64,
                        nbK::Int64=nbK_default, m_test::Float64=m_field)


    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + v^2/2.0
    L = r*vt





    # x = r*costheta*cosI
    # y = r*sintheta

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

        dvPar_vr_phi = 0.0
        dvParSq_vr_phi = 0.0
        dvPerSq_vr_phi = 0.0
        dvPar_vt_phi = 0.0
        dvParSq_vt_phi = 0.0
        dvPerSq_vt_phi = 0.0

        for iphi=1:nbK
            phi = 2.0*pi*(iphi-0.5)/nbK
            sinphi, cosphi = sincos(phi)

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
                # Lp = r * sqrt(v^2+w^2-2.0*v*w*cos(varphi)-(vr/v*(v-w1)+vt/v*(-w2))^2)

                v1p = v-w1
                v2p = -w2
                v3p = -w3

                Lp = r * sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )


                # wx = cosI*vX/v*w1 - sinI*w2 + cosI*vY/v
                # wy = vY/v*w1 - vX/v*w3
                # Lzp = Lz + y*wx - x*wy




                F, dFdEp, dFdLp = _FdF(Ep,Lp)

                dEpdvr = vr - vr/v*w*cosvarphi
                dEpdvt = vt - vt/v*w*cosvarphi
                dLpdvr = r^2*w*sinvarphi*sinphi/v * (vt+vr/v*w*sinvarphi*sinphi-vt/v*w*cosvarphi)/Lp
                dLpdvt = r^2  * (vt+vr/v*w*sinvarphi*sinphi-vt/v*w*cosvarphi)/Lp

                dFdvr = dEpdvr*dFdEp + dLpdvr*dFdLp
                dFdvt = dEpdvt*dFdEp + dLpdvt*dFdLp

                dvPar_vr_w += dFdvr
                dvPar_vt_w += dFdvt

                dvParSq_vr_w += F
                dvPerSq_vr_w += F
                dvParSq_vt_w += F
                dvPerSq_vt_w += F

            end

            dvPar_vr_w *= wmax
            dvPar_vt_w *= wmax

            dvParSq_vr_w *= wmax
            dvPerSq_vr_w *= wmax
            dvParSq_vt_w *= wmax
            dvPerSq_vt_w *= wmax

            dvPar_vr_phi += dvPar_vr_w
            dvPar_vt_phi += dvPar_vt_w

            dvParSq_vr_phi += dvParSq_vr_w
            dvPerSq_vr_phi += dvPerSq_vr_w
            dvParSq_vt_phi += dvParSq_vt_w
            dvPerSq_vt_phi += dvPerSq_vt_w

        end

        dvPar_vr += 2.0*sinvarphi*cosvarphi*dvPar_vr_phi
        dvPar_vt += 2.0*sinvarphi*cosvarphi*dvPar_vt_phi

        dvParSq_vr += -3.0*vr/v*sinvarphi^3*cosvarphi* dvParSq_vr_phi
        dvPerSq_vr += vr/v*sinvarphi*cosvarphi*(3.0*sinvarphi^2-2.0)*dvPerSq_vr_phi
        dvParSq_vt += -3.0*vt/v*sinvarphi^3*cosvarphi* dvParSq_vt_phi
        dvPerSq_vt += vt/v*sinvarphi*cosvarphi*(3.0*sinvarphi^2-2.0)*dvPerSq_vt_phi


    end

    # dvPar *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvParSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvPerSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)

    pref = 2.0*pi^2/nbK^3
    pref *= 2.0*pi*_G^2*logCoulomb

    dvPar_vr *= -pref * (m_field + m_test)
    dvPar_vt *= -pref * (m_field + m_test)
    dvParSq_vr *= 2.0*pref * m_field
    dvParSq_vt *= 2.0*pref * m_field
    dvPerSq_vr *= 2.0*pref * m_field
    dvPerSq_vt *= 2.0*pref * m_field

    return dvPar_vr, dvPar_vt, dvParSq_vr, dvParSq_vt, dvPerSq_vr, dvPerSq_vt
end


function localVelChange3DGradNum(r::Float64, vr::Float64, vt::Float64,
                        m_field::Float64, eps::Float64=10^(-5),
                        nbK::Int64=nbK_default, m_test::Float64=m_field)


    vrp = vr + eps*_v0
    vrm = vr - eps*_v0
    vtp = vt + eps*_v0
    vtm = vt - eps*_v0


    dvPar_rp, dvPar2_rp, dvPerp2_rp = localVelChange3D(r,0.0,vrp,vt,0.0,m_field,0.0,nbK,m_test)
    dvPar_rm, dvPar2_rm, dvPerp2_rm = localVelChange3D(r,0.0,vrm,vt,0.0,m_field,0.0,nbK,m_test)
    dvPar_tp, dvPar2_tp, dvPerp2_tp = localVelChange3D(r,0.0,vr,vtp,0.0,m_field,0.0,nbK,m_test)
    dvPar_tm, dvPar2_tm, dvPerp2_tm = localVelChange3D(r,0.0,vr,vtm,0.0,m_field,0.0,nbK,m_test)

    dvPar_vr = (dvPar_rp-dvPar_rm)/(2.0*eps*_v0)
    dvPar_vt = (dvPar_tp-dvPar_tm)/(2.0*eps*_v0)

    dvParSq_vr = (dvPar2_rp-dvPar2_rm)/(2.0*eps*_v0)
    dvParSq_vt = (dvPar2_tp-dvPar2_tm)/(2.0*eps*_v0)
    dvPerSq_vr = (dvPerp2_rp-dvPerp2_rm)/(2.0*eps*_v0)
    dvPerSq_vt = (dvPerp2_tp-dvPerp2_tm)/(2.0*eps*_v0)




    return dvPar_vr, dvPar_vt, dvParSq_vr, dvParSq_vt, dvPerSq_vr, dvPerSq_vt
end
