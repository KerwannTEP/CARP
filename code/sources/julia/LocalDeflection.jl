##################################################
# Computation of the local velocity deflections
##################################################


function localVelChange3D(r::Float64, theta::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbK::Int64=nbK_default, m_test::Float64=m_field)


    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + 0.5*vSq
    L = r*vt
    Lz = L*cosI
    sinI = sqrt(abs(1.0 - cosI^2))

    sintheta, costheta = sincos(theta)

    # x = r*costheta*cosI
    # y = r*sintheta

    dvPar = 0.0
    dvParSq = 0.0
    dvPerSq = 0.0

    for ivarphi=1:nbK
        varphi = pi*(ivarphi-0.5)/nbK
        sinvarphi, cosvarphi = sincos(varphi)

        wmax = v*cosvarphi + sqrt(vSq*cosvarphi^2 - 2.0*E)

        dvPar_phi = 0.0
        # dvParSq_phi = 0.0
        # dvPerSq_phi = 0.0

        dvSq_phi = 0.0

        for iphi=1:nbK
            phi = 2.0*pi*(iphi-0.5)/nbK
            sinphi, cosphi = sincos(phi)

            dvPar_w = 0.0
            # dvParSq_w = 0.0
            # dvPerSq_w = 0.0

            dvSq_w = 0.0



            for iw=1:nbK
                w = wmax*(iw-0.5)/nbK

                w1 = w*cosvarphi
                w2 = w*sinvarphi*cosphi
                w3 = w*sinvarphi*sinphi

                Ep = E + 0.5*w^2 - v*w1
                # Lp = r * sqrt(v^2+w^2-2.0*v*w*cos(varphi)-(vr/v*(v-w1)+vt/v*(-w2))^2)

                v1p = v-w1
                v2p = -w2
                v3p = -w3

                Lp = r * sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )


                # wx = cosI*vX/v*w1 - sinI*w2 + cosI*vY/v
                # wy = vY/v*w1 - vX/v*w3
                # Lzp = Lz + y*wx - x*wy

                Lzp = r*(v1p*vt_v-v3p*vr_v)*cosI + r*v2p*sintheta*sinI


                Frot = _Frot(Ep,Lp,Lzp,alpha)

                dvPar_w += Frot
                # dvParSq_w += w*Frot
                # dvPerSq_w += w*Frot

                dvSq_w += w*Frot

            end

            dvPar_w *= wmax
            # dvParSq_w *= wmax
            # dvPerSq_w *= wmax

            dvSq_w *= wmax

            dvPar_phi += dvPar_w
            # dvParSq_phi += dvParSq_w
            # dvPerSq_phi += dvPerSq_w

            # dvParSq_phi += dvSq_w
            # dvPerSq_phi += dvSq_w

            dvSq_phi += dvSq_w

        end



        dvPar += 2.0*sinvarphi*cosvarphi*dvPar_phi
        # dvParSq += sinvarphi^3*dvParSq_phi
        # dvPerSq += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_phi

        dvParSq += sinvarphi^3*dvSq_phi
        dvPerSq += sinvarphi*(1.0+cosvarphi^2)*dvSq_phi

    end

    # dvPar *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvParSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvPerSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)

    pref = 2.0*pi^2/nbK^3
    pref *= 2.0*pi*_G^2*logCoulomb

    dvPar *= -pref * (m_field + m_test)
    dvParSq *= 2.0*pref * m_field
    dvPerSq *= 2.0*pref * m_field

    return dvPar, dvParSq, dvPerSq
end




# function localOrbitalChange(r::Float64, theta::Float64, vr::Float64, vt::Float64,
#                         cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
#                         nbK::Int64=nbK_default, m_test::Float64=m_field)
#
#     v = sqrt(vr^2 + vt^2)
#     E = psi(r) + v^2/2.0
#     L = r*vt
#     Lz = L*cosI
#     vr_v = vr/v
#     vt_v = vt/v
#
#     #dvPar, dvParSq, dvPerpSq = localVelChange3DNew(r,theta,vr,vt,cosI,m_field,alpha,nbK,m_test)
#     dvPar, dvParSq, dvPerpSq = localVelChange3D(r,theta,vr,vt,cosI,m_field,alpha,nbK,m_test)
#
#     termSq = vt_v^2*dvParSq + 0.5*vr_v^2*dvPerpSq
#
#
#     dE = v*dvPar + 0.5*dvParSq + 0.5*dvPerpSq
#     dE2 = v^2*dvParSq
#     #dL = L/v*dvPar + r^2/(4.0*L)*dvPerpSq
#     dL = r*(vt_v*dvPar + 0.25/vt*dvPerpSq)
#     #dL2 = L^2/v^2*dvParSq + 0.5*(r^2-L^2/v^2)*dvPerpSq
#     dL2 = r^2*(termSq)
#     dEdL = L*dvParSq
#
#     dLz = Lz/v*dvPar
#     #dLz2 = (Lz/L)^2*(L^2/v^2*dvParSq + 0.5*(r^2-L^2/v^2)*dvPerpSq) + 0.5*r^2*sin(theta)^2*(1.0-(Lz/L)^2)*dvPerpSq
#     #dLz2 = (cosI)^2*(L^2/v^2*dvParSq + 0.5*(r^2*vr^2/v^2)*dvPerpSq) + 0.5*r^2*sin(theta)^2*(1.0-(cosI)^2)*dvPerpSq
#     dLz2 = r^2*((cosI)^2*(termSq) + 0.5*sin(theta)^2*(1.0-(cosI)^2)*dvPerpSq)
#
#     dEdLz = Lz*dvParSq
#     #dLdLz = Lz/L*(L^2/v^2*dvParSq + 0.5*(r^2-L^2/v^2)*dvPerpSq)
#     dLdLz = r^2*cosI*(termSq)
#
#     return dE, dL, dLz, dE2, dL2, dLz2, dEdL, dEdLz, dLdLz
# end
