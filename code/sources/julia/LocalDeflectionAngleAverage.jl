
# return int dtheta/2pi (dvPar, dvParSq, dvPerpSq, sin^2(theta)*dvPerpSq)

function localVelChange3DAngleAverage(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbAvrTh::Int64=nbAvrTh_default, nbK::Int64=nbK_default,
                        m_test::Float64=m_field)


    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + 0.5*vSq
    L = r*vt
    Lz = L*cosI
    sinI = sqrt(1.0 - cosI^2)



    # x = r*costheta*cosI
    # y = r*sintheta

    dvPar = 0.0
    dvParSq = 0.0
    dvPerSq = 0.0
    sinSqdvPerSq = 0.0

    for ivarphi=1:nbK
        varphi = pi*(ivarphi-0.5)/nbK
        sinvarphi, cosvarphi = sincos(varphi)

        wmax = v*cosvarphi + sqrt(vSq*cosvarphi^2 - 2.0*E)

        dvPar_phi = 0.0
        dvParSq_phi = 0.0
        dvPerSq_phi = 0.0
        sinSqdvPerSq_phi = 0.0

        for iphi=1:nbK
            phi = 2.0*pi*(iphi-0.5)/nbK
            sinphi, cosphi = sincos(phi)

            dvPar_w = 0.0
            dvParSq_w = 0.0
            dvPerSq_w = 0.0
            sinSqdvPerSq_w = 0.0



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

                Ftot = _F(Ep,Lp)

                sum_g = 0.0
                sum_sinSqg = 0.0

                for ith=1:nbAvrTh
                    theta = 2*pi/nbAvrTh * (ith-0.5)

                    sintheta = sin(theta)
                    Lzp = r*(v1p*vt_v-v3p*vr_v)*cosI + r*v2p*sintheta*sinI
                    gRot = _gRot(Lzp/Lp)
                    sum_g += gRot
                    sum_sinSqg += sintheta^2*gRot


                end

                dvPar_w += Ftot * (1.0 + alpha*sum_g/nbAvrTh)
                dvParSq_w += w*Ftot * (1.0 + alpha*sum_g/nbAvrTh)
                dvPerSq_w += w*Ftot * (1.0 + alpha*sum_g/nbAvrTh)
                sinSqdvPerSq_w += w*Ftot * (0.5 + alpha*sum_sinSqg/nbAvrTh)

            end

            dvPar_w *= wmax
            dvParSq_w *= wmax
            dvPerSq_w *= wmax
            sinSqdvPerSq_w *= wmax



            dvPar_phi += dvPar_w
            dvParSq_phi += dvParSq_w
            dvPerSq_phi += dvPerSq_w
            sinSqdvPerSq_phi += sinSqdvPerSq_w


        end



        dvPar += 2.0*sinvarphi*cosvarphi*dvPar_phi
        dvParSq += sinvarphi^3*dvParSq_phi
        dvPerSq += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_phi
        sinSqdvPerSq += sinvarphi*(1.0+cosvarphi^2)*sinSqdvPerSq_phi


    end

    # dvPar *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvParSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvPerSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)

    pref = 2.0*pi^2/nbK^3
    pref *= 2.0*pi*_G^2*logCoulomb

    dvPar *= -pref * (m_field + m_test)
    dvParSq *= 2.0*pref * m_field
    dvPerSq *= 2.0*pref * m_field
    sinSqdvPerSq *= 2.0*pref * m_field

    return dvPar, dvParSq, dvPerSq, sinSqdvPerSq
end



function localVelChange3DAngleAverage_nbKseparate(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbAvrTh::Int64=nbAvrTh_default, nbw::Int64=nbw_default,
                        nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                        m_test::Float64=m_field)


    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + 0.5*vSq
    L = r*vt
    Lz = L*cosI
    sinI = sqrt(1.0 - cosI^2)



    # x = r*costheta*cosI
    # y = r*sintheta

    dvPar = 0.0
    dvParSq = 0.0
    dvPerSq = 0.0
    sinSqdvPerSq = 0.0

    for ivarphi=1:nbw
        varphi = pi*(ivarphi-0.5)/nbw
        sinvarphi, cosvarphi = sincos(varphi)

        wmax = v*cosvarphi + sqrt(vSq*cosvarphi^2 - 2.0*E)

        dvPar_phi = 0.0
        dvParSq_phi = 0.0
        dvPerSq_phi = 0.0
        sinSqdvPerSq_phi = 0.0

        for iphi=1:nbphi
            phi = 2.0*pi*(iphi-0.5)/nbphi
            sinphi, cosphi = sincos(phi)

            dvPar_w = 0.0
            dvParSq_w = 0.0
            dvPerSq_w = 0.0
            sinSqdvPerSq_w = 0.0



            for iw=1:nbvarphi
                w = wmax*(iw-0.5)/nbvarphi

                w1 = w*cosvarphi
                w2 = w*sinvarphi*cosphi
                w3 = w*sinvarphi*sinphi

                Ep = E + 0.5*w^2 - v*w1
                # Lp = r * sqrt(v^2+w^2-2.0*v*w*cos(varphi)-(vr/v*(v-w1)+vt/v*(-w2))^2)

                v1p = v-w1
                v2p = -w2
                v3p = -w3

                Lp = r * sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )

                Ftot = _F(Ep,Lp)

                sum_g = 0.0
                sum_sinSqg = 0.0

                for ith=1:nbAvrTh
                    theta = 2*pi/nbAvrTh * (ith-0.5)

                    sintheta = sin(theta)
                    Lzp = r*(v1p*vt_v-v3p*vr_v)*cosI + r*v2p*sintheta*sinI
                    gRot = _gRot(Lzp/Lp)
                    sum_g += gRot
                    sum_sinSqg += sintheta^2*gRot


                end

                dvPar_w += Ftot * (1.0 + alpha*sum_g/nbAvrTh)
                dvParSq_w += w*Ftot * (1.0 + alpha*sum_g/nbAvrTh)
                dvPerSq_w += w*Ftot * (1.0 + alpha*sum_g/nbAvrTh)
                sinSqdvPerSq_w += w*Ftot * (0.5 + alpha*sum_sinSqg/nbAvrTh)

            end

            dvPar_w *= wmax
            dvParSq_w *= wmax
            dvPerSq_w *= wmax
            sinSqdvPerSq_w *= wmax



            dvPar_phi += dvPar_w
            dvParSq_phi += dvParSq_w
            dvPerSq_phi += dvPerSq_w
            sinSqdvPerSq_phi += sinSqdvPerSq_w


        end



        dvPar += 2.0*sinvarphi*cosvarphi*dvPar_phi
        dvParSq += sinvarphi^3*dvParSq_phi
        dvPerSq += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_phi
        sinSqdvPerSq += sinvarphi*(1.0+cosvarphi^2)*sinSqdvPerSq_phi


    end

    # dvPar *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvParSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvPerSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)

    pref = 2.0*pi^2/(nbw*nbphi*nbvarphi)
    pref *= 2.0*pi*_G^2*logCoulomb

    dvPar *= -pref * (m_field + m_test)
    dvParSq *= 2.0*pref * m_field
    dvPerSq *= 2.0*pref * m_field
    sinSqdvPerSq *= 2.0*pref * m_field

    return dvPar, dvParSq, dvPerSq, sinSqdvPerSq
end

# Exact theta-average for g=sign
function localVelChange3DAngleAverageExact(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbK::Int64=nbK_default,
                        m_test::Float64=m_field)


    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + 0.5*vSq
    L = r*vt
    Lz = L*cosI
    sinI = sqrt(abs(1.0 - cosI^2))



    # x = r*costheta*cosI
    # y = r*sintheta

    dvPar = 0.0
    dvParSq = 0.0
    dvPerSq = 0.0
    sinSqdvPerSq = 0.0

    for ivarphi=1:nbK
        varphi = pi*(ivarphi-0.5)/nbK
        sinvarphi, cosvarphi = sincos(varphi)

        wmax = v*cosvarphi + sqrt(abs(vSq*cosvarphi^2 - 2.0*E))

        dvPar_phi = 0.0
        dvParSq_phi = 0.0
        dvPerSq_phi = 0.0
        sinSqdvPerSq_phi = 0.0

        for iphi=1:nbK
            phi = 2.0*pi*(iphi-0.5)/nbK
            sinphi, cosphi = sincos(phi)

            dvPar_w = 0.0
            dvParSq_w = 0.0
            dvPerSq_w = 0.0
            sinSqdvPerSq_w = 0.0



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

                Ftot = _F(Ep,Lp)

                nu = -(v1p*vt_v-v3p*vr_v)*cosI

                sum_g =  0.0
                sum_sinSqg = 0.0

                if (v2p == 0.0)
                    sum_g = sign(-nu)
                    sum_sinSqg = 0.5*sign(-nu)
                else
                    mu = nu/(v2p*sinI)
                    if (abs(mu) >= 1.0)
                        sum_g = sign(-nu)
                        sum_sinSqg = sign(-nu)/2.0
                    else
                        sum_g = -2.0*asin(nu/(abs(v2p)*sinI))/pi
                        sum_sinSqg = 0.5*(-2.0*asin(nu/(abs(v2p)*sinI))/pi + 2.0/pi*(nu/(abs(v2p)*sinI))*sqrt(1.0-(nu/(abs(v2p)*sinI))^2))
                    end
                end



                dvPar_w += Ftot * (1.0 + alpha*sum_g)
                dvParSq_w += w*Ftot * (1.0 + alpha*sum_g)
                dvPerSq_w += w*Ftot * (1.0 + alpha*sum_g)
                sinSqdvPerSq_w += w*Ftot * (0.5 + alpha*sum_sinSqg)

            end

            dvPar_w *= wmax
            dvParSq_w *= wmax
            dvPerSq_w *= wmax
            sinSqdvPerSq_w *= wmax



            dvPar_phi += dvPar_w
            dvParSq_phi += dvParSq_w
            dvPerSq_phi += dvPerSq_w
            sinSqdvPerSq_phi += sinSqdvPerSq_w


        end



        dvPar += 2.0*sinvarphi*cosvarphi*dvPar_phi
        dvParSq += sinvarphi^3*dvParSq_phi
        dvPerSq += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_phi
        sinSqdvPerSq += sinvarphi*(1.0+cosvarphi^2)*sinSqdvPerSq_phi


    end

    # dvPar *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvParSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvPerSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)

    pref = 2.0*pi^2/nbK^3
    pref *= 2.0*pi*_G^2*logCoulomb

    dvPar *= -pref * (m_field + m_test)
    dvParSq *= 2.0*pref * m_field
    dvPerSq *= 2.0*pref * m_field
    sinSqdvPerSq *= 2.0*pref * m_field

    return dvPar, dvParSq, dvPerSq, sinSqdvPerSq
end

function localVelChange3DAngleAverageExact_nbKseparate(r::Float64, vr::Float64,
                        vt::Float64, cosI::Float64, m_field::Float64,
                        alpha::Float64=alphaRot, nbw::Int64=nbw_default,
                        nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                        m_test::Float64=m_field)


    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + 0.5*vSq
    L = r*vt
    Lz = L*cosI
    sinI = sqrt(abs(1.0 - cosI^2))



    # x = r*costheta*cosI
    # y = r*sintheta

    dvPar = 0.0
    dvParSq = 0.0
    dvPerSq = 0.0
    sinSqdvPerSq = 0.0

    for ivarphi=1:nbvarphi
        varphi = pi*(ivarphi-0.5)/nbvarphi
        sinvarphi, cosvarphi = sincos(varphi)

        wmax = v*cosvarphi + sqrt(abs(vSq*cosvarphi^2 - 2.0*E))

        dvPar_phi = 0.0
        dvParSq_phi = 0.0
        dvPerSq_phi = 0.0
        sinSqdvPerSq_phi = 0.0

        for iphi=1:nbphi
            phi = 2.0*pi*(iphi-0.5)/nbphi
            sinphi, cosphi = sincos(phi)

            dvPar_w = 0.0
            dvParSq_w = 0.0
            dvPerSq_w = 0.0
            sinSqdvPerSq_w = 0.0



            for iw=1:nbw
                w = wmax*(iw-0.5)/nbw

                w1 = w*cosvarphi
                w2 = w*sinvarphi*cosphi
                w3 = w*sinvarphi*sinphi

                Ep = E + 0.5*w^2 - v*w1
                # Lp = r * sqrt(v^2+w^2-2.0*v*w*cos(varphi)-(vr/v*(v-w1)+vt/v*(-w2))^2)

                v1p = v-w1
                v2p = -w2
                v3p = -w3

                Lp = r * sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )

                Ftot = _F(Ep,Lp)

                nu = -(v1p*vt_v-v3p*vr_v)*cosI

                sum_g =  0.0
                sum_sinSqg = 0.0

                if (v2p == 0.0)
                    sum_g = sign(-nu)
                    sum_sinSqg = 0.5*sign(-nu)
                else
                    mu = nu/(v2p*sinI)
                    if (abs(mu) >= 1.0)
                        sum_g = sign(-nu)
                        sum_sinSqg = sign(-nu)/2.0
                    else
                        sum_g = -2.0*asin(nu/(abs(v2p)*sinI))/pi
                        sum_sinSqg = 0.5*(-2.0*asin(nu/(abs(v2p)*sinI))/pi + 2.0/pi*(nu/(abs(v2p)*sinI))*sqrt(1.0-(nu/(abs(v2p)*sinI))^2))
                    end
                end



                dvPar_w += Ftot * (1.0 + alpha*sum_g)
                dvParSq_w += w*Ftot * (1.0 + alpha*sum_g)
                dvPerSq_w += w*Ftot * (1.0 + alpha*sum_g)
                sinSqdvPerSq_w += w*Ftot * (0.5 + alpha*sum_sinSqg)

            end

            dvPar_w *= wmax
            dvParSq_w *= wmax
            dvPerSq_w *= wmax
            sinSqdvPerSq_w *= wmax



            dvPar_phi += dvPar_w
            dvParSq_phi += dvParSq_w
            dvPerSq_phi += dvPerSq_w
            sinSqdvPerSq_phi += sinSqdvPerSq_w


        end



        dvPar += 2.0*sinvarphi*cosvarphi*dvPar_phi
        dvParSq += sinvarphi^3*dvParSq_phi
        dvPerSq += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_phi
        sinSqdvPerSq += sinvarphi*(1.0+cosvarphi^2)*sinSqdvPerSq_phi


    end

    # dvPar *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvParSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvPerSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)

    pref = 2.0*pi^2/(nbw*nbvarphi*nbphi)
    pref *= 2.0*pi*_G^2*logCoulomb

    dvPar *= -pref * (m_field + m_test)
    dvParSq *= 2.0*pref * m_field
    dvPerSq *= 2.0*pref * m_field
    sinSqdvPerSq *= 2.0*pref * m_field

    return dvPar, dvParSq, dvPerSq, sinSqdvPerSq
end




# cosI-drift coeff
# analytical
# w-integral becomes zeta-integral
function localDriftCosIAngleAverageExact_nbKseparate(r::Float64, vr::Float64,
                        vt::Float64, cosI::Float64, m_field::Float64,
                        alpha::Float64=alphaRot, nbzeta::Int64=nbw_default,
                        nbvarphi::Int64=nbvarphi_default, nbphi::Int64=nbphi_default,
                        m_test::Float64=m_field)


    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + 0.5*vSq
    L = r*vt
    Lz = L*cosI
    sinI = sqrt(abs(1.0 - cosI^2))

    drift = 0.0


    for ivarphi=1:nbvarphi
        varphi = pi*(ivarphi-0.5)/nbvarphi
        sinvarphi, cosvarphi = sincos(varphi)

        wmax = v*cosvarphi + sqrt(abs(vSq*cosvarphi^2 - 2.0*E))

        drift_phi = 0.0

        for iphi=1:nbphi
            phi = 2.0*pi*(iphi-0.5)/nbphi
            sinphi, cosphi = sincos(phi)

            drift_zeta = 0.0

            xwmax = -(vt-vt_v*wmax*cosvarphi+vr_v*wmax*sinvarphi*sinphi)*cosI/(wmax*sinI*abs(sinvarphi*cosphi))
            #wmin = -vt*cosI/(vr_v*sinvarphi*sinphi*cosI-vt_v*cosvarphi*cosI-sign(cosI)*sinI*abs(sinvarphi*cosphi))
            xinf = (vt_v*cosvarphi-vr_v*sinvarphi*sinphi)*cosI/(sinI*abs(sinvarphi*cosphi))
            # x(w) must be between -1 and 1
            # wmin = w(x=-1)



            # cos I > 0
            # x(0) = -infty
            # w(x) is increasing, i,e, x(w) is increasing
            if (cosI > 0.0)

                vanish_int = true # integrand vanishes

                if (-1.0 < xinf <= 1.0)
                    if (xwmax > -1.0)
                        xmin = -1.0
                        xmax = xwmax
                        vanish_int = false
                    end
                elseif (xinf > 1.0)
                    if (-1.0 < xwmax <= 1.0)
                        xmin = -1.0
                        xmax = xwmax
                        vanish_int = false
                    elseif (xwmax > 1.0)
                        xmin = -1.0
                        xmax = 1.0
                        vanish_int = false
                    end
                end

                if (vanish_int == false)
                    zetamin = acos(xmax)
                    zetamax = acos(xmin)

                    #println(xinf)

                    for izeta=1:nbzeta
                        zeta = zetamin + (zetamax-zetamin)/nbzeta*(izeta-0.5)
                        sinzeta, coszeta = sincos(zeta)
                        x = cos(zeta)
                        w = -vt*cosI/(vr_v*sinvarphi*sinphi*cosI-vt_v*cosvarphi*cosI+sinI*abs(sinvarphi*cosphi)*x)

                        w1 = w*cosvarphi
                        w2 = w*sinvarphi*cosphi
                        w3 = w*sinvarphi*sinphi

                        Ep = E + 0.5*w^2 - v*w1

                        v1p = v-w1
                        v2p = -w2
                        v3p = -w3

                        Lp = r * sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )

                        Ftot = _F(Ep,Lp)

                        den = vr_v*sinvarphi*sinphi*cosI - vt_v*cosvarphi*cosI+x*sinI*abs(sinvarphi*cosphi)

                        num1 = vt*coszeta*sinzeta^2*abs(cosphi*sinvarphi)*cosI^2*sinI
                        num2 = coszeta^2*(vt-vt_v*w*cosvarphi+vr_v*w*sinphi*sinvarphi)

                        drift_zeta += Ftot/den * (num1/den^2 + num2/den)*(zetamax-zetamin)/nbzeta
                    end

                end
            elseif (cosI < 0.0)
                # cos I < 0
                # x(0) = +infty
                # w(x) is decreasing, i,e, x(w) is decreasing
                vanish_int = true # integrand vanishes

                if (-1.0 < xinf <= 1.0)
                    if (xwmax < 1.0)
                        xmin = xwmax
                        xmax = 1.0
                        vanish_int = false
                    end
                elseif (xinf <= -1.0)
                    if (-1.0 < xwmax <= 1.0)
                        xmin = xwmax
                        xmax = 1.0
                        vanish_int = false
                    elseif (xwmax <= -1.0)
                        xmin = -1.0
                        xmax = 1.0
                        vanish_int = false
                    end
                end

                if (vanish_int == false)
                    zetamin = acos(xmax)
                    zetamax = acos(xmin)

                    for izeta=1:nbzeta
                        zeta = zetamin + (zetamax-zetamin)/nbzeta*(izeta-0.5)
                        sinzeta, coszeta = sincos(zeta)
                        x = cos(zeta)
                        w = -vt*cosI/(vr_v*sinvarphi*sinphi*cosI-vt_v*cosvarphi*cosI+sinI*abs(sinvarphi*cosphi)*x)

                        w1 = w*cosvarphi
                        w2 = w*sinvarphi*cosphi
                        w3 = w*sinvarphi*sinphi

                        Ep = E + 0.5*w^2 - v*w1

                        v1p = v-w1
                        v2p = -w2
                        v3p = -w3

                        Lp = r * sqrt(v2p^2 + (vr_v*v3p-vt_v*v1p)^2 )

                        Ftot = _F(Ep,Lp)

                        den = vr_v*sinvarphi*sinphi*cosI - vt_v*cosvarphi*cosI+x*sinI*abs(sinvarphi*cosphi)

                        num1 = vt*coszeta*sinzeta^2*abs(cosphi*sinvarphi)*cosI^2*sinI
                        num2 = coszeta^2*(vt-vt_v*w*cosvarphi+vr_v*w*sinphi*sinvarphi)

                        drift_zeta += Ftot/den * (num1/den^2 + num2/den)*(zetamax-zetamin)/nbzeta

                    end

                end

            end

            drift_phi += drift_zeta

        end

        drift += drift_phi*sinvarphi*(1.0+cosvarphi^2)
    end

    pref = 2.0*pi^2/(nbvarphi*nbphi)
    pref *= 2.0*pi*_G^2*logCoulomb

    drift *= 2.0*pref * m_field *(-r*alpha*abs(cosI))/(2.0*L*pi)

    return drift
end


















# Exact theta-average for g=sign
function localVelChange3DAngleAverageExactTest(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbK::Int64=nbK_default,
                        m_test::Float64=m_field)


    vSq =  vr^2 + vt^2
    v = sqrt(vSq)
    vr_v = vr/v
    vt_v = vt/v

    E = psi(r) + 0.5*vSq
    L = r*vt
    Lz = L*cosI
    sinI = sqrt(abs(1.0 - cosI^2))



    # x = r*costheta*cosI
    # y = r*sintheta

    dvPar = 0.0
    dvParSq = 0.0
    dvPerSq = 0.0
    sinSqdvPerSq = 0.0
    cosdvPerSq = 0.0

    for ivarphi=1:nbK
        varphi = pi*(ivarphi-0.5)/nbK
        sinvarphi, cosvarphi = sincos(varphi)

        wmax = v*cosvarphi + sqrt(abs(vSq*cosvarphi^2 - 2.0*E))

        dvPar_phi = 0.0
        dvParSq_phi = 0.0
        dvPerSq_phi = 0.0
        sinSqdvPerSq_phi = 0.0
        cosdvPerSq_phi = 0.0

        for iphi=1:nbK
            phi = 2.0*pi*(iphi-0.5)/nbK
            sinphi, cosphi = sincos(phi)

            dvPar_w = 0.0
            dvParSq_w = 0.0
            dvPerSq_w = 0.0
            sinSqdvPerSq_w = 0.0
            cosdvPerSq_w = 0.0



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

                Ftot = _F(Ep,Lp)

                nu = -(v1p*vt_v-v3p*vr_v)*cosI

                sum_g =  0.0
                sum_sinSqg = 0.0
                sum_cosg = 0.0

                if (v2p == 0.0)
                    sum_g = sign(-nu)
                    sum_sinSqg = 0.5*sign(-nu)
                    sum_cosg = 0.0
                else
                    mu = nu/(v2p*sinI)
                    if (abs(mu) >= 1.0)
                        sum_g = sign(-nu)
                        sum_sinSqg = sign(-nu)/2.0
                        sum_cosg = 0.0
                    else
                        sum_g = -2.0*asin(nu/(abs(v2p)*sinI))/pi
                        sum_sinSqg = 0.5*(-2.0*asin(nu/(abs(v2p)*sinI))/pi + 2.0/pi*(nu/(abs(v2p)*sinI))*sqrt(1.0-(nu/(abs(v2p)*sinI))^2))
                        sum_cosg = 0.5*(2.0/pi*(nu/(abs(v2p)*sinI))*sqrt(1.0-(nu/(abs(v2p)*sinI))^2))
                    end
                end



                dvPar_w += Ftot * (1.0 + alpha*sum_g)
                dvParSq_w += w*Ftot * (1.0 + alpha*sum_g)
                dvPerSq_w += w*Ftot * (1.0 + alpha*sum_g)
                sinSqdvPerSq_w += w*Ftot * (0.5 + alpha*sum_sinSqg)
                cosdvPerSq_w += w*Ftot * (0.5 + alpha*sum_cosg)

            end

            dvPar_w *= wmax
            dvParSq_w *= wmax
            dvPerSq_w *= wmax
            sinSqdvPerSq_w *= wmax
            cosdvPerSq_w *= wmax



            dvPar_phi += dvPar_w
            dvParSq_phi += dvParSq_w
            dvPerSq_phi += dvPerSq_w
            sinSqdvPerSq_phi += sinSqdvPerSq_w
            cosdvPerSq_phi += cosdvPerSq_w


        end



        dvPar += 2.0*sinvarphi*cosvarphi*dvPar_phi
        dvParSq += sinvarphi^3*dvParSq_phi
        dvPerSq += sinvarphi*(1.0+cosvarphi^2)*dvPerSq_phi
        sinSqdvPerSq += sinvarphi*(1.0+cosvarphi^2)*sinSqdvPerSq_phi
        cosdvPerSq += sinvarphi*(1.0+cosvarphi^2)*cosdvPerSq_phi


    end

    # dvPar *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvParSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)
    # dvPerSq *= pi/nbK*2.0*pi/nbK*(1.0/nbK)

    pref = 2.0*pi^2/nbK^3
    pref *= 2.0*pi*_G^2*logCoulomb

    dvPar *= -pref * (m_field + m_test)
    dvParSq *= 2.0*pref * m_field
    dvPerSq *= 2.0*pref * m_field
    sinSqdvPerSq *= 2.0*pref * m_field
    cosdvPerSq *= 2.0*pref * m_field

    return dvPar, dvParSq, dvPerSq, sinSqdvPerSq, cosdvPerSq
end



function localVelChange3DAngleAverage_compare(r::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, m_field::Float64, alpha::Float64=alphaRot,
                        nbAvrTh::Int64=nbAvrTh_default, nbK::Int64=nbK_default,
                        m_test::Float64=m_field)

    dvPar_n, dvParSq_n, dvPerSq_n, sinSqdvPerSq_n = localVelChange3DAngleAverage(r,vr,vt,cosI,m_field,alpha,nbAvrTh,nbK,m_test)
    dvPar, dvParSq, dvPerSq, sinSqdvPerSq = localVelChange3DAngleAverageExact(r,vr,vt,cosI,m_field,alpha,nbAvrTh,nbK,m_test)

    println("Relative error :")
    println("dvPar        : ",(dvPar_n-dvPar)/dvPar)
    println("dvParSq      : ",(dvParSq_n-dvParSq)/dvParSq)
    println("dvParSq      : ",(dvPerSq_n-dvPerSq)/dvPerSq)
    println("sinSqdvParSq : ",(sinSqdvPerSq_n-sinSqdvPerSq)/sinSqdvPerSq)

end