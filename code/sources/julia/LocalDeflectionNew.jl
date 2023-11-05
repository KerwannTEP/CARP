function _wmaxNew(theta::Float64, varphi::Float64, phi::Float64, E::Float64,
            vr::Float64, vt::Float64)

    b = sin(varphi) * (vr*cos(phi-theta) + vt*sin(phi-theta))
    discr = b^2 - 2.0*E

    return -b + sqrt(discr)
end

function _h_New(r::Float64, theta::Float64, vr::Float64, vt::Float64,
                    cosI::Float64, sinI::Float64, alpha::Float64=alphaRot,
                    nbK::Int64=nbK_default)

    vSq = vr^2 + vt^2
    E = psi(r) + vSq/2.0
    L = r*vt
    Lz = L*cosI

    h = 0.0

    for ivarphi=1:nbK
        varphi = pi*(ivarphi-0.5)/nbK

        h_phi = 0.0

        for iphi=1:nbK
            phi = 2.0*pi*(iphi-0.5)/nbK
            wmax = _wmaxNew(theta,varphi,phi,E,vr,vt)

            h_w = 0.0

            for iw=1:nbK
                w = wmax*(iw-0.5)/nbK

                v_dot_w = w*sin(varphi)*(vr*cos(phi-theta)+vt*sin(phi-theta))

                Ep  = E + w^2/2.0 - v_dot_w
                Lp  = r*sqrt(vSq + w^2 - 2.0*v_dot_w - (vr-w*sin(varphi)*cos(phi-theta))^2)
                Lzp = Lz - r*w*sin(varphi)*sin(phi-theta)*cosI - r*w*sin(theta)*cos(varphi)*sinI

                Frotp = _Frot(Ep,Lp,Lzp,alpha)

                h_w += w*Frotp

            end
            h_w *= wmax/nbK

            h_phi += h_w

        end
        h_phi *= 2.0*pi/nbK

        h += sin(varphi)*h_phi

    end

    h *= pi/nbK

    return h

end

function RosenbluthPotentialsNew(r::Float64, theta::Float64, vr::Float64, vt::Float64,
                    cosI::Float64, sinI::Float64, alpha::Float64=alphaRot,
                    nbK::Int64=nbK_default)

    vSq = vr^2 + vt^2
    E = psi(r) + vSq/2.0
    L = r*vt
    Lz = L*cosI

    dhdvr = 0.0
    dhdvt = 0.0
    dgdvt = 0.0
    d2gdvr2 = 0.0
    d2gdvrvt = 0.0
    d2gdvt2 = 0.0

    for ivarphi=1:nbK
        varphi = pi*(ivarphi-0.5)/nbK

        dhdvr_phi = 0.0
        dhdvt_phi = 0.0
        dgdvt_phi = 0.0
        d2gdvr2_phi = 0.0
        d2gdvrvt_phi = 0.0
        d2gdvt2_phi = 0.0

        for iphi=1:nbK
            phi = 2.0*pi*(iphi-0.5)/nbK
            wmax = _wmaxNew(theta,varphi,phi,E,vr,vt)

            dhdvr_w = 0.0
            dhdvt_w = 0.0
            dgdvt_w = 0.0
            d2gdvr2_w = 0.0
            d2gdvrvt_w = 0.0
            d2gdvt2_w = 0.0

            for iw=1:nbK
                w = wmax*(iw-0.5)/nbK

                v_dot_w = w*sin(varphi)*(vr*cos(phi-theta)+vt*sin(phi-theta))

                Ep  = E + w^2/2.0 - v_dot_w
                Lp  = r*sqrt(vSq + w^2 - 2.0*v_dot_w - (vr-w*sin(varphi)*cos(phi-theta))^2)
                Lzp = Lz - r*w*sin(varphi)*sin(phi-theta)*cosI - r*w*sin(theta)*cos(varphi)*sinI

                Frotp = _Frot(Ep,Lp,Lzp,alpha)

                dhdvr_w += Frotp
                dhdvt_w += Frotp
                dgdvt_w += w^2*Frotp
                d2gdvr2_w += w*Frotp
                d2gdvrvt_w += w*Frotp
                d2gdvt2_w += w*Frotp

            end
            dhdvr_w *= wmax/nbK
            dhdvt_w *= wmax/nbK
            dgdvt_w *= wmax/nbK
            d2gdvr2_w *= wmax/nbK
            d2gdvrvt_w *= wmax/nbK
            d2gdvt2_w *= wmax/nbK

            dhdvr_phi += cos(phi-theta)*dhdvr_w
            dhdvt_phi += sin(phi-theta)*dhdvt_w
            dgdvt_phi += sin(phi-theta)*dgdvt_w
            d2gdvr2_phi += (1.0 - sin(varphi)^2*cos(phi-theta)^2)*d2gdvr2_w
            d2gdvrvt_phi += (1.0 - 0.5*sin(varphi)^2*sin(2.0*(phi-theta)))*d2gdvrvt_w
            d2gdvt2_phi += (1.0 - sin(varphi)^2*sin(phi-theta)^2)*d2gdvt2_w

        end
        dhdvr_phi *= 2.0*pi/nbK
        dhdvt_phi *= 2.0*pi/nbK
        dgdvt_phi *= 2.0*pi/nbK
        d2gdvr2_phi *= 2.0*pi/nbK
        d2gdvrvt_phi *= 2.0*pi/nbK
        d2gdvt2_phi *= 2.0*pi/nbK

        dhdvr += sin(varphi)^2*dhdvr_phi
        dhdvt += sin(varphi)^2*dhdvt_phi
        dgdvt += sin(varphi)^2*dgdvt_phi
        d2gdvr2 += sin(varphi)*d2gdvr2_phi
        d2gdvrvt += sin(varphi)*d2gdvrvt_phi
        d2gdvt2 += sin(varphi)*d2gdvt2_phi

    end

    dhdvr *= -pi/nbK
    dhdvt *= -pi/nbK
    dgdvt *= pi/nbK
    d2gdvr2 *= pi/nbK
    d2gdvrvt *= pi/nbK
    d2gdvt2 *= pi/nbK

    return dhdvr, dhdvt, dgdvt, d2gdvr2, d2gdvrvt, d2gdvt2

end


function RosenbluthPotentialsNewXYZ(r::Float64, theta::Float64, vr::Float64, vt::Float64,
                    cosI::Float64, sinI::Float64, alpha::Float64=alphaRot,
                    nbK::Int64=nbK_default)

    vSq = vr^2 + vt^2
    E = psi(r) + vSq/2.0
    L = r*vt
    Lz = L*cosI

    dhdv1 = 0.0


    for ivarphi=1:nbK
        varphi = pi*(ivarphi-0.5)/nbK

        dhdv1_phi = 0.0


        for iphi=1:nbK
            phi = 2.0*pi*(iphi-0.5)/nbK
            wmax = _wmaxNew(theta,varphi,phi,E,vr,vt)

            dhdv1_w = 0.0

            for iw=1:nbK
                w = wmax*(iw-0.5)/nbK

                v_dot_w = w*sin(varphi)*(vr*cos(phi-theta)+vt*sin(phi-theta))

                Ep  = E + w^2/2.0 - v_dot_w
                Lp  = r*sqrt(vSq + w^2 - 2.0*v_dot_w - (vr-w*sin(varphi)*cos(phi-theta))^2)
                Lzp = Lz - r*w*sin(varphi)*sin(phi-theta)*cosI - r*w*sin(theta)*cos(varphi)*sinI

                Frotp = _Frot(Ep,Lp,Lzp,alpha)

                dhdv1_w += Frotp


            end
            dhdv1_w *= wmax/nbK


            dhdv1_phi += cos(phi-theta)*dhdv1_w


        end
        dhdv1_phi *= 2.0*pi/nbK


        dhdv1 += sin(varphi)^2*dhdv1_phi


    end

    dhdv1 *= -pi/nbK

    return dhdv1

end





function localVelChangeNew(r::Float64, theta::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, sinI::Float64, m_field::Float64,
                        alpha::Float64=alphaRot, nbK::Int64=nbK_default,
                        m_test::Float64=m_field)

    dhdvr, dhdvt, dgdvt, d2gdvr2, d2gdvrdvt, d2gdvt2 = RosenbluthPotentialsNew(r,theta,vr,vt,cosI,sinI,alpha,nbK)

    cst       = 4.0*PI*_G^2*logCoulomb
    v         = sqrt(vr^2+vt^2)

    # if (vt>0)
    vr_v   = vr/v
    vt_v   = vt/v

    dvPar   = cst*(m_field+m_test) *(vr_v*dhdvr+vt_v*dhdvt)
    dvPar2  = cst* m_field         *(vr_v^2*d2gdvr2+(2*vr_v*vt_v)*d2gdvrdvt
                                    +vt_v^2*d2gdvt2)
    dvPerp2 = cst* m_field         *(vt_v^2*d2gdvr2-(2*vr_v*vt_v)*d2gdvrdvt
                                    +vr_v^2*d2gdvt2+(1/vt)*dgdvt)


    return dvPar, dvPar2, dvPerp2
end


function localVelChangeNewXYZ(r::Float64, theta::Float64, vr::Float64, vt::Float64,
                        cosI::Float64, sinI::Float64, m_field::Float64,
                        alpha::Float64=alphaRot, nbK::Int64=nbK_default,
                        m_test::Float64=m_field)

    dhdv1 = RosenbluthPotentialsNewXYZ(r,theta,vr,vt,cosI,sinI,alpha,nbK)

    cst       = 4.0*PI*_G^2*logCoulomb
    v         = sqrt(vr^2+vt^2)

    # if (vt>0)
    vr_v   = vr/v
    vt_v   = vt/v

    dvPar   = cst*(m_field+m_test) * vr_v * dhdv1



    return dvPar
end
