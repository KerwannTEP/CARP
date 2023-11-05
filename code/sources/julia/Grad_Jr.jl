# Compute the gradient of Jr(E,L) by sampling linearly in u-anomaly
# This should be done for orbits which are not too radial
function grad_Jr_E_L_Linear(rp::Float64, ra::Float64, nbu::Int64 = 300)
    E, L = E_L_from_rp_ra(rp,ra)
    sp, sa = sp_sa_from_rp_ra(rp,ra)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)


    dJrdE = 0.0
    dJrdL = 0.0
    d2JrdE2 = 0.0
    d2JrdEL = 0.0
    d2JrdL2 = 0.0

    for iu=1:nbu
        u = -1.0 + (2/nbu) * (iu-0.5)
        su = s_from_u_sma_ecc(u,sma,ecc)
        ru = r_from_s(su)
        xu = ru/_b
        jac = Theta(u,sp,sa)

        djacdE, djacdL, dsdL = djac_and_ds(u,sp,sa)

        dJrdE += jac
        dJrdL += jac/ru^2
        d2JrdE2 += djacdE
        d2JrdEL += djacdL
        d2JrdL2 += jac/ru^2 + L/ru^4*(djacdL*ru^2 - 2.0*_b^2*jac*su*dsdL)

    end

    dJrdE *= (2/nbu)*1/pi
    dJrdL *= (2/nbu)*(-L/pi)
    d2JrdE2 *= (2/nbu)*1/pi
    d2JrdEL *= (2/nbu)*1/pi
    d2JrdL2 *= (2/nbu)*(-1/pi)

    return dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2
end

# Compute the gradient of Jr(E,L) by sampling semi-logarithmically
# Integral is done in two part: 
# - near pericentre (log u-sampling)
# - away from pericentre (linear u-sampling)
# This should be done for orbits which are nearly radial
function grad_Jr_L_logIntegral(rp::Float64, ra::Float64, nbv::Int64 = nbu0, eps::Float64=10^(-5),
            Lcutoff::Float64=0.0005)

    E, L = E_L_from_rp_ra(rp,ra)

    if (L >= Lcutoff)
        sp, sa = sp_sa_from_rp_ra(rp,ra)
        sma, ecc = sma_ecc_from_sp_sa(sp,sa)
        uminPlusOne = eps*(sp^2-1.0)^2/(3.0*(sa-1.0))

        vmin = log(eps) + 4.0*log(rp) - (log(3.0) + log(sa-1.0))



        dJrdL = 0.0
        d2JrdL2 = 0.0

        for iv=1:nbv
            v = vmin + (log(2.0)-vmin)/nbv * (iv-0.5)
            u = exp(v)-1.0
            su = s_from_u_sma_ecc(u,sma,ecc)
            ru = r_from_s(su)
            xu = ru/_b
            jac = Theta(u,sp,sa)

            djacdE, djacdL, dsdL = djac_and_ds(u,sp,sa)
            expv = exp(v)

            dJrdL += jac*expv/ru^2

            d2JrdL2 += jac*expv/ru^2 + L*expv/ru^4*(djacdL*ru^2 - 2.0*_b^2*jac*su*dsdL)


        end

        dJrdL *= (log(2.0)-vmin)/nbv * (-L/pi)
        d2JrdL2 *= (log(2.0)-vmin)/nbv * (-1/pi)

        # Leftover integral
        u = uminPlusOne/2.0
        su = s_from_u_sma_ecc(u,sma,ecc)
        ru = r_from_s(su)
        jac = Theta(u,sp,sa)
        djacdE, djacdL, dsdL = djac_and_ds(u,sp,sa)

        dJrdL += (uminPlusOne)*jac/ru^2
        d2JrdL2 += (uminPlusOne)*(jac/ru^2 + L/ru^4*(djacdL*ru^2 - 2.0*_b^2*jac*su*dsdL))


        return dJrdL, d2JrdL2

    else # near-radial orbits
        # maybe we should see what is the exact limit for L=0
        # compute at cutoff for now


        spcut, sacut = sp_sa_from_E_L(E,Lcutoff)
        rpcut, racut = rp_ra_from_sp_sa(spcut,sacut)
        smacut, ecccut = sma_ecc_from_sp_sa(spcut,sacut)

        uminPlusOnecut = eps*(spcut^2-1.0)^2/(3.0*(sacut-1.0))

        rpcut, racut = rp_ra_from_sp_sa(spcut,sacut)

        vmincut = log(eps) + 4.0*log(Lcutoff)-2.0*log(E+1.0)-2.0*log(2.0) - log(3.0) - log(sacut-1.0)

        dJrdLcut = 0.0

        for iv=1:nbv
            v = vmincut + (log(2.0)-vmincut)/nbv * (iv-0.5)
            u = exp(v)-1.0
            su = s_from_u_sma_ecc(u,smacut,ecccut)
            ru = r_from_s(su)
            xu = ru/_b
            jac = Theta(u,spcut,sacut)


            dJrdLcut += jac*exp(v)/ru^2

        end

        dJrdLcut *= (log(2.0)-vmincut)/nbv * (-Lcutoff/pi)

        # Leftover integral
        u = uminPlusOnecut/2.0 - 1.0
        su = s_from_u_sma_ecc(u,smacut,ecccut)
        ru = r_from_s(su)
        jac = Theta(u,spcut,sacut)

        dJrdLcut += (uminPlusOnecut)*jac/ru^2

        # Interpolate linearly dJrdL

        dJrdL = -0.5 + L/Lcutoff * (dJrdLcut+0.5)
        d2JrdL2 = (dJrdLcut+0.5)/Lcutoff

        return dJrdL, d2JrdL2


    end
end

# Wrapping the computation of the gradient of Jr(E,L) for all orbits
function grad_Jr_E_L(E::Float64, L::Float64, nbu::Int64 = nbu0, eps::Float64=10^(-5), Lcutoff::Float64=0.0005)
    sp, sa = sp_sa_from_E_L(E,L)
    rp, ra = rp_ra_from_sp_sa(sp,sa)

    saCutoff = sp + 0.005
    _, raCutoff = rp_ra_from_sp_sa(sp,saCutoff)
    if (sa < saCutoff) # near circular

        # compute exact circular
        dJrdE_circ = 1/_Omega0 * sp^(5/2)/sqrt(3.0+sp^2)
        dJrdL_circ = -sp/sqrt(3.0+sp^2)
        d2JrdE2_circ = -(3.0*sp^(7/2)*(-5.0+35.0*sp^2+16.0*sp^4+2.0*sp^6))/(2.0*_Omega0*_E0*(3.0+sp^2)^(7/2))
        d2JrdEL_circ = (15.0*sp^2*(sp^2-1.0))/(2.0*_E0*(3.0+sp^2)^(7/2))
        d2JrdL2_circ = -3.0/(2.0*_L0) * sp^(1/2)*(5.0+7.0*sp^2+4.0*sp^4)/(3.0+sp^2)^(7/2)

        # compute value at cutoff sa

        if (L >= 0.05)
            dJrdECut, dJrdLCut, d2JrdE2Cut, d2JrdELCut, d2JrdL2Cut = grad_Jr_E_L_Linear(rp,raCutoff,nbu)
        else
            dJrdECut, _, d2JrdE2Cut, d2JrdELCut, _ = grad_Jr_E_L_Linear(rp,raCutoff,nbu)
            dJrdLCut, d2JrdL2Cut = grad_Jr_L_logIntegral(rp,raCutoff,nbu,eps,Lcutoff)
        end

        # taylor expansion near circular orbit (sp,sp)
        # smooth junction at cutoff sa=sp+0.0001

        dJrdE = dJrdE_circ+(dJrdECut-dJrdE_circ)*(sa-sp)/(saCutoff-sp)
        dJrdL = dJrdL_circ+(dJrdLCut-dJrdL_circ)*(sa-sp)/(saCutoff-sp)
        d2JrdE2 = d2JrdE2_circ+(d2JrdE2Cut-d2JrdE2_circ)*(sa-sp)/(saCutoff-sp)
        d2JrdEL = d2JrdEL_circ+(d2JrdELCut-d2JrdEL_circ)*(sa-sp)/(saCutoff-sp)
        d2JrdL2 = d2JrdL2_circ+(d2JrdL2Cut-d2JrdL2_circ)*(sa-sp)/(saCutoff-sp)

        return dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2

    else


        if (L >= 0.05)

            return grad_Jr_E_L_Linear(rp,ra,nbu)
        else

            dJrdE, _, d2JrdE2, d2JrdEL, _ = grad_Jr_E_L_Linear(rp,ra,nbu)
            dJrdL, d2JrdL2 = grad_Jr_L_logIntegral(rp,ra,nbu,eps,Lcutoff)
            return dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2
        end
    end
end
