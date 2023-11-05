
function alpha_beta_from_E_L(E::Float64, L::Float64, nbu::Int64 = nbu0)
    sp, sa = sp_sa_from_E_L(E,L)
    rp, ra = rp_ra_from_sp_sa(sp,sa)
    return alpha_beta_from_rp_ra(rp,ra,nbu)
end

function alpha_beta_from_rp_ra(rp::Float64, ra::Float64, nbu::Int64 = nbu0)
    E, L = E_L_from_rp_ra(rp,ra)
    sp, sa = sp_sa_from_rp_ra(rp,ra)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)

    inv_alpha = 0.0
    beta = 0.0
    for iu=1:nbu
        u = -1.0 + (2/nbu) * (iu-0.5)
        su = s_from_u_sma_ecc(u,sma,ecc)
        ru = r_from_s(su)
        jac = Theta(u,sp,sa)
        inv_alpha += jac
        if (rp != 0.0) # not radial
            beta += jac/ru^2
        end
    end
    inv_alpha *= (2/nbu) * (_Omegam/pi)
    if (L != 0.0) # not radial
        beta *= (2/nbu) * (L/pi)
    else
        beta = 0.5 # approximate by the limit. should be find a better taylor expansion?
    end

    return 1/inv_alpha, beta

end

function beta_from_E_L_logIntegral(E::Float64, L::Float64, nbv::Int64 = nbu0, eps::Float64=10^(-5),
            Lcutoff::Float64=0.00005)
    sp, sa = sp_sa_from_E_L(E,L)
    rp, ra = rp_ra_from_sp_sa(sp,sa)
    return beta_from_rp_ra_logIntegral(rp,ra,nbv,eps,Lcutoff)
end

function beta_from_rp_ra_logIntegral(rp::Float64, ra::Float64, nbv::Int64 = nbu0, eps::Float64=10^(-5),
            Lcutoff::Float64=0.00005)
    E, L = E_L_from_rp_ra(rp,ra)

    if (L >= Lcutoff)
        sp, sa = sp_sa_from_rp_ra(rp,ra)
        sma, ecc = sma_ecc_from_sp_sa(sp,sa)
        uminPlusOne = eps*(sp^2-1.0)^2/(3.0*(sa-1.0))

        vmin = log(eps) + 4.0*log(rp) - (log(3.0) + log(sa-1.0))


        beta = 0.0


        for iv=1:nbv
            v = vmin + (log(2.0)-vmin)/nbv * (iv-0.5)
            u = exp(v)-1.0
            su = s_from_u_sma_ecc(u,sma,ecc)
            ru = r_from_s(su)
            jac = Theta(u,sp,sa)
            if (rp != 0.0) # not radial
                beta += jac*exp(v)/ru^2
            end
        end
        beta *= (log(2.0)-vmin)/nbv * (L/pi)

        # Leftover integral
        u = uminPlusOne/2.0
        su = s_from_u_sma_ecc(u,sma,ecc)
        ru = r_from_s(su)
        beta += (uminPlusOne)*L*Theta(u,sp,sa)/(pi*ru^2)

        return beta

    elseif (L != 0.0)

        # compute beta(Lcutoff)

        spcut, sacut = sp_sa_from_E_L(E,Lcutoff)
        rpcut, racut = rp_ra_from_sp_sa(spcut,sacut)
        smacut, ecccut = sma_ecc_from_sp_sa(spcut,sacut)

        uminPlusOne = eps*(spcut^2-1.0)^2/(3.0*(sacut-1.0))

        vmin = log(eps) + 4.0*log(rpcut) - (log(3.0) + log(sacut-1.0))


        betacut = 0.0
        for iv=1:nbv
            v = vmin + (log(2.0)-vmin)/nbv * (iv-0.5)
            u = exp(v)-1.0
            su = s_from_u_sma_ecc(u,smacut,ecccut)
            ru = r_from_s(su)
            jac = Theta(u,spcut,sacut)
            if (rp != 0.0) # not radial
                betacut += jac*exp(v)/ru^2
            end
        end
        betacut *= (log(2.0)-vmin)/nbv * (Lcutoff/pi)

        # Leftover integral
        u = uminPlusOne/2.0
        su = s_from_u_sma_ecc(u,smacut,ecccut)
        ru = r_from_s(su)
        betacut += (uminPlusOne)*Lcutoff*Theta(u,spcut,sacut)/(pi*ru^2)

        # use linear taylor expansion near L=0
        # (beta-0.5)/(betacut-0.5) = L/Lcut, hence
        # beta = 0.5 + L/Lcut * (betacut-0.5)

        beta = 0.5 + (betacut-0.5)*L/Lcutoff

        return beta

    else # radial orbits
        beta = 0.5 # approximate by the limit. should be find a better taylor expansion?

        return beta
    end


end
