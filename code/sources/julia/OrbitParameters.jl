
##########################################
# Computes the change of variables formulae
##########################################

# CAUTION: NOT ALL sma, ecc are available
# sma must be greater than 1.0
# ecc must be between 0.0 and 1.0-1/a

function sma_ecc_from_sp_sa(sp::Float64, sa::Float64)
    sma = (sa + sp)/2
    ecc = (sa - sp)/(sa + sp)
    return sma, ecc
end

function sp_sa_from_sma_ecc(sma::Float64, ecc::Float64)
    sp = sma * (1.0 - ecc) # must be > 1.0 <=> only ecc < 1.0-1/a allowed
    sa = sma * (1.0 + ecc)
    return sp, sa
end

##########################################
##########################################

function smar_eccr_from_rp_ra(rp::Float64, ra::Float64)
    sma = (ra + rp)/2
    ecc = (ra - rp)/(ra + rp)
    return sma, ecc
end

function rp_ra_from_smar_eccr(smar::Float64, eccr::Float64)
    rp = smar * (1.0 - eccr) # must be > 1.0 <=> only ecc < 1.0-1/a allowed
    ra = smar * (1.0 + eccr)
    return rp, ra
end

##########################################
##########################################

function sp_sa_from_rp_ra(rp::Float64, ra::Float64)
    sp = sqrt(1.0 + (rp/_b)^2)
    sa = sqrt(1.0 + (ra/_b)^2)
    return  sp, sa
end

function rp_ra_from_sp_sa(sp::Float64, sa::Float64)
    rp = _b * sqrt(abs(sp^2 - 1.0))
    ra = _b * sqrt(abs(sa^2 - 1.0))
    return  rp, ra
end

##########################################
##########################################

function E_L_from_sp_sa(sp::Float64, sa::Float64)

    E = _E0/sp - _E0*(sa^2-1.0)/(sa*sp*(sa+sp))
    L = _L0*sqrt(2.0*(sp^2-1.0)*(sa^2-1.0)/(sa*sp*(sa+sp)))

    return E, L
end

function E_L_from_rp_ra(rp::Float64, ra::Float64)

    sp, sa = sp_sa_from_rp_ra(rp,ra)

    return E_L_from_sp_sa(sp,sa)

end

##########################################
##########################################

# Use bissection
# We know that 1 < sp < sc and sc < sa < E0/E for bound non-circular non-radial orbits
# upper bound should be E0/E + 1 in order to get a proper bracket
# should find some proper cutoff for quasi circular orbits
function sp_sa_from_E_L(E::Float64, L::Float64)
    tE = E/_E0
    tL = L/_L0

    if (tE >= 0.0) # bound orbits
        sc = _sc(tE)
        if (L >= _Lc(E)) # circular orbit (inequality to take care of small Float errors)
            return sc, sc
        elseif (L == 0.0) # radial orbit
            return 1.0, 1.0/tE
        else # arbitrary orbit
            fct = (s -> tE*s^3 - s^2 + (0.5*tL^2-tE)*s + 1.0)
            if (abs(fct(sc)) <= 10^(-10)) # quasi circular
                return sc, sc
            else # not quasi-circular

                # find sp
                if (abs(fct(1.0)) <= 10.0*eps(Float64))
                    sp = 1.0
                else
                    sp = bisection(fct,1.0,sc)
                end

                #find sa
                if (abs(fct(1.0/tE)) <= 10.0*eps(Float64))
                    sa = 1.0/tE
                else
                    sa = bisection(fct,sc,1.0/tE+1)
                end

                return sp, sa
            end
        end
    else # unbound orbits
        if (L == 0.0) # radial orbits
            return 1.0, Inf
        else # arbitrary orbit
            fct = (s -> tE*s^3 - s^2 + (0.5*tL^2-tE)*s + 1.0)
            sp = bisection(fct,1.0,1.0/(0.5*tL^2-tE))
            return sp, Inf

        end
    end


end

##########################################
##########################################

function Jr_from_rp_ra(rp::Float64, ra::Float64, nbu::Int64 = nbu0)
    E, L = E_L_from_rp_ra(rp,ra)
    sp, sa = sp_sa_from_rp_ra(rp,ra)
    sma, ecc = sma_ecc_from_sp_sa(sp,sa)

    sum = 0.0
    for iu=1:nbu
        u = -1.0 + (2/nbu) * (iu-0.5)
        su = s_from_u_sma_ecc(u,sma,ecc)
        ru = r_from_s(su)
        vrSq = _vrSq(ru,E,L)
        jac = Theta(u,sp,sa)
        sum += vrSq*jac
    end
    sum *= (2/nbu) * (1.0/pi)

    return sum
end

##########################################
##########################################

function Jr_from_E_L(E::Float64, L::Float64, nbu::Int64 = nbu0)

    if (L > _Lc(E))
        return Inf
    else
        sp, sa = sp_sa_from_E_L(E,L)
        sma, ecc = sma_ecc_from_sp_sa(sp,sa)
        sum = 0.0
        for iu=1:nbu
            u = -1.0 + (2/nbu) * (iu-0.5)
            su = s_from_u_sma_ecc(u,sma,ecc)
            ru = r_from_s(su)
            vrSq = vrSq_from_E_L(ru,E,L)
            jac = Theta(u,sp,sa)
            sum += vrSq*jac
        end
        sum *= (2.0/nbu) * (1.0/pi)
        return sum
    end
end

function E_from_Jr_L(Jr::Float64, L::Float64, nbu::Int64=nbu0, eps::Float64=4.0*10^(-15))
    En = -0.5 # E0
    delta = 0.5 # initial precision
    JrE = 0.0
    if (Jr==0.0 && L==0.0)
        return -1.0
    else
        while (delta>eps) # while E is not a solution precise enough
            delta /= 2
            if (L > _Lc(En)) # if (En,L) is not an orbit, we move towards the accepted (E,L) region
                En += delta
            else # if within the accepted (E,L) region, use the fact that Jr(E,L) increases with E
                JrE = Jr_from_E_L(En,L,nbu)
                if (JrE < Jr)
                    En += delta
                else
                    En -= delta
                end
            end
        end
        return En
    end
end

##########################################
##########################################

function alpha_beta_from_E_L(E::Float64, L::Float64, nbu::Int64 = 300, eps::Float64=10^(-5), Lcutoff::Float64=0.00005)
    sp, sa = sp_sa_from_E_L(E,L)
    rp, ra = rp_ra_from_sp_sa(sp,sa)
    if (L >= 0.05)
        return alpha_beta_from_rp_ra_Linear(rp,ra,nbu)

    else
        alpha, _ = alpha_beta_from_rp_ra_Linear(rp,ra,nbu)
        beta = beta_from_rp_ra_logIntegral(rp,ra,nbu,eps,Lcutoff)
        return alpha, beta
    end
end

function E_L_from_alpha_beta(alpha::Float64, beta::Float64, nbu::Int64 = nbu0,
            eps::Float64=4.0*10^(-10), nbStepMax::Int64=10)

    alphac = alpha_circ(beta)

    if (alpha == alphac)
        return E_L_from_beta_Circ(beta)
    else
        if (beta == 0.5)
            return E_L_from_alpha_beta_Radial(alpha)
        else
            return E_L_from_alpha_beta_Arbitrary(alpha,beta,alphac,nbu,eps,nbStepMax)
        end
    end
end


##########################################
# Compactified radius
##########################################

function xi_from_r(r::Float64)
    x = r/_b
    return (x^2 - 1.0)/(x^2 + 1.0)
end

function r_from_xi(xi::Float64)
    return _b * sqrt((1.0 + xi)/(1.0 - xi))
end
