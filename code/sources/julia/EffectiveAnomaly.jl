##########################################
# Orbit-average functions
##########################################

function Theta(u::Float64, sp::Float64, sa::Float64)



    A = sp * (u+2.0)*(u-1.0)^2 - sa*(u-2.0)*(u+1.0)^2
    B = sp * (u^3-3.0*u+6.0) - sa*(u^3-3.0*u-6.0)


    return (1.0/_Omega0 * 3.0/(4.0*sqrt(2.0)) * sqrt(sa*sp*(sa+sp))/sqrt(4.0-u^2)
            * A^(1.5)/sqrt(sa*sp*A + B))
end

function s_from_u_sma_ecc(u::Float64, sma::Float64, ecc::Float64)
    fu = u * (1.5 - 0.5*u^2)
    return sma * (1.0 + ecc * fu)
end

function r_from_s(s::Float64)
    return _b * sqrt(abs(s^2 - 1.0))
end

##########################################
# Gradients
##########################################


function _djacdsp(u::Float64, sp::Float64, sa::Float64)

    u2 = u*u
    u3 = u*u*u
    u4 = u2*u2
    u5 = u2*u3
    u6 = u3*u3
    uMinus1_2 = (-1 + u)*(-1 + u)
    uPlus1_2 = (1 + u)*(1 + u)
    uMinus1_4 = uMinus1_2*uMinus1_2


    num = (3 *(4 *sp^3 *uMinus1_2 *(12 - 3 *u2 + 2 *u3 + u4) +
           sa *sp^2 *(108 - 120 *u - 27 *u2 + 40 *u3 + 18 *u4 - 3 *u6 +
           3 *sp^2 *uMinus1_4 *(2 + u)^2) -
           2 *sa^2 *sp *(-6 - 3 *u + u3) *(6 - 3 *u + u3 +
           sp^2 *(2 - 3 *u + u3)) -
           sa^3 *(-2 + u) *uPlus1_2 *(6 + 3* u - u3 + sp^2 *(6 - 3 *u + u3))))

    den = (8 *sqrt(2)* sp *(sa +
           sp) *sqrt(-(((-4 + u2) *(sa^2 *sp *(-2 + u) *uPlus1_2 -
           sp *(6 - 3 *u + u3) +
           sa *(-6 - 3 *u + u3 - sp^2 *(2 - 3* u + u3))))/(
           sa *sp *(sa + sp) *(sa *(-2 + u) *uPlus1_2 -
           sp *(2 - 3 *u + u3)))))* (sa^2* sp *(2 + 3 *u - u3) +
           sp *(6 - 3* u + u3) + sa *(6 + 3*u - u3 + sp^2 *(2 - 3 *u + u3))))

    return 1.0/_Omega0 * num/den
end

function _djacdsa(u::Float64, sp::Float64, sa::Float64)

    u2 = u*u
    u3 = u*u*u
    u4 = u2*u2
    u5 = u2*u3
    u6 = u3*u3
    sa4 = sa*sa*sa*sa
    uPlus1_2 = (1 + u)*(1 + u)
    uPlus1_4 = uPlus1_2*uPlus1_2

    num = (3 *(-4 + u2) *(3 *sa4 *sp *(-2 + u)^2 *uPlus1_4 +
           sp^3 *(-1 + u)^2 *(12 - 3 *u2 + 2 *u3 + u4) -
           2 *sa *sp^2 *(-36 + 9 *u2 - 6 *u4 + u6) -
           2 *sa^3 *(-2 + u) *(1 + u)^2 *(12 + 6 *u - 2 *u3 +
           sp^2 *(6 - 3 *u + u3)) +
           sa^2 *sp *(108 + 120 *u - 27 *u2 - 40*u3 + 18 *u4 - 3 *u6 -
           sp^2 *(-12 + 12 *u + 9 *u2 - 4 *u3 - 6 *u4 + u6))))

    den = (8 *sqrt(2)*
           sa^2 *sp *(sa + sp)^2 *(sa *(-2 + u) *(1 + u)^2 -
           sp *(2 - 3 *u + u3)) *(-(((-4 + u2) *(sa^2 *sp *(-2 + u) *(1 + u)^2 -
           sp *(6 - 3* u + u3) +
           sa *(-6 - 3* u + u3 - sp^2 *(2 - 3 *u + u3))))/(
           sa *sp* (sa + sp) *(sa *(-2 + u) *(1 + u)^2 - sp *(2 - 3 *u + u3)))))^(3/2))

    return 1.0/_Omega0 * num/den
end


function djac_and_ds(u::Float64, sp::Float64, sa::Float64)
    rp, ra = rp_ra_from_sp_sa(sp,sa)
    E, L = E_L_from_rp_ra(rp,ra)

    if (sp == 1.0) # Take the limit rp -> 0+
        dPHIdsa, dPHIdLa = dpsieff_s(sa,L)
        dspdE = 0.0
        dsadE = 1.0/dPHIdsa
        dspdL = 0.0
        dsadL = -dPHIdLa/dPHIdsa
    else
        dPHIdsp, dPHIdLp = dpsieff_s(sp,L)
        dPHIdsa, dPHIdLa = dpsieff_s(sa,L)
        dspdE = 1.0/dPHIdsp
        dsadE = 1.0/dPHIdsa
        dspdL = -dPHIdLp/dPHIdsp
        dsadL = -dPHIdLa/dPHIdsa
    end

    fu = u * (1.5 - 0.5*u^2)

    djacdsp = _djacdsp(u,sp,sa)
    djacdsa = _djacdsa(u,sp,sa)

    dsdsp = (1.0 - fu)/2.0
    dsdsa = (1.0 + fu)/2.0

    djacdE = dsadE*djacdsa + dspdE*djacdsp
    djacdL = dsadL*djacdsa + dspdL*djacdsp

    dsdE = dsadE*dsdsa+dspdE*dsdsp
    dsdL = dsadL*dsdsa+dspdL*dsdsp

    return djacdE, djacdL, dsdL, dsdE
end
