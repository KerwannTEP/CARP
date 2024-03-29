##########################################################
# Compute orbital parameters E and L from frequency parameters alpha and beta
##########################################################

# This tests if an orbit is near-circular
function test_near_circular(alpha::Float64, alphac::Float64, sc::Float64)
    d2JrdE2_circ = -(3.0*sc^(7/2)*(-5.0+35.0*sc^2+16.0*sc^4+2.0*sc^6))/(2.0*_Omega0*_E0*(3.0+sc^2)^(7/2))
    dE = (alpha-alphac)/(_Omegam*alphac^2*d2JrdE2_circ)
    return abs(dE/_E0) < 10^(-10)
end

# This computes the next step in E and L in the Newton inversion
function step_E_L(E::Float64, L::Float64, alpha::Float64, alpha_guess::Float64,
                beta::Float64, beta_guess::Float64)

    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = grad_Jr_E_L(E,L)

    dalphadE = -_Omegam*alpha_guess^2 * d2JrdE2
    dalphadL = -_Omegam*alpha_guess^2 * d2JrdEL
    dbetadE = -d2JrdEL
    dbetadL = -d2JrdL2

    detJac = dalphadE*dbetadL-dbetadE*dalphadL

    dAlpha = alpha-alpha_guess
    dBeta = beta-beta_guess

    stepE = (dbetadL*dAlpha-dalphadL*dBeta)/detJac
    stepL = (-dbetadE*dAlpha+dalphadE*dBeta)/detJac

    return stepE, stepL
end


# This computes the next values of E and L in any position in (E,L) space
function next_E_L(E::Float64, L::Float64, stepE::Float64, stepL::Float64)
    if (E+stepE >= 0.0)
        E = E/2.0
    elseif (E+stepE <= _E0)
        E = (E+_E0)/2.0
    else
        E = E + stepE
    end


    if (L+stepL > _Lc(E))
        L = (L+_Lc(E))/2.0
    elseif (L+stepL <= 0.0)
        L = L/2.0
    else
        L = L + stepL
    end
    return E, L
end

# This computes E and L from alpha for a radial orbit
function E_L_from_alpha_beta_Radial(alpha::Float64)
    # L=0.0
    # bissection on alpha_beta_from_E_L(E,0.0)
    # bisection(fun, xl, xu)
    L = 0.0
    Emin = -1.0 # alpha=1.0
    Emax = -1.0
    while (alpha_beta_from_E_L(Emax,L)[1] > alpha)
        Emax /= 2.0
    end

    E = bisection(E->alpha_beta_from_E_L(E,L)[1]-alpha, Emin, Emax)

    alpha_guess, beta_guess = alpha_beta_from_E_L(E,L)

    return E, L, alpha-alpha_guess, 0.0, 1
end

# This computes E and L from alpha for a circular orbit
function E_L_from_beta_Circ(beta::Float64)
    sc = beta*sqrt(3.0)/sqrt(1.0-beta^2)
    rp, ra = rp_ra_from_sp_sa(sc,sc)
    E, L = E_L_from_rp_ra(rp,ra)
    return E, L, 0.0, 0.0, 1
end

# This computes (E,L) from (alpha,beta)
function E_L_from_alpha_beta_Arbitrary(alpha::Float64, beta::Float64, alphac::Float64,
            nbu::Int64 = 300, eps::Float64=4.0*10^(-10), nbStepMax::Int64=10)

    # initial guess (E, L)
    sc = beta*sqrt(3.0)/sqrt(abs(1.0-beta^2))
    rc, _ = rp_ra_from_sp_sa(sc,sc)
    Ec, Lc = E_L_from_rp_ra(rc,rc)
    E, L = Ec, Lc

    alpha_guess, beta_guess = alphac, beta

    # Is the orbit near circular ?
    #dAlpha = dAlpha/dE dE
    # => dE = dAlpha / (dAlpha/dE) evaluated at circular orbits with Lc
    # dAlpha/dE = -Omega_m alphac^2 dJr2/dE2 : closed expression for circular orbits

    if (test_near_circular(alpha,alphac,sc))
        return E, L, alpha-alpha_guess, 0.0, 0
    else # Begin Newton's method
        # https://en.wikipedia.org/wiki/Newton%27s_method#k_variables,_k_functions

        nbStep = 0

        while (sqrt((alpha - alpha_guess)^2 + (beta_guess - beta)^2) > eps)

            stepE, stepL = step_E_L(E,L,alpha,alpha_guess,beta,beta_guess)
            E, L = next_E_L(E,L,stepE,stepL)

            alpha_guess, beta_guess = alpha_beta_from_E_L(E,L,nbu)
            nbStep += 1

            if (nbStep > nbStepMax || (abs(stepE)<abs(_E0*eps) && abs(stepL)<abs(_L0*eps)) )
                break
            end

        end

        return E, L, alpha - alpha_guess, beta - beta_guess, nbStep
    end
end
