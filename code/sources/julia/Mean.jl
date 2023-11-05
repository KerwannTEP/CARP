##################################################
# Non-dimensional quantities
##################################################

function _tE(E::Float64)
    return E/_E0
end

function _tL(L::Float64)
    return L/_L0
end

function _lz(Lz::Float64, L::Float64)
    return Lz/L
end

function sincosI(Lz::Float64, L::Float64)
    lz = _lz(Lz,L)
    cosI = lz
    sinI = sqrt(1.0 - lz^2)
    return sinI, cosI
end

##################################################
# Coordinates
##################################################

function _xyz(r::Float64, th::Float64, Lz::Float64, L::Float64)
    sinI, cosI = sincosI(Lz,L)
    x = r*cos(th)*cosI
    y = r*sin(th)
    z = r*cos(th)*sinI
    return x, y, z
end

function _vxvyvz(vr::Float64, vt::Float64, th::Float64, Lz::Float64, L::Float64)
    sinI, cosI = sincosI(Lz,L)
    vx = cosI * (vr*cos(th) - vt*sin(th))
    vy = vr*sin(th) + vt*cos(th)
    vz = sinI * (vr*cos(th) - vt*sin(th))
    return vx, vy, vz
end

##################################################
# Potential
##################################################

function psi(r::Float64)
    return _E0/sqrt(1.0+(r/_b)^2)
end

function psiEff(r::Float64, L::Float64)
    return psi(r) + L^2/(2.0*r^2)
end

function dpsidr(r::Float64)
    return -(_E0/_b)*(r/_b)/(1.0+(r/_b)^2)^(3/2)
end

function dpsiEffdr(r::Float64, L::Float64)
    return dpsidr(r) - L^2/r^3
end

##########################################
# Mean field effective potential in s-variable
##########################################

function psieff_s(s::Float64, L::Float64)
    if (L != 0.0)
       return _E0/s + L^2/(2.0*_b^2*(s^2-1.0))
   else
       return _E0/s
   end
end

function dpsieff_s(s::Float64, L::Float64)
    if (L != 0.0)
        ds = -_E0/s^2 - s*L^2/(_b^2*(s^2-1.0)^2)
        dL = L/(_b^2*(s^2-1.0))
        return ds, dL
    else
        ds = -_E0/s^2
        dL = 0.0
        return ds, dL
    end
end

##########################################
##########################################

function vr_from_E_L(r::Float64, E::Float64, L::Float64)
    return sqrt(abs(2.0*(E - psi(r)) - L^2/r^2))
end

function vt_from_E_L(r::Float64, E::Float64, L::Float64)
    return L/r
end

function vrSq_from_E_L(r::Float64, E::Float64, L::Float64)
    return 2.0*(E - psi(r)) - L^2/r^2
end

function vSq_from_E_L(r::Float64, E::Float64, L::Float64)
    return 2.0*(E - psi(r))
end

##################################################
# Compute Lc(E) the angular momentum per unit mass of a circular orbit
##################################################

function _tEta(s::Float64, tE::Float64)
    return (s^2-1)*(1/s - tE)
end

function _sc(tE::Float64)
    t1 = 1.0/(6.0*tE)
    t2 = (1.0+54.0*tE^2)/(216.0*tE^3)
    t3 = (1.0/(4.0*tE))*sqrt(1.0+1.0/(27.0*tE^2))
    return t1+cbrt(t2+t3)+cbrt(t2-t3)
end

function _rc(tE::Float64)
    sc = _sc(tE)
    return _b*sqrt(sc^2-1.0)
end

function _Lc(E::Float64)
    tE = E/_E0
    if (tE == 1.0)
        return 0.0
    else
        return _L0*sqrt(abs(2.0*_tEta(_sc(tE),tE)))
    end
end




function bisection(fun, xl::Float64, xu::Float64, tolx::Float64=4.0*eps(Float64), tolf::Float64=4.0*eps(Float64), iterMAX::Int64=50)
    if (xl > xu)
        xl, xu = xu, xl # Ordering the input arguments
    end
    #####
    fl, fu = fun(xl), fun(xu) # Evaluating the function on the bracket
    #####
    if (abs(fl) <= tolf) # We have already found a solution on the left bracket
        return xl # Returning the left bracket
    end
    #####
    if (abs(fu) <= tolf) # We have already found a solution on the right bracket
        return xu # Returning the right bracket
    end
    #####
    @assert fl*fu < 0.0 (xl,xu,fl,fu)  #"bisection: NOT A BRACKET"
    #####
    iter = 0 # Counter for the iterations
    #####
    while true # Bisection loop
        #####
        xm = (xl+xu)*0.5 # Middle value
        #####
        if ((abs(xu-xl) <= tolx) || (iter > iterMAX)) # The considered bracket is smaller than the tolerance, or we have made too many iterations
            return xm # Returning the middle value
        end
        #####
        fm = fun(xm) # Value at the midpoint
        #####
        iter += 1 # Updating the counter of iterations
        #####
        if (abs(fm) <= tolf) # The middle value is below the threshold
            return xm # Returning the middle value
        end
        #####
        # Otherwise, we iterate the bisection
        if (fm*fl < 0.0) # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
    end
end
