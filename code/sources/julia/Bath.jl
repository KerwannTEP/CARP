using HypergeometricFunctions
using SpecialFunctions
using StaticArrays # To have access to static arrays
using Interpolations # To have access to interpolation functions

##################################################
# Distribution function in (E,L) for a Plummer sphere
##################################################

function _Hq1(x::Float64, q::Float64=qCalc)
    if (x <= 1)
        pref = 1/(GAMMA_ca)
        HG   = H_1(x)
        return pref*HG
    else
        pref = 1/(GAMMA_db*GAMMA_bc)
        HG   = H_2(1/x)
        return pref*1/sqrt(x)*HG
    end
end

function _Hqm6(x::Float64, q::Float64=qCalc)
    if (x <= 1)
        pref = 1/(GAMMA_ca)
        HG   = H_1(x)
        return pref*HG
    else
        pref = 1/(GAMMA_db*GAMMA_bc)
        HG   = H_2(1/x)
        return pref*x^3*HG
    end
end

function _Hq(x::Float64, q::Float64=qCalc)
    if (x <= 1)
        pref = 1/(GAMMA_ca)
        HG   = H_1(x)
        return pref*HG
    else
        pref = 1/(GAMMA_db*GAMMA_bc)
        HG   = H_2(1/x)
        return pref*x^(-q/2)*HG
    end
end

if (qCalc == 1.0)
    const _H = _Hq1
elseif (qCalc == -6.0)
    const _H = _Hqm6
else
    const _H = _Hq
end


function _tFq0(tE::Float64, tL::Float64)
    if (tE < 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    end
    return 3.0/(7.0*PI^3) * (2.0*tE)^(3)*sqrt(2.0*tE)

end

function _tFq2(tE::Float64, tL::Float64)
    if (tE < 0.0 || tL < 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    end
    # If E and L are positive
    x = tL^2/(2.0*tE)
    if (x <= 1)
        return 6.0/(2.0*PI)^3 * (2.0*tE - tL^2)^(3/2)
    end

    return 0.0

end

function _tFq1(tE::Float64, tL::Float64)
    if (tE < 0.0 || tL < 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    end
    # If E and L are positive
    x = tL^2/(2.0*tE)
    return (3.0*GAMMA_6q/(2.0*(2.0*PI)^(5/2)) *
               tE*tE*sqrt(tE) * _H(x,1.0))


end

function _tFqm6(tE::Float64, tL::Float64)
    if (tE < 0.0 || tL < 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    end
    # If E and L are positive
    x = tL^2/(2.0*tE)
    return (3.0*GAMMA_6q/(2.0*(2.0*PI)^(5/2)) *
               tE^2*tE^2*tE^2*tE^2*tE*sqrt(tE) * _H(x,-6.0))

end

function _tFq(tE::Float64, tL::Float64)
    if (tE < 0.0 || tL < 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    end
    # If E and L are positive
    x = tL^2/(2.0*tE)
    return (3.0*GAMMA_6q/(2.0*(2.0*PI)^(5/2)) *
               tE^(3.5-qCalc) * _H(x,qCalc))

end


if (qCalc == 0.0)
    const _tF = _tFq0
elseif (qCalc == 2.0)
    const _tF = _tFq2
elseif (qCalc == 1.0)
    const _tF = _tFq1
elseif (qCalc == -6.0)
    const _tF = _tFqm6
else
    const _tF = _tFq
end

function _F(E::Float64, L::Float64)
    tE = _tE(E)
    tL = _tL(L)
    DF  = _tF(tE,tL)
    return _M*_F0*DF
end


##################################################
# Adding rotation to the DF
##################################################


function _Frot(E::Float64, L::Float64, Lz::Float64, alpha::Float64=alphaRot)
    Ftot = _F(E,L)
    Frot = Ftot*(1.0 + alpha*sign(Lz/L))
    return Frot
end

# Normalized to M = int dJr dL dcosI _Frot_cosI(Jr,L,cosI)
function _Frot_cosI(E::Float64, L::Float64, cosI::Float64, alpha::Float64=alphaRot)
    Ftot = _F(E,L)
    Frot = L*Ftot*(1.0 + alpha*sign(cosI))
    return Frot
end

##################################################
# Interpolate the hypergeometric function
# HG1   = _₂F₁(qCalc/2,qCalc-3.5,1.0,x)
# HG2   = _₂F₁(qCalc/2,qCalc/2,4.5-qCalc/2,1/x)

# valeur en x=1: https://homepage.tudelft.nl/11r49/documents/wi4006/hyper.pdf equation 4
# 2F1(a,b,c,1) = gamma(c)gamma(c-a-b)/(gamma(c-a)gamma(c-b))

# Valeur en x=1
# HG1   = _₂F₁(qCalc/2,qCalc-3.5,1.0,1) = gamma(4.5-1.5*qCalc)/(gamma(1-qCalc/2)gamma(4.5-qCalc))
# HG2   = _₂F₁(qCalc/2,qCalc/2,4.5-qCalc/2,1) = gamma(4.5-qCalc/2)gamma(4.5-1.5*qCalc)/(gamma(4.5-qCalc)gamma(4.5-qCalc))
##################################################

function getHGInt(nbxInt::Int64=1000)
    xminInt = 0.0
    xmaxInt = 1.0
    #####
    rangexInt = range(xminInt,length=nbxInt,xmaxInt)
    tabxInt = collect(rangexInt)
    tabHG1Int = zeros(Float64,nbxInt)
    tabHG2Int = zeros(Float64,nbxInt)

    #####
    for indx=2:nbxInt-1
        xloc = tabxInt[indx]
        hg1loc = _₂F₁(qCalc/2,qCalc-3.5,1.0,xloc)
        hg2loc = _₂F₁(qCalc/2,qCalc/2,4.5-qCalc/2,xloc)
        tabHG1Int[indx] = hg1loc
        tabHG2Int[indx] = hg2loc
    end
    # x=0
    tabHG1Int[1] = 1.0
    tabHG2Int[1] = 1.0

    # x=1
    tabHG1Int[nbxInt] = gamma(4.5-1.5*qCalc)/(gamma(1-qCalc/2)*gamma(4.5-qCalc))
    tabHG2Int[nbxInt] = gamma(4.5-qCalc/2)*gamma(4.5-1.5*qCalc)/(gamma(4.5-qCalc)*gamma(4.5-qCalc))
    #####
    intHG1 = Interpolations.scale(interpolate(tabHG1Int, BSpline(Cubic(Line(OnGrid())))),rangexInt)
    intHG2 = Interpolations.scale(interpolate(tabHG2Int, BSpline(Cubic(Line(OnGrid())))),rangexInt)
    #####
    return [intHG1, intHG2]
end


const tabHyperGeoInt = SVector{2}(getHGInt())

# HG1   = _₂F₁(qCalc/2,qCalc-3.5,1.0,x)`
const H_1 = tabHyperGeoInt[1]



# HG2   = _₂F₁(qCalc/2,qCalc/2,4.5-qCalc/2,1/x)`
const H_2 = tabHyperGeoInt[2]
