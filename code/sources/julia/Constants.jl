##################################################
# Useful constants
##################################################

using SpecialFunctions

# Numerical Constants
# Calculation done with with Julia

"""
    PI

`π` constant, set equal to `3.1415926535897`.
Computed by `Julia`.
"""
const PI = 3.1415926535897

##################################################
# Constants of the problem
##################################################

"""
    nbGlobularCluster

Number of stars within the Plummer cluster.
"""
const nbGlobularCluster = 10^5


"""
    parameterCoulomb

Coulomb logarithm of the cluster.
Set as `ln Λ = ln(0.11 N)` where `N` is the number of stars in the cluster.

See `Heggie et Hut (2003)`.
"""
const logCoulomb = log(0.11*nbGlobularCluster)

"""
    _G

Newton's gravitational constant.
"""
const _G = 1.0

"""
    _M

Total mass of the globular cluster.
"""
const _M = 1.0

"""
    _b

Plummer radius of the globular cluster.
"""
const _b = 1.0


##################################################
# Dimensional units (Plummer potential)
##################################################

"""
    _v0

Velocity unit of the Plummer cluster.
"""
const _v0 = sqrt(_G*_M/_b)

"""
    _E0

Energy unit of the Plummer cluster.
"""
const _E0 = -_G*_M/_b

"""
    _L0

Action unit of the Plummer cluster.
"""
const _L0 = sqrt(_G*_M*_b)

"""
    _F0

Distribution function unit of the Plummer cluster.
"""
const _F0 = (_G*_M*_b)^(-3/2)

"""
    _Omega0

Frequency unit of the Plummer cluster.
"""
const _Omega0 = sqrt(_G*_M/_b^3)

"""
    _rho0

Density unit of the Plummer cluster.
"""
const _rho0 = _M/(4*PI*_b^3/3)

"""
    m_field

Individual mass of the globular cluster's stars
"""
const m_field = _M/nbGlobularCluster


##################################################
# Constant used to compute H(x,q) and the DF
##################################################

"""
    GAMMA_ca

`gamma(c-a) = gamma(9/2 - qCalc)`.
"""
const GAMMA_ca = gamma(9/2 - qCalc)

"""
    GAMMA_db

`gamma(d-b) = gamma(1 - qCalc/2)`.
"""
const GAMMA_db = gamma(1 - qCalc/2)

"""
    GAMMA_bc

`gamma(b-c) = gamma(9/2 - qCalc/2)`.
"""
const GAMMA_bc = gamma(9/2 - qCalc/2)

"""
    GAMMA_6q

`gamma(6-q) = gamma(6-qCalc)`.
"""
const GAMMA_6q = gamma(6-qCalc)