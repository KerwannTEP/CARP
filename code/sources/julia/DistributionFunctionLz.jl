

function F_Lz(Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            Jrmax::Float64=3.0*_L0, nbJr::Int64=30,
            DeltaLmax::Float64=2.0*_L0, nbL::Int64=20)

    nbJrL = nbJr*nbL
    tab_JrL = zeros(Float64,2,nbJrL)

    sum = Threads.Atomic{Float64}(0.0)

    iGrid = 1
    for iJr=1:nbJr, iL=1:nbL
        Jr = Jrmax* (iJr-0.5)/nbJr
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        tab_JrL[1,iGrid], tab_JrL[2,iGrid] = Jr, L
        iGrid += 1
    end

    Threads.@threads for iGrid=1:nbJrL
        Jr = tab_JrL[1,iGrid]
        L = tab_JrL[2,iGrid]
        E = E_from_Jr_L(Jr,L)

        Frot = _Frot(E,L,Lz)
        #dfdt = dFdtOptiExactSign(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,epsEff,m_test)


        Threads.atomic_add!(sum,Frot)

    end

    sum[] *=  (2*pi)^3*Jrmax/nbJr* DeltaLmax/nbL

    return sum[]
end

function dFdt_Lz(Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            Jrmax::Float64=2.0*_L0, nbJr::Int64=20,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-2))

    nbJrL = nbJr*nbL
    tab_JrL = zeros(Float64,2,nbJrL)

    sum = Threads.Atomic{Float64}(0.0)

    iGrid = 1
    for iJr=1:nbJr, iL=1:nbL
        Jr = Jrmax* (iJr-0.5)/nbJr
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        tab_JrL[1,iGrid], tab_JrL[2,iGrid] = Jr, L
        iGrid += 1
    end

    Threads.@threads for iGrid=1:nbJrL
        Jr = tab_JrL[1,iGrid]
        L = tab_JrL[2,iGrid]
        E = E_from_Jr_L(Jr,L)

        #Frot = _Frot(E,L,Lz)
        dfdt, _ = dFdtOptiExactSign(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*L),m_test)


        Threads.atomic_add!(sum,dfdt)

    end

    sum[] *=  (2*pi)^3*Jrmax/nbJr* DeltaLmax/nbL

    return sum[]
end


function dFdtNoRot_Lz(Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            Jrmax::Float64=2.0*_L0, nbJr::Int64=20,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0,
            m_test::Float64=m_field, eps::Float64=10^(-2))

    nbJrL = nbJr*nbL
    tab_JrL = zeros(Float64,2,nbJrL)

    sum = Threads.Atomic{Float64}(0.0)

    iGrid = 1
    for iJr=1:nbJr, iL=1:nbL
        Jr = Jrmax* (iJr-0.5)/nbJr
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        tab_JrL[1,iGrid], tab_JrL[2,iGrid] = Jr, L
        iGrid += 1
    end

    Threads.@threads for iGrid=1:nbJrL
        Jr = tab_JrL[1,iGrid]
        L = tab_JrL[2,iGrid]
        E = E_from_Jr_L(Jr,L)

        #Frot = _Frot(E,L,Lz)
        dfdt = dFdt_NoRot(Jr,L,Lz,m_field,nbK,nbAvr,nbu,m_test,min(eps,0.01*Jr,0.01*L))


        Threads.atomic_add!(sum,dfdt)

    end

    sum[] *=  (2*pi)^3*Jrmax/nbJr* DeltaLmax/nbL

    return sum[]
end


###########################################
# 2D Jr,Lz
###########################################


function dFdt_JrLz(Jr::Float64, Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20, eps::Float64=10^(-2),
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)

    sum = 0.0

    for iL=1:nbL
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        E = E_from_Jr_L(Jr,L)

        #Frot = _Frot(E,L,Lz)
        dfdt, _ = dFdtOptiExactSign(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*L),m_test)


        sum += dfdt

    end

    sum *=  (2*pi)^3*DeltaLmax/nbL

    return sum
end


function dFdt_JrLz_Par(Jr::Float64, Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20, eps::Float64=10^(-2),
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)

    sum = 0.0

    for iL=1:nbL
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        E = E_from_Jr_L(Jr,L)

        #Frot = _Frot(E,L,Lz)
        dfdt = dFdtOptiExactSignPar(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*L),m_test)


        sum += dfdt

    end

    sum *=  (2*pi)^3*DeltaLmax/nbL

    return sum
end

function dFdt_JrLz_DecompPar(Jr::Float64, Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20, eps::Float64=10^(-2),
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)

    sumJr_Jr = 0.0
    sumL_L = 0.0
    sumLz_Lz = 0.0

    for iL=1:nbL
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        E = E_from_Jr_L(Jr,L)

        #Frot = _Frot(E,L,Lz)
        fJr_Jr, fL_L, fLz_Lz = dFdtOptiExactSignParDecomp(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*L),m_test)


        sumJr_Jr += fJr_Jr
        sumL_L += fL_L
        sumLz_Lz += fLz_Lz

    end

    sumJr_Jr *=  (2*pi)^3*DeltaLmax/nbL
    sumL_L *=  (2*pi)^3*DeltaLmax/nbL
    sumLz_Lz *=  (2*pi)^3*DeltaLmax/nbL

    return sumJr_Jr, sumL_L, sumLz_Lz
end


function FluxLz_JrLz_DecompPar(Jr::Float64, Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20, eps::Float64=10^(-2),
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)

    sumLz_Lz = 0.0
    sumJrLz_JrLz = 0.0
    sumLLz_LLz = 0.0
    sumLzLz_LzLz = 0.0
    sumJrLz_Lz = 0.0
    sumLLz_Lz = 0.0

    for iL=1:nbL
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        E = E_from_Jr_L(Jr,L)

        #Frot = _Frot(E,L,Lz)
        dLz_Lz, dJrLz_JrLz,dLLz_LLz,dLzLz_LzLz,dJrLz_Lz,dLLz_Lz = FluxLzOptiExactSignParDecomp(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*L),m_test)


        sumLz_Lz += dLz_Lz
        sumJrLz_JrLz += dJrLz_JrLz
        sumLLz_LLz += dLLz_LLz
        sumLzLz_LzLz += dLzLz_LzLz
        sumJrLz_Lz += dJrLz_Lz
        sumLLz_Lz += dLLz_Lz

    end


    sumLz_Lz *=  (2*pi)^3*DeltaLmax/nbL
    sumJrLz_JrLz *=  (2*pi)^3*DeltaLmax/nbL
    sumLLz_LLz *=  (2*pi)^3*DeltaLmax/nbL
    sumLzLz_LzLz *=  (2*pi)^3*DeltaLmax/nbL
    sumJrLz_Lz *=  (2*pi)^3*DeltaLmax/nbL
    sumLLz_Lz*=  (2*pi)^3*DeltaLmax/nbL

    return sumLz_Lz, sumJrLz_JrLz,sumLLz_LLz ,sumLzLz_LzLz,sumJrLz_Lz ,sumLLz_Lz
end

# Compute flux2D components
# Jr : from integrated diff coeff
# Lz : idem
# returns (dJr, dJrJr, dJrLz, dLz, dLzLz)
function Components_Flux2D_JrLz_Par(Jr::Float64, Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)

    sumJr = 0.0
    sumJrJr = 0.0
    sumJrLz = 0.0
    sumLz = 0.0
    sumLzLz = 0.0

    for iL=1:nbL
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        E = E_from_Jr_L(Jr,L)

        #Frot = _Frot(E,L,Lz)
        dJr, dL, dLz, dJrJr, dLL, dLzLz, dJrL, dJrLz, dLLz = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)
        Frot = _Frot(E,L,Lz,alpha)

        sumJr += dJr*Frot
        sumJrJr += dJrJr*Frot
        sumJrLz += dJrLz*Frot
        sumLz += dLz*Frot
        sumLzLz += dLzLz*Frot

    end

    sumJr *=  (2*pi)^3*DeltaLmax/nbL
    sumJrJr *=  (2*pi)^3*DeltaLmax/nbL
    sumJrLz *=  (2*pi)^3*DeltaLmax/nbL
    sumLz *=  (2*pi)^3*DeltaLmax/nbL
    sumLzLz *=  (2*pi)^3*DeltaLmax/nbL


    return sumJr, sumJrJr, sumJrLz, sumLz, sumLzLz
end



# integrate flux components over L
function Flux_JrLz(Jr::Float64, Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20, eps::Float64=10^(-2),
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)

    integratefluxJr = 0.0
    integratefluxL = 0.0
    integratefluxLz = 0.0
    integratefluxLzBound = 0.0

    for iL=1:nbL
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        E = E_from_Jr_L(Jr,L)

        #Frot = _Frot(E,L,Lz)
        fJr, fL, fLz = fluxOptiExactSign(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*L),m_test)

        integratefluxJr += fJr
        integratefluxLz += fLz



    end

    # boundary
    # Check at Lz = 0.999 L instead of L
    # To get the correct approximate value at |Lz|=L
    fJr, fL, fLz = fluxOptiExactSign(Jr,abs(Lz),0.999*Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*abs(Lz)),m_test)
    integratefluxL = -fL
    integratefluxLzBound = sign(Lz)*fLz


    integratefluxJr *=  (2*pi)^3*DeltaLmax/nbL
    integratefluxL *=  (2*pi)^3
    integratefluxLz *=  (2*pi)^3*DeltaLmax/nbL
    integratefluxLzBound *=  (2*pi)^3


    return integratefluxJr, integratefluxL, integratefluxLz, integratefluxLzBound
end



# integrate flux components over L
function Flux_JrLz_Par(Jr::Float64, Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20, eps::Float64=10^(-2),
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)

    integratefluxJr = 0.0
    integratefluxL = 0.0
    integratefluxLz = 0.0
    integratefluxLzBound = 0.0

    for iL=1:nbL
        L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
        E = E_from_Jr_L(Jr,L)

        #Frot = _Frot(E,L,Lz)
        fJr, fL, fLz = fluxOptiExactSignPar(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*L),m_test)

        integratefluxJr += fJr
        integratefluxLz += fLz



    end

    # boundary
    # Check at Lz = 0.999 L instead of L
    # To get the correct approximate value at |Lz|=L
    fJr, fL, fLz = fluxOptiExactSignPar(Jr,abs(Lz),0.999*Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*abs(Lz)),m_test)
    integratefluxL = -fL
    integratefluxLzBound = sign(Lz)*fLz


    integratefluxJr *=  (2*pi)^3*DeltaLmax/nbL
    integratefluxL *=  (2*pi)^3
    integratefluxLz *=  (2*pi)^3*DeltaLmax/nbL
    integratefluxLzBound *=  (2*pi)^3


    return integratefluxJr, integratefluxL, integratefluxLz, integratefluxLzBound
end


# integrate flux components over Jr,L
function Flux_Lz_Par(Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20,
            JrLmax::Float64=1.5*_L0, nbJr::Int64=20, eps::Float64=10^(-2),
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)


    integratefluxLz = 0.0

    for iJr=1:nbJr
        Jr = Jrmax*(iJr-0.5)/nbJr
        for iL=1:nbL
            L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
            E = E_from_Jr_L(Jr,L)

            #Frot = _Frot(E,L,Lz)
            fJr, fL, fLz = fluxOptiExactSignPar(Jr,L,Lz,m_field,alpha,nbAvr,nbK,nbu,min(eps,0.01*Jr,0.01*L),m_test)
            integratefluxLz += fLz
        end
    end

    integratefluxLz *=  (2*pi)^3*DeltaLmax/nbL*Jrmax/nbJr



    return integratefluxLz
end

# Return Dz*Fz, Dzz*Fz
function DiffCoeff_Lz_Par(Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20,
            JrLmax::Float64=1.5*_L0, nbJr::Int64=20,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)


    Dz = 0.0
    Dzz = 0.0
    Fz = 0.0

    for iJr=1:nbJr
        Jr = Jrmax*(iJr-0.5)/nbJr
        for iL=1:nbL
            L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
            E = E_from_Jr_L(Jr,L)

            Frot = _Frot(E,L,Lz,alpha)
            avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz = orbitAverageActionCoeffsOptiExactPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbK,nbu,m_test)

            Dz += avrDLz*Frot
            Dzz += avrDLzLz*Frot

            Fz += Frot

        end
    end

    Dz *=  (2*pi)^3*DeltaLmax/nbL*Jrmax/nbJr
    Dzz *=  (2*pi)^3*DeltaLmax/nbL*Jrmax/nbJr

    Fz *=  (2*pi)^3*DeltaLmax/nbL*Jrmax/nbJr

    return Dz, Dzz, Fz
end


# Return Dz*Fz, Dzz*Fz
function DiffCoeff_Lz_SmoothPar(Lz::Float64, m_field::Float64, alpha::Float64=alphaRot,
            DeltaLmax::Float64=1.5*_L0, nbL::Int64=20,
            JrLmax::Float64=1.5*_L0, nbJr::Int64=20,
            nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default, nbAvrTh::Int64=nbAvrTh_default,
            nbu::Int64 = nbu0, m_test::Float64=m_field)


    Dz = 0.0
    Dzz = 0.0
    Fz = 0.0

    for iJr=1:nbJr
        Jr = Jrmax*(iJr-0.5)/nbJr
        for iL=1:nbL
            L = abs(Lz) + DeltaLmax*(iL-0.5)/nbL
            E = E_from_Jr_L(Jr,L)

            Frot = _Frot(E,L,Lz,alpha)
            avrDJr, avrDL, avrDLz, avrDJrJr, avrDLL, avrDLzLz, avrDJrL, avrDJrLz, avrDLLz = orbitAverageActionCoeffsOptiPar(Jr,L,Lz/L,m_field,alpha,nbAvr,nbAvrTh,nbK,nbu,m_test)

            Dz += avrDLz*Frot
            Dzz += avrDLzLz*Frot

            Fz += Frot

        end
    end

    Dz *=  (2*pi)^3*DeltaLmax/nbL*Jrmax/nbJr
    Dzz *=  (2*pi)^3*DeltaLmax/nbL*Jrmax/nbJr

    Fz *=  (2*pi)^3*DeltaLmax/nbL*Jrmax/nbJr

    return Dz, Dzz, Fz
end
