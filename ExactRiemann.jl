module ExactRiemann

export Exact_Riemann

struct Const
    Y0 ::Float64
    ρ0 ::Float64
    Γ0 ::Float64
    μ  ::Float64
    a0 ::Float64
    s0 ::Float64
end
    

struct Var
    ρ ::Float64
    u ::Float64
    p ::Float64
    sxx::Float64
end


struct Shock
    cs::Int 
end

struct Rare
    cs::Int
    
end
###  cs=1 One elastic wave; cs=2 One plastic wave ; cs=3 One plastic wave and one plastic wave

function CaseSelect(sxx₀::Float64,ρ₀::Float64, ρ₁::Float64,con::Const)
    Y0 = con.Y0
     sxx₁ =  sxx₀ -4/3*con.μ *(log(ρ₁) - log(ρ₀))
       
        if ρ₁ > ρ₀   #Shock
                if sxx₀ ≤ -2/3*Y0            
                    case = Shock(2) # Sᵖ
                elseif sxx₁ ≤ -2/3*Y0
                    case = Shock(3) # SᵉSᵖ
                else
                    case = Shock(1) # Sᵉ
                end
            else   #Rarefaction
                if sxx₀ ≥ 2/3*Y0
                    case = Rare(2) #Rᵖ
                elseif sxx₁ ≥  2/3*Y0
                    case = Rare(3) # RᵉRᵖ
                else
                    case = Rare(1) #Rᵉ
                end
            end
    return case
end

function Exact_Riemann(UL::Array{Float64,1}, UR::Array{Float64,1}, conL::Const, conR::Const )
        
    
    var₀L = Var(UL[1],UL[2],UL[3],UL[4])
    var₀R = Var(UR[1],UR[2],UR[3],UR[4])

    TOL = 2.e-4
    i = 1
    U = zeros(Float64, 2)
    
    ρL, uL, pL, sxxL = UL[1:4]
    ρR, uR, pR, sxxR = UR[1:4]
    ρLStar = ρL
    ρRStar = ρR
    varL = var₀L
    varR = var₀R
    caseL = Shock(1)
    caseR = Shock(1)
    
    while true  #Iteration  
        
        caseL.cs == 3 ?  varL = Ṽar(var₀L,caseL,conL,1) : varL = var₀L
        caseR.cs == 3 ?  varR = Ṽar(var₀R,caseR,conR,2) : varR = var₀R
        
        J = [[-fR¹(ρRStar,varR, conR, caseR) fL¹(ρLStar,varL, conL, caseL)]; 
             [-gR¹(ρRStar,varR, conR, caseR) gL¹(ρLStar,varL, conL, caseL)]]
    
        U = [ρRStar, ρLStar]
    
        F = [fL(ρLStar,varL, conL, caseL)- fR(ρRStar,varR, conR, caseR), 
             gL(ρLStar,varL, conL, caseL)- gR(ρRStar,varR, conR, caseR)]
        @show UL,UR
     c = max(abs(F[1]),abs(F[2])) 
    
      if c <= TOL
            break
        end
    #   @show F
        U = U - J\F
        ρRStar = U[1]
        ρLStar = U[2]
        
        caseL = CaseSelect(sxxL,ρL,ρLStar,conL)
        caseR = CaseSelect(sxxR,ρR,ρRStar,conR)
   end
   
    sxxLStar = Sxx(ρLStar, varL, conL)
    sxxRStar = Sxx(ρRStar, varR, conR)
    
    ### 重新求解
    
    pLStar = fp(ρLStar, varL, conL,caseL)
    s_Star = fu(ρLStar, varL, conL,caseL,1)
    pRStar = fp(ρRStar, varR, conR,caseR)    
    ŨL = [varL.ρ, varL.u, varL.p, varL.sxx]
    ŨR = [varR.ρ, varR.u, varR.p, varR.sxx]
    ULStar = [ρLStar, s_Star, pLStar,sxxLStar]
    URStar = [ρRStar, s_Star, pRStar,sxxRStar]
    @show caseL,caseR      
    return  ŨL, ULStar,URStar, ŨR, caseL,caseR  
end

function Ṽar(var0::Var, case, con::Const,LoR::Int)

    Y0, μ= con.Y0, con.μ
    sxx,ρ = var0.sxx,var0.ρ

        if typeof(case) == Rare 
            sxx1 = 2/3*Y0
            ρ1 = ρ*exp(-Y0/(2μ)+(3sxx)/(4μ))
        elseif typeof(case) == Shock
            sxx1 = -2/3*Y0
            ρ1 = ρ*exp(Y0/(2μ)+(3sxx)/(4μ))
        end
        p1 = fp(ρ1, var0, con, case)
        u1 = fu(ρ1, var0, con, case,LoR)
        var1= Var(ρ1, u1, p1, sxx1)
    return var1
end

function Sxx(ρ, var0::Var, con::Const)
    sxx0,ρ0 = var0.sxx,var0.ρ
    Y0   = con.Y0
    
    sxx = sxx0 -4/3*con.μ *(log(ρ) - log(ρ0))
    if abs(sxx) ≥ 2/3*Y0
        sxx = 2/3*Y0*sign(sxx)
    end
    return sxx
end

function fp(ρ::Float64, var0::Var, con::Const, case::Rare)
    sxx₀,ρ₀,p₀ =var0.sxx, var0.ρ,var0.p
    ρ0,Γ0,a0,Y0= con.ρ0,con.Γ0,con.a0,con.Y0 
    λ₁ = ρ0*Γ0

            
           p = p₀*exp(λ₁/ρ₀ - λ₁/ρ) + exp(-λ₁/ρ)*GaussIntegral(ρ->
           f3(ρ,con,var0,λ₁,a0)*exp(λ₁/ρ), ρ₀, ρ,5)

           return p 
    
end
  
f3(ρ::Float64,con::Const,var0::Var,λ₁::Float64,a0::Float64) =a0^2*fηη(ρ, con) -λ₁*Sxx(ρ, var0,con)/ρ^2

function fu(ρ::Float64, var0::Var, con::Const,case::Rare,LoR::Int)
    ρ₀,u₀ = var0.ρ, var0.u
    
    
    if LoR == 1 
        u = u₀ - GaussIntegral(ρ->f2(ρ,con,var0,case), ρ₀, ρ, 7)
    else
         u = u₀ + GaussIntegral(ρ->f2(ρ,con,var0,case), ρ₀, ρ, 7)   
    end
    
    return u
end
function f2(ρ::Float64,con::Const,var0::Var,case::Rare)
    a0,ρ0,Γ0,μ,Y0= con.a0, con.ρ0, con.Γ0, con.μ,con.Y0
    λ₁ = ρ0*Γ0
    sxx₀ = var0.sxx
    S2 = f3(ρ,con,var0,λ₁,a0) +λ₁*fp(ρ,var0,con,case)/ρ^2
 #   @show sxx₀, Y0
    if sxx₀ ≥ 2/3*Y0   ####!!!
       f2 = √(S2)/ρ
    #    @show f2
        return f2
    else 
        f2 = √(S2 + 4μ/(3ρ))/ρ
        
        return f2
    end  
end

function GaussIntegral(f::Function,x₀::Float64,x₁::Float64,order::Int)
    t₁= (x₁-x₀)/2
    t₂= (x₁+x₀)/2
    ω = zeros(Float64, 5)
    p = zeros(Float64, 5)
    
    if order == 1
        ω[1] = 2.0
        p[1] = 0.0
    elseif order == 3
        ω[1] = 1.0; ω[2] = 1.0
        p[1] = 1/√3.0; p[2] = -1/√3.0
    elseif order == 5
        ω[1] = 8.0/9; ω[2] = 5.0/9; ω[3] = 5.0/9
        p[1] = 0.0; p[2] = -√(3.0/5); p[3] = √(3.0/5)
    elseif order == 7
        ω[1] = (18+√30)/36; ω[2] = (18+√30)/36
        ω[3] = (18-√30)/36; ω[4] = (18-√30)/36
        p[1] = √(3/7-2/7*√(6/5)); p[1] = -√(3/7-2/7*√(6/5))
        p[3] = √(3/7+2/7*√(6/5)); p[4] = -√(3/7+2/7*√(6/5))
    end
    ∑ =sum( t₁*ω[i]*f(t₁*p[i]+t₂) for i in 1:floor(Int,order/2)+1)

    return ∑
end

function fp(ρ::Float64, var0::Var, con::Const,case::Shock)
   
    ρ₀,p₀,sxx₀ = var0.ρ,var0.p,var0.sxx
    ρ0,Γ0,a0 = con.ρ0,con.Γ0,con.a0
    
    σ₀ = -p₀+sxx₀
    e₀ =  PToe(ρ₀, p₀, con)
   # @show e₀
    t =  (ρ-ρ₀)/(ρ₀*ρ) ; c₀ = 1/(Γ0*ρ0); c₁= a0^2/Γ0
  
    p = (2(c₁*fη(ρ,con) + e₀) -t*(σ₀+Sxx(ρ,var0,con)))/(-t+2c₀)
 
    return p
end

function fu(ρ::Float64, var0::Var, con::Const,case::Shock,LoR::Int)
   
    ρ₀,u₀,p₀,sxx₀ = var0.ρ,var0.u,var0.p,var0.sxx
    σ₀ = -p₀+sxx₀
    
    sxx=Sxx(ρ,var0,con)
    p  = fp(ρ, var0, con, case)
    σ = -p+sxx
    
    t =  (ρ-ρ₀)/(ρ₀*ρ)
    if LoR==1
        u = u₀ - √(t*(σ₀ - σ))
    
    else
        u = u₀ + √(t*(σ₀ - σ))
    end
#    @show t, σ₀ - σ
    return u
end

function Sxx(ρ, var0::Var, con::Const)
    sxx0,ρ0 = var0.sxx,var0.ρ
    Y0   = con.Y0
    
    sxx = sxx0 -4/3*con.μ *(log(ρ) - log(ρ0))
    if abs(sxx) ≥ 2/3*Y0
        sxx = 2/3*Y0*sign(sxx)
    end
    return sxx
end

function fp(ρ::Float64, var0::Var, con::Const, case::Rare)
    sxx₀,ρ₀,p₀ =var0.sxx, var0.ρ,var0.p
    ρ0,Γ0,a0,Y0= con.ρ0,con.Γ0,con.a0,con.Y0 
    λ₁ = ρ0*Γ0

            
           p = p₀*exp(λ₁/ρ₀ - λ₁/ρ) + exp(-λ₁/ρ)*GaussIntegral(ρ->
           f3(ρ,con,var0,λ₁,a0)*exp(λ₁/ρ), ρ₀, ρ,5)

           return p 
    
end
  
f3(ρ::Float64,con::Const,var0::Var,λ₁::Float64,a0::Float64) =a0^2*fηη(ρ, con) -λ₁*Sxx(ρ, var0,con)/ρ^2

function fu(ρ::Float64, var0::Var, con::Const,case::Rare,LoR::Int)
    ρ₀,u₀ = var0.ρ, var0.u
    
    
    if LoR == 1 
        u = u₀ - GaussIntegral(ρ->f2(ρ,con,var0,case), ρ₀, ρ, 7)
    else
         u = u₀ + GaussIntegral(ρ->f2(ρ,con,var0,case), ρ₀, ρ, 7)   
    end
    
    return u
end
function f2(ρ::Float64,con::Const,var0::Var,case::Rare)
    a0,ρ0,Γ0,μ,Y0= con.a0, con.ρ0, con.Γ0, con.μ,con.Y0
    λ₁ = ρ0*Γ0
    sxx₀ = var0.sxx
    S2 = f3(ρ,con,var0,λ₁,a0) +λ₁*fp(ρ,var0,con,case)/ρ^2
 #   @show sxx₀, Y0
    if sxx₀ ≥ 2/3*Y0   ####!!!
       f2 = √(S2)/ρ
    #    @show f2
        return f2
    else 
        f2 = √(S2 + 4μ/(3ρ))/ρ
        
        return f2
    end
  
end

function GaussIntegral(f::Function,x₀::Float64,x₁::Float64,order::Int)
    t₁= (x₁-x₀)/2
    t₂= (x₁+x₀)/2
    ω = zeros(Float64, 5)
    p = zeros(Float64, 5)
    
    if order == 1
        ω[1] = 2.0
        p[1] = 0.0
    elseif order == 3
        ω[1] = 1.0; ω[2] = 1.0
        p[1] = 1/√3.0; p[2] = -1/√3.0
    elseif order == 5
        ω[1] = 8.0/9; ω[2] = 5.0/9; ω[3] = 5.0/9
        p[1] = 0.0; p[2] = -√(3.0/5); p[3] = √(3.0/5)
    elseif order == 7
        ω[1] = (18+√30)/36; ω[2] = (18+√30)/36
        ω[3] = (18-√30)/36; ω[4] = (18-√30)/36
        p[1] = √(3/7-2/7*√(6/5)); p[1] = -√(3/7-2/7*√(6/5))
        p[3] = √(3/7+2/7*√(6/5)); p[4] = -√(3/7+2/7*√(6/5))
    end
    ∑ =sum( t₁*ω[i]*f(t₁*p[i]+t₂) for i in 1:floor(Int,order/2)+1)

    return ∑
end

function Ṽar(var0::Var, case, con::Const,LoR::Int)

    Y0, μ= con.Y0, con.μ
    sxx,ρ = var0.sxx,var0.ρ

        if typeof(case) == Rare 
            sxx1 = 2/3*Y0
            ρ1 = ρ*exp(-Y0/(2μ)+(3sxx)/(4μ))
        elseif typeof(case) == Shock
            sxx1 = -2/3*Y0
            ρ1 = ρ*exp(Y0/(2μ)+(3sxx)/(4μ))
        end
        p1 = fp(ρ1, var0, con, case)
        u1 = fu(ρ1, var0, con, case,LoR)
        var1= Var(ρ1, u1, p1, sxx1)
    return var1
end


function fp(ρ::Float64, var0::Var, con::Const,case::Shock)
   
    ρ₀,p₀,sxx₀ = var0.ρ,var0.p,var0.sxx
    ρ0,Γ0,a0 = con.ρ0,con.Γ0,con.a0
    
    σ₀ = -p₀+sxx₀
    e₀ =  PToe(ρ₀, p₀, con)
   # @show e₀
    t =  (ρ-ρ₀)/(ρ₀*ρ) ; c₀ = 1/(Γ0*ρ0); c₁= a0^2/Γ0
  
    p = (2(c₁*fη(ρ,con) + e₀) -t*(σ₀+Sxx(ρ,var0,con)))/(-t+2c₀)
 
    return p
end

function fu(ρ::Float64, var0::Var, con::Const,case::Shock,LoR::Int)
   
    ρ₀,u₀,p₀,sxx₀ = var0.ρ,var0.u,var0.p,var0.sxx
    σ₀ = -p₀+sxx₀
    
    sxx=Sxx(ρ,var0,con)
    p  = fp(ρ, var0, con, case)
    σ = -p+sxx
    
    t =  (ρ-ρ₀)/(ρ₀*ρ)
    if LoR==1
        u = u₀ - √(t*(σ₀ - σ))
    
    else
        u = u₀ + √(t*(σ₀ - σ))
    end
#    @show t, σ₀ - σ
    return u
end

function Derivative(f::Function,x::Float64)
    ϵ₀ = 1e-8
    if abs(x) >= ϵ₀
        ϵ= ϵ₀*x
    else
        ϵ =ϵ₀
    end
    f¹(x) = (f(x+ϵ)-f(x))/ϵ
    return f¹(x)
end

    fL(ρ, varL, conL, caseL) = fu(ρ, varL, conL, caseL,1)
    fR(ρ, varR, conR, caseR) = fu(ρ, varR, conR, caseR,2)
    gL(ρ, varL, conL, caseL) = -fp(ρ,varL,conL,caseL) +Sxx( ρ, varL, conL)
    gR(ρ, varR, conR, caseR) = -fp(ρ,varR,conR,caseR) +Sxx( ρ, varR, conR)

    fL¹(ρ, varL, conL, caseL) = Derivative(ρ->fL(ρ, varL, conL, caseL),ρ)
    fR¹(ρ, varR, conR, caseR) = Derivative(ρ->fR(ρ, varR, conR, caseR),ρ)
    gL¹(ρ, varL, conL, caseL) = Derivative(ρ->gL(ρ, varL, conL, caseL),ρ)
    gR¹(ρ, varR, conR, caseR) = Derivative(ρ->gR(ρ, varR, conR, caseR),ρ)
 

function CHA(R0::Array{Float64,1}, R1::Array{Float64,1}, F::Array{Float64,1})
    
    ρₗ₀ = R0[1]; ρᵣ₀ = R0[2]
    ρₗ₁ = R1[1]; ρᵣ₁ = R1[2]
    
    CHA = max(abs(2(ρₗ₁-ρₗ₀)/(ρₗ₁+ρₗ₀)),
              abs(2(ρᵣ₁-ρᵣ₀)/(ρᵣ₁+ρᵣ₀)), F[1], F[2])
    return CHA
end   

function sound(uo::Array{Float64,1},con::Const,EoP::Int=1)
    a0  = con.a0
    ρ0  = con.ρ0
    Γ0  = con.Γ0
    Y0  = con.Y0
    μ   = con.μ
    ρ   = uo[1]
    uu  = uo[2]
    p   = uo[3]
    sxx = uo[4]
    a2  = a0^2*fηη(ρ,con)+p/ρ^2*ρ0*Γ0
  #  c = ftest(a2,ρ0,ρ,Γ0,sxx,μ,a0,fηη(ρ,con),p)
    if EoP == 2
        c=sqrt(a2-ρ0/ρ^2*Γ0*sxx)
        return c
    else
        c=sqrt(a2-ρ0/ρ^2*Γ0*sxx+4.0/3*μ/ρ)
        return c
    end
end

function sound(U::Const,con::Const,EoP::Int = 1)
    a0,ρ0,Γ0,μ,Y0 = con.a0, con.ρ0, con.Γ0, con.μ,con.Y0
    ρ,u,p,sxx    = U.ρ, U.u,U.p,U.sxx
    a2  = a0^2*fηη(ρ,con)+p/ρ^2*ρ0*Γ0 
    if EoP == 2
        c=sqrt(a2-ρ0/ρ^2*Γ0*sxx)
        return c
    else
        c=sqrt(a2-ρ0/ρ^2*Γ0*sxx+4.0/3*μ/ρ)
        return c
    end
end


function PToe(ρ::Float64,p::Float64,con::Const)
    c=con
    ei = (p-c.ρ0*c.a0^2*fη(ρ,c))/(c.ρ0*c.Γ0)
    return ei
end
function EToP(ρ::Float64,ei::Float64,con::Const)
    c=con
    p = c.ρ0*c.Γ0*ei+c.ρ0*c.a0^2*fη(ρ,c)
    return p
end

function fη(ρ::Float64,c::Const)
    η = ρ/c.ρ0
    fη=(η-1.0)*(η-c.Γ0*(η-1.0)/2.0)/(η-c.s0*(η-1))^2
end

function fηη(ρ::Float64,c::Const)
    η = ρ/c.ρ0
    fηη=(η+(c.s0-c.Γ0)*(η-1))/(η-c.s0*(η-1))^3
end

function HalfRiemann(para::Float64, U::Array{Float64,1}, con::Const, uorσ::Int)
   
    if uorσ == 1
        uStar = para
    else 
        σStar = para
    end
    case = Shock(1)
    
    ρ,u,p,sxx = U[1],U[2],U[3],U[4]
    var = Var(ρ,u,p,sxx)
    ρStar = ρ
    var1=var
    TOL = 1.e-4
    
    
    if uorσ==1  #with a given u
        while true
            case.cs == 3 ?  var1= Ṽar(var, case,con,2) : var1=var 
            F = fR(ρStar,var1, con, case) - para
            if abs(F) < TOL
                break
            end
            ρStar = ρStar - F/fR¹(ρStar,var1, con, case)
            case = CaseSelect(sxx, ρ,ρStar, con)
            
              @show case
        end
    else   # with a given σ
        
        while true
            case.cs == 3 ?  var1= Ṽar(var, case,con,2) : var1=var
           
            F = gR(ρStar,var, con, case) - para
            if abs(F) < TOL
                break
            end

            ρStar = ρStar - F/gR¹(ρStar,var, con, case)
            case = CaseSelect(sxx, ρ,ρStar, con)
            @show ρStar,case
        end
    end
    
##presolve 

    pStar = fp(ρStar, var1, con,case)
    uStar = fu(ρStar,var1,con,case,2)
    sxxStar=Sxx(ρStar, var1, con)
    UStar = [ρStar, uStar, pStar, sxxStar]
    
    Ũ = [var1.ρ, var1.u, var1.p, var1.sxx]
    
    return UStar, Ũ, case
end
   

function Solve(t::Float64,UL::Array{Float64,1},UR::Array{Float64,1},conL::Const,conR::Const)  


ŨL,  ULStar,URStar,ŨR, caseL,caseR = Exact_Riemann(UL,UR, conL, conR )

ρL, uL, pL, sxxL = UL
ρR, uR, pR, sxxR = UR

ULVar = Var(ρL, uL, pL, sxxL)
URVar=Var(ρR, uR, pR, sxxR)
    
ρLStar, uLStar, pLStar, sxxLStar = ULStar
ρRStar, uRStar, pRStar, sxxRStar = URStar

ULStarVar = Var(ρLStar, uLStar, pLStar, sxxLStar)
URStarVar = Var(ρRStar, uRStar, pRStar, sxxRStar) 
    
ρ̃L,ũL, p̃L, s̃xxL = ŨL
ρ̃R,ũR, p̃R, s̃xxR = ŨR
    
ŨLVar = Var(ρ̃L,ũL, p̃L, s̃xxL)
ŨRVar = Var(ρ̃R,ũR, p̃R, s̃xxR)

I=600
U = zeros(Float64, (I,4))
U1 = zeros(Float64, (100,4))
x1 = zeros(Float64,100)
x  = zeros(Float64,I)
    i =1
#for x0 in x        
#if ũL*t ≥ x0   -\infty -> u\tildeL*t

    
    if caseL == Shock(1) ||  caseL == Shock(2)
        sL = (ρLStar*uLStar - ρL*uL)/(ρLStar-ρL)
      #  if sL *t ≥ x0
            for i =1:100
            x[i] =(2-0.01*(i-1))*sL*t
            U[i,:] = UL
        end
        for i = 101:300
            x[i] = ((uLStar-sL)/200*(i-100) + sL)*t
            U[i,:] = ULStar
        end
        
        
    elseif caseL == Shock(3)
        sL = (ρLStar*uLStar - ρ̃L*ũL)/(ρLStar-ρ̃L)
        s̃L = (ρL*uL-ρ̃L*ũL)/(ρL-ρ̃L)
        
        for i =1:100
            x[i] =(2-0.01*(i-1))*s̃L*t
            U[i,:] = UL
        end
        
        for i = 101:200
            x[i] = ((-s̃L+sL)/100*(i-100) + s̃L)*t
            U[i,:] =  ŨL
        end
        for i = 201:300
           x[i] = ((uLStar-sL)/100*(i-200) + sL)*t
         U[i,:] = ULStar
        end
       
        
    elseif caseL == Rare(1) || caseL == Rare(2)
         caseL == Rare(1) ? EoP = 1 : EoP =2
        cL =sound(UL, conL,EoP)
        
        cLStar = sound(ULStar,conL,EoP)
        sLSlow = uL-cL
        sLFast = uLStar - cLStar
   
        for i =1:100
            x[i] =(2-0.01*(i-1))*sLSlow*t
            U[i,:] = UL
        end
        

        
        x[101:200],U[101:200,:] = RareSpeedToU(t,ULVar,ULStarVar,conL,1,caseL)
        
        for i =201:300
            x[i] = ((uLStar - sLFast)/100*(i-200) + sLFast)*t 
            U[i,:] = ULStar
        end 
        
        
    elseif caseL == Rare(3)
        
        cL =sound(UL, conL,1)    
        cLStar = sound(ULStar,conL,2)
        c̃L = sound(ŨL,conL,1)
            
        sLSlow = uL-cL
        sLMid  = ũL-c̃L 
        sLFast = uLStar - cLStar
       @show sLSlow, sLMid, sLFast 
        for i =1:100
            x[i] =(2-0.01*(i-1))*sLSlow*t
            U[i,:] = UL
        end
        
        x1,U1= RareSpeedToU(t, ULVar,ŨLVar,conL,1,Rare(1))
        for i= 1:50
            x[100+i] = x1[2i]
            U[100+i,:] = U1[2i,:]
        end
        
        x1,U1= RareSpeedToU(t, ŨLVar,ULStarVar,conL,1,Rare(2))
        for i= 1:50
            x[150+i] = x1[2i]
            U[150+i,:] = U1[2i,:]
        end
        
         for i =201:300
            x[i] = ((uLStar - sLFast)/100*(i-200) + sLFast)*t 
            U[i,:] = ULStar
        end 
        
    end
    
    
    
        
    if caseR == Shock(1) || caseR == Shock(2) 
        sR = (ρRStar*uRStar - ρR*uR)/(ρRStar-ρR)
      #  if sL *t ≥ x0
      
        for i = 301:500
            x[i] = ((-uRStar+sR)/200*(i-300) + uRStar)*t
            U[i,:] = URStar
        end
        for i =501:600
            x[i] =(1+0.01*(i-501))*sR*t
            U[i,:] = UR
        end
        
    elseif caseR == Shock(3)
        
        s̃R = (ρRStar*uRStar - ρ̃R*ũR)/(ρRStar-ρ̃R)
        
        sR = (ρR*uR-ρ̃R*ũR)/(ρR-ρ̃R)
        
        for i =501:600
            x[i] =(1+0.01*(i-501))*sR*t
            U[i,:] = UR
        end
        
        for i = 401:500
            x[i] = ((-s̃R+sR)/100*(i-400) + s̃R)*t
            U[i,:] =  ŨR
        end
        for i = 301:400
           x[i] = ((-uRStar+s̃R)/100*(i-300) + uRStar)*t
         U[i,:] = URStar
        end
       @show sR ,s̃R
        
    elseif caseR == Rare(1) || caseR == Rare(2)
        
            caseR ==Rare(1) ? EoP = 1 : EoP =2
        cR =sound(UR, conR,EoP)
        
        cRStar = sound(URStar,conR,EoP)
        sRFast = uR + cR
        sRSlow = uRStar + cRStar
   
        for i =501:600
            x[i] =(1+0.01*(i-501))*sRFast*t
            U[i,:] = UR
        end
        
        
        x[401:500],U[401:500,:] = RareSpeedToU(t,URStarVar,URVar,conR,2,caseR)
        
        for i =301:400
            x[i] = ((-uRStar + sRSlow)/100*(i-300) + uRStar)*t 
            U[i,:] = URStar
        end 
        
        
    elseif caseR == Rare(3)
        
        cR =sound(UR, conR,1)    
        cRStar = sound(URStar,conR,2)
        c̃R = sound(ŨR,conR,2)
            
        sRFast = uR+cR
        sRMid  = ũR+c̃R 
        sRSlow = uRStar + cRStar
        
        for i =501:600
            x[i] =(1+0.01*(i-501))*sRFast*t
            U[i,:] = UR
        end
        for i =301:400
            x[i] = ((-uRStar + sRSlow)/100*(i-300) + uRStar)*t 
            U[i,:] = URStar
        end 
        
        x1,U1= RareSpeedToU(t, URStarVar,ŨRVar,conR,2,Rare(2))  # plastic 
        for i= 1:50
            x[400+i] = x1[2i]
            U[400+i,:] = U1[2i,:]
        end 
        
        x1,U1= RareSpeedToU(t,ŨRVar, URVar,conR,2,Rare(1))      #elastic
        for i= 1:50
            x[450+i] = x1[2i]
            U[450+i,:] = U1[2i,:]
        end
                      
    end
    
return x,U
end       
  

function Solve(t::Float64,para::Float64,UR::Array{Float64,1},conR::Const,uorσ::Int)
   
    URStar, ŨR, caseR = HalfRiemann(para,UR, conR,uorσ)

ρR, uR, pR, sxxR = UR
URVar=Var(ρR, uR, pR, sxxR)
    
ρRStar, uRStar, pRStar, sxxRStar = URStar
URStarVar = Var(ρRStar, uRStar, pRStar, sxxRStar) 
    
ρ̃R,ũR, p̃R, s̃xxR = ŨR
ŨRVar = Var(ρ̃R,ũR, p̃R, s̃xxR)
 
I=300
U = zeros(Float64, (I,4))
U1 = zeros(Float64, (100,4))
x1 = zeros(Float64,100)
x  = zeros(Float64,I)
    i =1
    

    
if caseR == Shock(1) || caseR == Shock(2)
        sR = (ρRStar*uRStar - ρR*uR)/(ρRStar-ρR)
      #  if sL *t ≥ x0
    
        for i = 1:200
            x[i] = ((-uRStar+sR)/200*(i) + uRStar)*t
            U[i,:] = URStar
        end
        for i =201:300
            x[i] =(1+0.01*(i-201))*sR*t
            U[i,:] = UR
        end
        
    elseif caseR == Shock(3)
        
        s̃R = (ρRStar*uRStar - ρ̃R*ũR)/(ρRStar-ρ̃R)
        
        sR = (ρR*uR-ρ̃R*ũR)/(ρR-ρ̃R)
        
        for i =201:300
            x[i] =(1+0.01*(i-201))*sR*t
            U[i,:] = UR
        end
        
        for i = 101:200
            x[i] = ((-s̃R+sR)/100*(i-100) + s̃R)*t
            U[i,:] =  ŨR
        end
        for i = 1:100
           x[i] = ((-uRStar+s̃R)/100*(i) + uRStar)*t
         U[i,:] = URStar
        end
       @show sR ,s̃R
        
    elseif caseR == Rare(1)|| caseR == Rare(2)
        
        cR =sound(UR, conR,2)
        
        cRStar = sound(URStar,conR,2)
        sRFast = uR + cR
        sRSlow = uRStar + cRStar
   
        for i =201:300
            x[i] =(1+0.01*(i-201))*sRFast*t
            U[i,:] = UR
        end
        
        
        x[101:200],U[101:200,:] = RareSpeedToU(t,URStarVar,URVar,conR,2,caseR)
        
        for i =1:100
            x[i] = ((-uRStar + sRSlow)/100*(i) + uRStar)*t 
            U[i,:] = URStar
        end 
        
        
    elseif caseR == Rare(3)
        
        cR =sound(UR, conR,1)    
        cRStar = sound(URStar,conR,2)
        c̃R = sound(ŨR,conR,2)
            
        sRFast = uR+cR
        sRMid  = ũR+c̃R 
        sRSlow = uRStar + cRStar
        
        for i =201:300
            x[i] =(1+0.01*(i-201))*sRFast*t
            U[i,:] = UR
        end
        for i =1:100
            x[i] = ((-uRStar + sRSlow)/100*i + uRStar)*t 
            U[i,:] = URStar
        end 
        
        x1,U1= RareSpeedToU(t, URStarVar,ŨRVar,conR,2,Rare(2)) #plastic
        for i= 1:50
            x[100+i] = x1[2i]
            U[100+i,:] = U1[2i,:]
        end 
        
        x1,U1= RareSpeedToU(t,ŨRVar, URVar,conR,2,Rare(1))  # elastic 
        for i= 1:50
            x[150+i] = x1[2i]
            U[150+i,:] = U1[2i,:]
        end
     end
 

    
return x,U
end     
    

function RareSpeedToU(t::Float64, U0::Var, UStar::Var,con ::Const,LoR::Int, case::Rare)
    
    ρ1,u1,p1,sxx1 = U0.ρ, U0.u,U0.p,U0.sxx
    ρ2,u2,p2,sxx2 = UStar.ρ, UStar.u,UStar.p,UStar.sxx

    U = zeros(Float64, (100,4))
    x = zeros(Float64,100)
    EoP =1
    
    Δρ = (ρ2 -ρ1)/100
    
    for i in 1:100

        ρ = ρ1+Δρ*(i-1)
       
                p = fp(ρ, U0, con, case) 
                u = fu(ρ, U0, con, case,LoR)
               if case.cs == 1   #elastic
                    sxx = Sxx(ρ, U0, con) 
                    EoP =1
                elseif case.cs == 2 || case.cs == 3 # plastic
                    sxx = 2/3*con.Y0
                    EoP = 2
                end

        U[i,1:4] = [ρ, u, p, sxx]
        c = sound(U[i,:], con,EoP)
        if LoR == 1
        uc = u-c 
        else
        uc = u+c
        end
        x[i] = uc*t
 #       println(ρ,"   ",uc₀)
    end
    return x, U
end

end


