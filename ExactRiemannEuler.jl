

module Riemann

export W, SolveDistr

struct W
    ρ::Float64
    u::Float64
    p::Float64
end

function f(p::Float64,W₀::W,γ::Float64)
    ρ₀,u₀,p₀ = W₀.ρ,W₀.u,W₀.p
    
    A = 2/((γ +1 )*ρ₀)
    B = (γ-1)/(γ+1)*p₀
    
    if p>p₀ 
        return (p-p₀)*√(A/(p+B))
    else
        a₀ = √(γ*p₀/ρ₀)
        return 2a₀/(γ-1)*((p/p₀)^((γ-1)/2/γ) -1)
    end
end    

function ρₚ(p::Float64,W₀::W,γ::Float64)
    if p < W₀.p
       return W₀.ρ*(p/W₀.p)^(1/γ)
    else 
        return W₀.ρ*((γ-1)/(γ+1)+p/W₀.p)/((γ-1)/(γ+1)*p/W₀.p+1)
    end
end

function Derivative(f::Function,x::Float64)
    ϵ₀ = 1e-10
    if abs(x) >= ϵ₀
        ϵ= ϵ₀*x
    else
        ϵ =ϵ₀
    end
    f¹(x) = (f(x+ϵ)-f(x))/ϵ
    return f¹(x)
end

    fₗ(p) = f(p,Wₗ,γ)
    fᵣ(p) = f(p,Wᵣ,γ)
    F(p) = fₗ(p) + fᵣ(p)+Wᵣ.u-Wₗ.u
    F¹(p) = Derivative(p->F(p),p)

function RiemannSolver(Wₗ::W,Wᵣ::W,γ::Float64)    
   #     γ = 1.4
     p = (Wₗ.p + Wᵣ.p)/2
     TOL = 1.0e-4
   #  i = 1
    
    ϵ = 1.0e-8
    
    WₗStar = W(ϵ,0,ϵ)
    WᵣStar = W(ϵ,0,ϵ)
    
     aₗ = √(γ*Wₗ.p/Wₗ.ρ) 
     aᵣ = √(γ*Wᵣ.p/Wᵣ.ρ) 
     cri = aₗ/(γ-1) + aᵣ/(γ-1) +Wₗ.u-Wᵣ.u
    
    if cri >0
        while true 
            F₁ = F(p)
             p1 = p- F₁/F¹(p)
             if max(abs(F₁),2*abs(p1-p)/(p1+p)) < TOL
                break
            end 
        #    i += 1
            p=max(p1,TOL)
         #   @show p,F₁  
        end

        pStar = p
        uStar = Wₗ.u-fₗ(p)
        ρₗStar = ρₚ(p,Wₗ,γ)
        ρᵣStar = ρₚ(p,Wᵣ,γ)
        WₗStar = W(ρₗStar,uStar,pStar)
        WᵣStar = W(ρᵣStar,uStar,pStar)
    end
    
    return WₗStar,WᵣStar
end

function SolveDistr(t::Float64,Wₗ::W,Wᵣ::W,γ::Float64)
    
    WₗStar,WᵣStar = RiemannSolver(Wₗ,Wᵣ,γ)
    
    ρₗStar,uStar, pStar =WₗStar.ρ,WₗStar.u,WₗStar.p
    ρᵣStar = WᵣStar.ρ
    ρₗ,pₗ,uₗ = Wₗ.ρ,Wₗ.p,Wₗ.u
    ρᵣ,pᵣ,uᵣ = Wᵣ.ρ,Wᵣ.p,Wᵣ.u
    
    I = 1200
    U = zeros(Float64,(I,3))
    x = zeros(Float64,I)
    if pStar > pₗ  # Left shock
        sL = (ρₗStar*uStar - ρₗ*uₗ)/(ρₗStar-ρₗ)
        for i =1:200
            x[i] = (2-0.005*(i-1))*sL*t
            U[i,:] = [ρₗ,uₗ,pₗ]
        end
        for i = 201:600
            x[i] = ((uStar-sL)/400*(i-200) + sL)*t
            U[i,:] = [ρₗStar,uStar,pStar]
        end
    else  # Left Rarefaction
        cL = √(γ*pₗ/ρₗ)
        cLStar = √(γ*pStar/ρₗStar)
        sLSlow = uₗ - cL
        sLFast = uStar - cLStar
        for i = 1:200
            x[i] =(2-0.005*(i-1))*sLSlow*t
            U[i,:] = [ρₗ,uₗ,pₗ]
        end
        
        x[201:400],U[201:400,:] = RareSpeedToU(t,Wₗ,pStar,γ,1)
        for i =401:600
            x[i] = ((uStar - sLFast)/200*(i-400) + sLFast)*t 
            U[i,:] = [ρₗStar,uStar,pStar]
        end 
    end
    
    if pStar > pᵣ  # Right shock
        sR = (ρᵣStar*uStar - ρᵣ*uᵣ)/(ρᵣStar-ρᵣ)
      #  @show ρᵣStar, ρᵣ
        for i =601:1000
            x[i] = ((-uStar+sR)/400*(i-600) + uStar)*t
            U[i,:] = [ρᵣStar,uStar,pStar]
        end
        for i = 1001:1200
            x[i] = (1+0.005*(i-1001))*sR*t
            U[i,:] = [ρᵣ,uᵣ,pᵣ]
        end
    else  # Right Rarefaction
        cR = √(γ*pᵣ/ρᵣ)
        cRStar = √(γ*pStar/ρᵣStar)
        sRFast = uᵣ + cR
        sRSlow = uStar + cRStar
        
        for i = 1000:1200
           x[i] =(1+0.005*(i-1001))*sRFast*t
            U[i,:] = [ρᵣ,uᵣ,pᵣ]
        end
        
        x[801:1000],U[801:1000,:] = RareSpeedToU(t,Wᵣ,pStar,γ,2)
        for i =601:800
             x[i] = ((-uStar + sRSlow)/200*(i-600) + uStar)*t 
             U[i,:] = [ρᵣStar,uStar,pStar]
        end 
    end
    return x, U
end

function RareSpeedToU(t::Float64, W₀::W, p::Float64,γ::Float64,LoR::Int)
    ρ₀,u₀,p₀ = W₀.ρ, W₀.u,W₀.p
    
    x = zeros(Float64,200)
    U = zeros(Float64,(200,3))
    
    Δp = (p-p₀)/200
    
    for i in 1:200
        p = p₀ + Δp*(i-1)
        ρ = ρₚ(p, W₀, γ)
        c= √(γ*p/ρ)
        f₀ = f(p, W₀, γ)
        if LoR == 1
            u = u₀-f₀
            uc = u-c
        x[i] = uc*t
        U[i,1:3] = [ρ,u,p]
        else
            u = u₀+f₀
            uc = u+c
        x[200-i+1] = uc*t
        U[200-i+1,1:3] = [ρ,u,p]
        end
        
    end
    return x, U            
end

Wₗ = W(1.0,0.0,1.0)
Wᵣ = W(0.125,0.0,0.1)
γ = 1.4
x,U = SolveDistr(0.1,Wₗ,Wᵣ,γ)

end


