# Test based on example 1 of:
# Bürger, Mulet, Villada, A difussion Corrected Multiclass LWR Traffic model
# with Anticipation Lengths and Reaction times, Advances in Applied Mathematics
# and mechanics, 2013

# Parameters:
const CFL = 0.25
const Tend = 0.2
const L = 10.0 #mi
const ϕc = exp(-7/e)

#Functions:
Flux(ϕ::Vector) = VV(sum(ϕ))*ϕ.*Vmax
function JacF(ϕ::Vector)
  M = size(ϕ,1)
  F = zeros(M,M)
  Vϕ = VV(sum(ϕ))
  VPϕ = VP(sum(ϕ))
  for i =  1:M
    for j = 1:M
      F[i,j]=Vmax[i]*(((i==j)? Vϕ:0.0) + ϕ[i]*VPϕ)
    end
  end
  F
end

function BB(ϕ::Vector)
  if (sum(ϕ) < ϕc)
    0.0
  else
    M = size(ϕ,1)
    B = zeros(M,M)
    Vϕ = VV(sum(ϕ))
    VPϕ = VP(sum(ϕ))
    ϕVm = VPϕ*dot(ϕ,Vmax)
    for i = 1:M
     B[i,:] = -VPϕ*(Lmin[i]+τ[i]*(ϕVm+(Vmax-Vmax[i])*Vϕ))*ϕ[i]*Vmax[i]
    end
    B
  end
end
VV(ϕ) = (ϕ < ϕc) ? 1.0 : -e/7*log(ϕ)
VP(ϕ) = (ϕ < ϕc) ? 0.0 : -e/7*1/ϕ
#Setup initial Conditions
function setup_initial(N,M)
  # We use ghost cells
  dx = L/N
  xx = [i*dx+dx/2 for i in 0:(N-1)]
  uinit = zeros(N,M)
  for (i,x) in enumerate(xx)
    if (0.0<x<=0.1)
      uinit[i,:] = 10*x*[0.2,0.3,0.2,0.3]
    elseif 0.1<x<=0.9
      uinit[i,:] = [0.2,0.3,0.2,0.3]
    elseif 0.9<x<=1
      uinit[i,:] = -10*(x-1)*[0.2,0.3,0.2,0.3]
    else
      uinit[i,:] = 0.0
    end
  end
  return dx, xx, uinit
end

#Run Test
include("kt_scheme.jl")
const M = 4
const N = 500
dx, xx, uinit = setup_initial(N,M)
const Vmax = [60.0,55.0,50.0,45.0]
const Lmin = [0.03,0.03,0.03,0.03]
const τ = [0.0013,0.0011,0.0008,0.0006]
@time uu =  KT(uinit,dx,CFL,Tend, TVD_RK2, PERIODIC)
#@time uu2 = KT(uinit,dx,CFL,Tend, TVD_RK2, PERIODIC, 1,1)

#Plot
using(Plots)
plot(xx, uinit, line=(:dot,2))
plot(xx, uu, line=(:dot,2))
plot!(xx, [sum(uu[i,:]) for i=1:N],lab="ϕ")
#plot(xx, uu2, line=(:dot,2))
#plot!(xx, [sum(uu2[i,:]) for i=1:N],lab="ϕ")
