# Test based on example 7 of:
# Bürger, Mulet, Villada, A difussion Corrected Multiclass LWR Traffic model
# with Anticipation Lengths and Reaction times, Advances in Applied Mathematics
# and mechanics, 2013

# Parameters:
const CFL = 0.1
const Tend = 0.1
const L = 4.0 #mi
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
  ϕ1 = 0.12
  ϕ2 = 0.4
  δϕ = 0.01
  for (i,x) in enumerate(xx)
    uinit[i,:] = [ϕ1,ϕ2]+δϕ*(cosh(320/L*(x-5*L/16-2.3))^(-2)-0.25*cosh(40/L*(x-11*L/32-2.3))^(-2))
  end
  return dx, xx, uinit
end

#Run Test
include("kt_scheme.jl")
const M = 2
const N = 1000
dx, xx, uinit = setup_initial(N,M)
const Vmax = [80.0,30.0]
const Lmin = [0.01,0.01]
const τ = [0.00095,0.00075]
@time uu =  KT(uinit,dx,CFL,Tend, TVD_RK2, PERIODIC)

#Plot
using(Plots)
plot(xx, uinit, line=(:dot,2))
plot(xx, uu[:,1], line=(:dot,2))
plot!(xx, [sum(uu[i,:]) for i=1:N],lab="ϕ")
