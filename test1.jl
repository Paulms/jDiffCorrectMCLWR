# Parameters:
const CFL = 0.25
const Tend = 0.2
const L = 10.0 #mi
const ϕc = exp(-7/e)


#Functions:
Flux(ϕ::Vector) = VV(sum(ϕ))*ϕ.*Vmax
function BB(ϕ::Vector)
  if (sum(ϕ) < ϕc)
    0.0
  else
    M = size(ϕ,1)
    B = zeros(M,M)
    for i = 1:M
     B[i,:] = -VP(sum(ϕ))*(Lmin[i]+τ[i]*(VP(sum(ϕ))*dot(ϕ,Vmax)+(Vmax-Vmax[i])*VV(sum(ϕ))))*ϕ[i]*Vmax[i]
    end
    B
  end
end
VV(ϕ) = (sum(ϕ) < ϕc) ? 1.0 : -e/7*log(ϕ)
VP(ϕ) = (sum(ϕ) < ϕc) ? 0.0 : -e/7*1/ϕ
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
const N = 100
dx, xx, uinit = setup_initial(N,M)
const Vmax = [60,55,50,45]
const Lmin = [0.03,0.03,0.03,0.03]
const τ = [0.0013,0.0011,0.0008,0.0006]
@time uu =  KT(uinit,dx,CFL,N,M,Tend, TVD_RK2, PERIODIC)

#Plot
using(Plots)
plot(xx, uinit, line=(:dot,2))
