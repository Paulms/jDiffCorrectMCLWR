# Test based on example 2 of:
# Abeynaike et al, The experimental measurement and modelling of sedimentation
# and creaming for glycerol/biodiesel droplet dispersions, 2012

# 20% glycerol-rich phase

# Parameters:
const CFL = 0.4
const Tend = 300 #s
const L = 20e-3 #m
const ϕc = 0.1 #TODO: why?
const nn = 4.65 # index of Richardson and Zaki
const Do = 1e-7 #m²/s
const ρd = 1090 # kg/m³ dispersed fase (glycerol)
const ρc = 880  # kg/m³ continuous fase (biodiesel)
const gr = 9.81 #m/
const μc = 6.5e-3 # Pa s
const dd = [2, 23, 34, 50, 70, 100, 150, 200]*1e-6 #m
const ut = (ρd - ρc)*gr/18*μc*dd.^2*1e4

#Functions:
Flux(ϕ::Vector) = VV(sum(ϕ))*ϕ.*ut
VV(ϕ) = (1-ϕ)^nn
VP(ϕ) = -nn*(1-ϕ)^(nn-1)

function BB(ϕ::Vector)
  if (sum(ϕ) < ϕc)
    0.0
  else
    M = size(ϕ,1)
    Do*(1-sum(ϕ))^nn*eye(M)
  end
end

function JacF(ϕ::Vector)
  M = size(ϕ,1)
  F = zeros(M,M)
  Vϕ = VV(sum(ϕ))
  VPϕ = VP(sum(ϕ))
  for i =  1:M
    for j = 1:M
      F[i,j]=ut[i]*(((i==j)? Vϕ:0.0) + ϕ[i]*VPϕ)
    end
  end
  F
end

#Function for Entropy C schemes
function FluxN(ul::Vector, ur::Vector)
  VV(sum(ur))*ul.*ut
end

function kvisc(ul::Vector, ur::Vector)
  M = size(ul,1)
  K = zeros(M,M)
  C1 = Do * (1-sum((ul+ur)/2))^nn
  for i = 1:M
      K[i,i] = C1/(ut[i])
  end
  K
end

#Setup initial Conditions
function setup_initial(N,M)
  # We use ghost cells
  dx = L/N
  xx = [i*dx+dx/2 for i in 0:(N-1)]
  uinit = zeros(N,M)
  for (i,x) in enumerate(xx)
      uinit[i,:] = ϕo
  end
  return dx, xx, uinit
end

#Run Test
include("kt_scheme.jl")
include("ec_scheme.jl")
const M = 8
const N = 200
const ϕo = [0, 0.006, 0.018, 0.048, 0.08, 0.042, 0.006, 0]
dx, xx, uinit = setup_initial(N,M)
#@time uu =  KT(uinit,dx,CFL,Tend, TVD_RK2)
@time uu2 =  entropy_nonconservative(uinit,dx,CFL,Tend, TVD_RK2)

#Plot
using(Plots)
plot(xx, uinit, line=(:dot,2))
plot(xx, uu, line=(:dot,2))
plot!(xx, [sum(uu[i,:]) for i=1:N],lab="ϕ")
plot(xx, uu2, line=(:dot,2))
plot!(xx, [sum(uu2[i,:]) for i=1:N],lab="ϕ")
