# Test problem 5
# Shallow water system with flat bottom

# Parameters:
const CFL = 0.1
const Tend = 0.2
#const M = 8000       #Cells in reference solution
const μ = 1e-3
const gr = 9.8    #gravedad

#Functions:
BB(u) = μ*sum(vv(u).^2)*[1.0 0.0; 0.0 1.0]
vv(u) = [gr*u[1]-0.5*(u[2]/u[1])^2 (u[2]/u[1])]
FluxN(ul, ur) = [0.25*(ur[1]+ul[1])*(ur[2]/ur[1]+ul[2]/ul[1]); (0.5*gr*0.25*(ur[1]+ul[1])^2+0.25*0.5*(ur[1]+ul[1])*(ur[2]/ur[1]+ul[2]/ul[1])^2)]
kvisc(ul,ur) = μ*(sum(vv(ul).^2 + vv(ur).^2))/2.0*[1.0 0.0;0.0 1.0]

function JacF(u::Vector)
  h = u[1]
  q = u[2]
  F =[0.0 1.0;-q^2/h^2+gr*h 2*q/h]
end

#Setup initial Conditions
function setup_initial(N)
  dx = 10.0/N
  # (-5, 5)
  xx = [i*dx+dx/2-5.0 for i in (0:N-1)]
  uinit = zeros(N, 2)
  for i = 1:N
    if xx[i] < 0.0
      uinit[i,1] = 2.0
    else
     uinit[i,1] = 1.0
   end
  end
  return dx,xx, uinit
end

include("ec_scheme.jl")
N=500
dx, xx, uinit = setup_initial(N)
@time uu3 = entropy_nonconservative(uinit,dx,CFL,Tend) #ESNC

#Plot
using(Plots)
plot(xx, uinit[:,1], lab="ho",line=(:dot,2))
plot!(xx, uu3[:,1],lab="ESNC h")
plot(xx, uinit[:,2]./uinit[:,1], lab="uo",line=(:dot,2))
plot!(xx, uu3[:,2]./uu3[:,1],lab="ESNC u")
