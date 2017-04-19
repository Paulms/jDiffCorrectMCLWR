using ForwardDiff

@enum StepMethod FORWARD_EULER TVD_RK2 RK4
@enum BoundaryCondition ZERO_FLUX PERIODIC

function KT(uinit::AbstractArray,dx,CFL,Tend,tempSteps = FORWARD_EULER, boundary = ZERO_FLUX, Θ = 1, order=2)
  uu = copy(uinit)
  #Print progress
  percentage = 0
  limit = Tend/5
  println("Starting KT scheme")
  t = 0.0
  M = size(uinit,2)
  N = size(uinit,1)
  rhs1 = zeros(N,M)
  rhs2 = zeros(N,M)
  rhs3 = zeros(N,M)
  rhs4 = zeros(N,M)
  while t <= Tend
    uold = copy(uu)
    dt = cdt(uold, CFL, dx)
    if (tempSteps == FORWARD_EULER)
      update_KT(rhs1, uold, N, M, dx, Θ, boundary, order)
      uu = uold + dt*rhs1
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_KT(rhs1, uold, N, M, dx, Θ, boundary, order)
      #Second Step
      update_KT(rhs2, uold+dt*rhs1, N, M,dx, Θ, boundary, order)
      uu = 0.5*(uold + uold + dt*rhs1 + dt*rhs2)
    elseif (tempSteps == RK4)
      #FIRST Step
      update_KT(rhs1, uold, N, M,dx, boundary)
      #Second Step
      update_KT(rhs2, uold+dt/2*rhs1, N, M,dx, Θ, boundary, order)
      #Third Step
      update_KT(rhs3, uold+dt/2*rhs2, N, M,dx, Θ, boundary, order)
      #Fourth Step
      update_KT(rhs4, uold+dt*rhs3, N, M,dx, dt, Θ, boundary, order)
      uu = uold + dt/6 * (rhs1+2*rhs2+2*rhs3+rhs4)
    end
    # Print Progress
    if (t > limit)
      percentage = percentage + 20
      limit = limit + Tend/5
      println(percentage, "% completed")
    end
    t = t+dt
  end
  println("Completed...")
  return uu
end

function update_KT(rhs, uold, N, M, dx, dt,Θ, boundary, order)
  if (order == 1)
    update_KT1(rhs, uold, N, M, dx, boundary)
  else
    update_KT2(rhs, uold, N, M, dx, dt,Θ, boundary)
  end
end

#Flux must be like flux(x::Vector) = x.^2
Jf = x -> ForwardDiff.jacobian(Flux,x)
function fluxρ(uj::Vector)
  abs(eigvals(Jf(uj)))
end

function cdt(u::AbstractArray, CFL, dx)
  maxρ = 0
  maxρB = 0
  N = size(u,1)
  for i in 1:N
    maxρ = max(maxρ, maximum(fluxρ(u[i,:])))
    maxρB = max(maxρB, maximum(abs(eigvals(BB(u[i,:])))))
  end
  CFL/(1/dx*maxρ+1/(2*dx^2)*maxρB)
end

function minmod(a,b,c)
  if (a > 0 && b > 0 && c > 0)
    min(a,b,c)
  elseif (a < 0 && b < 0 && c < 0)
    max(a,b,c)
  else
    zero(a)
  end
end
function minmod(a,b)
  0.5*(sign(a)+sign(b))*min(abs(a),abs(b))
end

function update_KT1(rhs, uold, N, M, dx, boundary)
  # Impose boundary condition
  ss = 0
  if boundary == PERIODIC
    ss = 1
    N = N + 2   #Create ghost cells
    utemp = copy(uold)
    uold = zeros(N,M)
    uold[2:N-1,:] = utemp
    uold[1,:] = utemp[N-2,:]
    uold[N,:] = utemp[1,:]
  end
  # Slope vector
  ∇u = zeros(N,M)
  for j = 2:N-1
    ∇u[j,:] = minmod.((uold[j,:] - uold[j-1,:])/dx,(uold[j+1,:] - uold[j,:])/dx)
  end
  uplus = uold[2:N,:] - dx/2*∇u[2:N,:]
  uminus = uold[1:N-1,:] + dx/2*∇u[1:N-1,:]
  for j = 1:(N-1)
    aa[j]=max(fluxρ(uminus[j,:]),fluxρ(uplus[j,:]))
  end

  # Numerical Fluxes
  hh = zeros(N-1,M)
  for j = 1:N-1
    hh[j,:] = 0.5*(Flux(uplus[j,:])+Flux(uminus[j,:]))-0.5*(aa[j]*(uplus[j]-uminus[j]))
  end
  ∇u_ap = (uold[2:N,:]-uold[1:N-1,:])/dx
  # Diffusion
  pp = zeros(N-1,M)
  for j = 1:N-1
    pp[j,:] = 0.5*(BB(uold[j+1,:])+BB(uold[j,:]))*∇u_ap[j,:]'
  end
  hhleft = 0; hhright = 0; ppleft = 0; ppright = 0
  if boundary == PERIODIC
    hhleft = hh[1]; ppleft = pp[1]
    hhright = hh[N]; ppright = pp[N]
  end

  j = 1 + ss
  rhs[j-ss] = - 1/dx * (hh[j] -hhleft - (pp[j]-ppleft))
  for j = (2+ss):(N-1-ss)
    rhs[j-ss] = - 1/dx * (hh[j]-hh[j-1]-(pp[j]-pp[j-1]))
  end
  j = N-ss
  rhs[j-ss] =  -1/dx*(hhright-hh[j-1]-(ppright - pp[j-1]))
end


function update_KT2(rhs, uold, N, M, dx, dt,Θ, boundary)
  ss = 0
  if boundary == PERIODIC
    ss = 1
    N = N + 2   #Create ghost cells
    utemp = copy(uold)
    uold = zeros(N,M)
    uold[2:N-1,:] = utemp
    uold[1,:] = utemp[N-2,:]
    uold[N,:] = utemp[1,:]
  end

  #Compute diffusion
  λ = dt/dx
  #update vector
  # 1. slopes
  ∇u = zeros(N,M)
  for i = 1:M
    for j = 2:(N-1)
      ∇u[j,i] = minmod(Θ*(uold[j,i]-uold[j-1,i]),(uold[j+1,i]-uold[j-1,i])/2,Θ*(uold[j+1,i]-uold[j,i]))
    end
  end
  # Local speeds of propagation
  uminus = uold[1:N-1,:]+0.5*∇u[1:N-1,:]
  uplus = uold[2:N,:]-0.5*∇u[2:N,:]
  for j = 1:(N-1)
    aa[j]=max(fluxρ(uminus[j,:]),fluxρ(uplus[j,:]))
  end

  #Flux slopes
  u_l = zeros(N-1,M)
  u_r = zeros(N-1,M)
  for i = 1:M
    for j = 2:N
      u_l[j-1,i] = uold[j-1,i] + (0.5-λ*aa[j-1])*∇u[j-1,i]
      u_r[j-1,i] = uold[j,i] - (0.5-λ*aa[j-1])*∇u[j,i]
    end
  end
  ∇f_l = zeros(N-1,M)
  ∇f_r = zeros(N-1,M)

  for j = 2:(N-2)
    ∇f_l[j,:] = minmod.(Θ*(Flux(u_l[j,:])-Flux(u_l[j-1,:])),(Flux(u_l[j+1,:])-Flux(u_l[j-1,:]))/2,
    Θ*(Flux(u_l[j+1,:])-Flux(u_l[j,:])))
    ∇f_r[j,:] = minmod.(Θ*(Flux(u_r[j,:])-Flux(u_r[j-1,:])),(Flux(u_r[j+1,:])-Flux(u_r[j-1,:]))/2,
    Θ*(Flux(u_r[j+1,:])-Flux(u_r[j,:])))
  end

  # Predictor solution values
  Φ_l = u_l - λ/2*∇f_l
  Φ_r = u_r - λ/2*∇f_r

  # Aproximate cell averages
  Ψr = zeros(N-1,M)
  Ψ = zeros(N,M)
  for j = 1:(N-1)
    if (aa[j] != 0)
      Ψr[j,:] = 0.5*(uold[j,:]+uold[j+1,:])+(1-λ*aa[j])/4*(∇u[j,:]-∇u[j+1,:])-1/(2*aa[j])*
      (Flux(Φ_r[j,:])-Flux(Φ_l[j,:]))
    else
      Ψr[j,:] = 0.5*(uold[j,:]+uold[j+1,:])
    end
  end
  for j = 2:(N-1)
    Ψ[j,:] = uold[j,:] - λ/2*(aa[j]-aa[j-1])*∇u[j,:]-λ/(1-λ*(aa[j]+aa[j-1]))*
    (Flux(Φ_l[j,:])-Flux(Φ_r[j-1,:]))
  end

  # Discrete derivatives
  ∇Ψ = zeros(N-1,M)
  for j = 2:(N-2)
    ∇Ψ[j,:]=2\dx*minmod.(Θ*(Φ_r[j,:]-Ψ[j,:])/(1+λ*(aa[j]-aa[j-1])),
    (Ψ[j+1,:]-Ψ[j,:])/(2+λ*(2*aa[j]-aa[j-1]-aa[j+1])),
    Θ*(Ψ[j+1,:]-Φ_r[j,:])/(1+λ*(aa[j]-aa[j+1])))
  end

  # Numerical Fluxes
  hh = zeros(N-1,M)
  for j = 1:(N-1)
    hh[j] = 0.5*(Flux(Φ_r[j,:])+Flux(Φ_l[j,:]))-0.5*(uold[j+1,:]-uold[j,:])*aa[j]+
    aa[j]*(1-λ*aa[j])/4*(∇u[j+1,:]+∇u[j,:]) + λ*dx/2*(aa[j])^2*∇Ψ[j,:]
  end
  ∇u_ap = (uold[2:N,:]-uold[1:N-1,:])/dx
  # Diffusion
  pp = zeros(N-1,M)
  for j = 1:N-1
    pp[j,:] = 0.5*(BB(uold[j+1,:])+BB(uold[j,:]))*∇u_ap[j,:]'
  end
  hhleft = 0; hhright = 0; ppleft = 0; ppright = 0
  if boundary == PERIODIC
    hhleft = hh[1]; ppleft = pp[1]
    hhright = hh[N]; ppright = pp[N]
  end
  j = 1 + ss
  rhs[j-ss] = - 1/dx * (hh[j] -hhleft - (pp[j]-ppleft))
  for j = (2+ss):(N-1-ss)
    rhs[j-ss] = - 1/dx * (hh[j]-hh[j-1]-(pp[j]-pp[j-1]))
  end
  j = N-ss
  rhs[j-ss] =  -1/dx*(hhright-hh[j-1]-(ppright - pp[j-1]))
end
