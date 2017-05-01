include("schemes_common.jl")

function entropy_nonconservative(uinit::AbstractArray,dx,CFL,Tend; tempSteps = FORWARD_EULER,
  boundary = ZERO_FLUX,  ϵ = 0.0, Extra_Viscosity = false)
  uu = copy(uinit)
  #Print progress
  percentage = 0
  limit = Tend/5
  println("Starting Entropy Stable Non conservative scheme...")
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
      update_enc(rhs1, uold, N, M, dx, dt, boundary, ϵ, Extra_Viscosity)
      uu = uold + dt*rhs1
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_enc(rhs1, uold, N, M, dx, dt, boundary, ϵ, Extra_Viscosity)
      #Second Step
      update_enc(rhs2, uold+dt*rhs1, N, M,dx, dt, boundary, ϵ, Extra_Viscosity)
      uu = 0.5*(uold + uold + dt*rhs1 + dt*rhs2)
    elseif (tempSteps == RK4)
      #FIRST Step
      update_enc(rhs1, uold, N, M,dx, dt,  boundary, ϵ, Extra_Viscosity)
      #Second Step
      update_enc(rhs2, uold+dt/2*rhs1, N, M,dx, dt, boundary, ϵ, Extra_Viscosity)
      #Third Step
      update_enc(rhs3, uold+dt/2*rhs2, N, M,dx, dt, boundary, ϵ, Extra_Viscosity)
      #Fourth Step
      update_enc(rhs4, uold+dt*rhs3, N, M,dx, dt,  boundary, ϵ, Extra_Viscosity)
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

function update_enc(rhs, uold, N, M,dx, dt,  boundary, ϵ, Extra_Viscosity)
  @boundary_header
  # Numerical Fluxes
  hh = zeros(N-1,M)
  for j = 1:(N-1)
    hh[j,:] = FluxN(uold[j,:], uold[j+1,:]) + ϵ/dx*(Extra_Viscosity ? uold[j+1,:]-uold[j,:] : 0.0)
  end
  # Diffusion
  pp = zeros(N-1,M)
  for j = 1:N-1
    pp[j,:] = 1/dx*(kvisc(uold[j,:],uold[j+1,:])*(uold[j+1,:]-uold[j,:]))
  end
  @boundary_update
  @update_rhs
end
