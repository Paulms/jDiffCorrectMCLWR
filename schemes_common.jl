@enum StepMethod FORWARD_EULER TVD_RK2 RK4
@enum BoundaryCondition ZERO_FLUX PERIODIC
#using ForwardDiff

#Flux must be like flux(x::Vector) = x.^2
#Jf = x -> ForwardDiff.jacobian(Flux,x)
@inline function fluxρ(uj::Vector)
  #maximum(abs(eigvals(Jf(uj))))
  maximum(abs(eigvals(JacF(uj))))
end

function cdt(u::AbstractArray, CFL, dx)
  maxρ = 0
  maxρB = 0
  N = size(u,1)
  for i in 1:N
    maxρ = max(maxρ, fluxρ(u[i,:]))
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

#Macros
macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

@def boundary_header begin
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
end

@def boundary_update begin
  hhleft = 0; hhright = 0; ppleft = 0; ppright = 0
  if boundary == PERIODIC
    hhleft = hh[1,:]; ppleft = pp[1,:]
    hhright = hh[N-1,:]; ppright = pp[N-1,:]
  end
end

@def update_rhs begin
  j = 1 + ss
  rhs[j-ss,:] = - 1/dx * (hh[j,:] -hhleft - (pp[j,:]-ppleft))
  for j = (2+ss):(N-1-ss)
    rhs[j-ss,:] = - 1/dx * (hh[j,:]-hh[j-1,:]-(pp[j,:]-pp[j-1,:]))
  end
  j = N-ss
  rhs[j-ss,:] =  -1/dx*(hhright-hh[j-1,:]-(ppright - pp[j-1,:]))
end
