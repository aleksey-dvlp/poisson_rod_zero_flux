# poisson_rod_zero_flux

Numerical solution (using FDM - finite difference method) of the equation with the boundaty conditions
u'' + u' = x*exp(x)
u(0) = 5
du(1)/dx = 0

Discretization (where c is a numerical representation of the analitical function u):

n = 1:
c[2] -2*c[1]+c[0]      c[2] - c[0]
-----------------   + ------------- = h*(exp(h)),  c[0] = 5
      h*h                 2*h

n = 2, N-1

c[n+1] -2*c[n]+c[n-1]      c[n+1] - c[n-1]
---------------------  + ------------------ = n*h*(exp(n*h))
      h*h                         2*h

n = N
 2*(c[N-1]-c[N])    
----------------  = N*h*(exp(N*h))
      h*h        
means that the second boundary condition (du(1)/dx=0) approximezed by
adding a fictive node at the coordinate Nh+h (with the enumeration nubmer N+1) and putting there a value
c[N+1] = c[N-1]

The scheme has a second order convergence.
