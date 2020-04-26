import sympy as sp
import Functions as fn
sp.init_printing(use_latex='mathjax', pretty_print=True)
D = fn.mass_matrix()
V, U, Jm, r, R, km = fn.V, fn.U, fn.Jm, fn.r, fn.km, fn.R
J, M = sp.zeros(6, 6), sp.zeros(6, 6)
for i in range(6):
    J[i,i] = Jm[i] * r[i]**2
    U[i] = V[i]*r[i]*km[i]/R[i]
M = D + J

    
    
