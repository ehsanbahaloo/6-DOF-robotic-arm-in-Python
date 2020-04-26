import sympy as sp
import Functions as fn
sp.init_printing(use_latex='mathjax', pretty_print=True)
qDot = sp.Matrix(sp.symbols('qDot1:7'))
xDot = sp.Matrix(sp.symbols('xDot1:7'))
qDot = fn.geo_jacobian(6).inv()*xDot