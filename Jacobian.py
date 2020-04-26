import sympy as sp
import Functions as fn
q = fn.q
sp.init_printing(use_latex='mathjax', pretty_print=True)
D_H = fn.D_H
Jacobian_of_end_effector = fn.geo_jacobian(6)
jacobian_of_4th_joint = fn.geo_jacobian(4)
analytic_jacobian_of_end_effector = fn.anal_jacobian(6)
determinant_of_J = sp.trigsimp(sp.det(Jacobian_of_end_effector[0:3,0:3])*\
                               sp.det(Jacobian_of_end_effector[3:6,3:6]))

                        
                    
             
