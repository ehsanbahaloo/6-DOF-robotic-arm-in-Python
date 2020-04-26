from sympy import init_printing, pi
import Functions as fn
init_printing(use_latex='mathjax', pretty_print=True)
alpha = fn.alpha
a = fn.a
d = fn.d
q = fn.q
D_H = { alpha[0]:      0, a[0]:      0,   d[0]:  0.75,   q[0]:         q[0],
        alpha[1]:  -pi/2, a[1]:   0.35,   d[1]:     0,   q[1]:   -pi/2+q[1],
        alpha[2]:      0, a[2]:   1.25,   d[2]:     0,   q[2]:         q[2],
        alpha[3]:  -pi/2, a[3]: -0.055,   d[3]:   1.1,   q[3]:         q[3],
        alpha[4]:   pi/2, a[4]:      0,   d[4]:     0,   q[4]:         q[4],
        alpha[5]:  -pi/2, a[5]:      0,   d[5]:     0,   q[5]:         q[5]}

T0_3 = fn.H_joints_to_base(D_H)[2]
T0_6 = fn.H_joints_to_base(D_H)[5]
T3_6 = fn.H_joints_to_base(D_H)[6]

