import sympy as sp
import Functions as fn
sp.init_printing(use_latex='mathjax', pretty_print=True)
q = fn.q
D_H = fn.D_H
T0_6_known = (fn.H_joints_to_base(D_H)[5].subs({q[0]:-sp.pi/3,q[1]:sp.pi/4,q[2]:sp.pi/4,\
              q[3]:sp.pi/3,q[4]:sp.pi/3,q[5]:sp.pi/3})).evalf()#giving position and orientation of ee
Xc = T0_6_known[0,3]#extracting X, Y, Z of ee (Position of ee)
Yc = T0_6_known[1,3]
Zc = T0_6_known[2,3]
R  = T0_6_known[0:3,0:3]#extracting orientation of ee
theta = []
theta[0:3] = fn.first_three_angles(Xc, Yc, Zc)#calculate first three joints angle due to position of ee
R0_3 = fn.H_joints_to_base(D_H)[2][0:3,0:3].subs({q[0]: theta[0], q[1]: theta[1], q[2]: theta[2]})#decoupling
R_3_6 = R0_3.T*R
theta[3:6] = fn.last_three_angles(R_3_6)

print(fn.Matrix_Equality((fn.H_joints_to_base(D_H)[5].subs({q[0]:-sp.pi/3,q[1]:sp.pi/4,q[2]:sp.pi/4,\
              q[3]:sp.pi/3,q[4]:sp.pi/3,q[5]:sp.pi/3})).evalf(),(fn.H_joints_to_base(D_H)[5].\
              subs({q[0]:theta[0],q[1]:theta[1],q[2]:theta[2],\
              q[3]:theta[3],q[4]:theta[4],q[5]:theta[5]})).evalf()))#cheks whether computed angles are True
