#from Trajectory_Planning import q, traj_wva, traj_va
#from Inverse_Kinematics import IK
#from DenavHart import T0_6
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import Functions as fn
q = fn.q
D_H = fn.D_H
T0_6 = fn.H_joints_to_base(D_H)[5]
pii = sp.pi
P0 = (T0_6.subs({q[0]:-pii/3, q[1]:pii/4, q[2]:pii/4, q[3]:pii/3, q[4]:pii/3, q[5]:pii/3})).evalf()
Pf = (T0_6.subs({q[0]: pii/3, q[1]:pii/4, q[2]:pii/4, q[3]:pii/3, q[4]:pii/3, q[5]:pii/3})).evalf()

X_line = fn.traj_line(P0, Pf)#x, y, z of trajectory
X_quintic = fn.traj_quintic(P0, Pf)


theta_line = fn.Joint_Space_Traj(X_line[0], X_line[1], X_line[2])
theta_quintic = fn.Joint_Space_Traj(X_quintic[0], X_quintic[3], X_quintic[6])
t = np.linspace(0,10,100)
for i in range(3):
    plt.figure(i)
    plt.plot(t, theta_quintic[i])
    plt.xlabel('time')
    plt.ylabel('theta '+str(i+1))
for i in range(3):
    plt.figure(i)
    plt.plot(t, theta_line[i])
    plt.xlabel('time')
    plt.ylabel('theta '+str(i+1))
    

        
        

