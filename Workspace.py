import sympy as sp
import matplotlib.pyplot as plt
import Functions as fn
Q_range = sp.Matrix([[-185, -45, -210, -350, -125, -350],
                     [ 185,  85,   65,  350,  125,  350]])
q = fn.q
w = fn.reach_wspace(Q_range)
x, y, z = [], [], []
for i in range(len(w)):
    x.append(w[i][0,0])
    y.append(w[i][1,0])
    z.append(w[i][2,0])
fig = plt.figure()    
plt.scatter(x, z)
plt.show()