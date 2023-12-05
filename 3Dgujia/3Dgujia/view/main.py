import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

fig=plt.figure()
ax = fig.add_subplot(projection="3d")
# cnt = 0
All = np.loadtxt('../stl.txt', delimiter=' ')
All = All.T
print(All[0].size)

#ax.scatter(All[0][::100], All[1][::100], All[2][::100], c='r')

All = np.loadtxt('../out.txt', delimiter=' ')
All = All.T
print(All[0].size)
ax.scatter(All[0][::], All[1][::], All[2][::], c='b')

All = np.loadtxt('../cp.txt', delimiter=' ')
All = All.T
print(All[0].size)
#ax.scatter(All[0][::], All[1][::], All[2][::], c='b')

All = np.loadtxt('../qiepian.txt', delimiter=' ')
All = All.T
print(All[0].size)
ax.scatter(All[0][::], All[1][::], All[2][::], c='g')

# x = []
# y = []
# z = []
# cnt = 0
# while True:
#     r = input()
#     if r == 'e':
#         break
#     R = list(map(float, r.split()))
#     x.append(R[0])
#     y.append(R[1])
#     z.append(R[2])
#
# xx = np.array(x)
# yy = np.array(y)
# zz = np.array(z)
#
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.scatter(xx, yy, zz, c='b')
#
plt.show()
input()
