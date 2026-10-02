import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Example data
t = np.linspace(0, 10, 100)
x = np.sin(t)
y = np.cos(t)
z = t

# Create a 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Plot a 3D line
ax.plot(x, y, z, label="3D Line")
ax.scatter(x, y, z, color="red", s=20)  # add points

ax.set_xlabel("X axis")
ax.set_ylabel("Y axis")
ax.set_zlabel("Z axis")
ax.legend()

plt.show()