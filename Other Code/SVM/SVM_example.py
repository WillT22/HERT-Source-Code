import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVC

# Set seed for reproducibility
np.random.seed(42)
n_samples = 100

# --- 1. Data Generation (First Quadrant Only: theta in [0, pi/2]) ---
# Class 0: Inner arc (r in [0, 4])
r0 = np.random.uniform(0, 4, n_samples) + np.random.normal(0, 0.15, n_samples)
r0 = np.clip(r0, 0, None)  # Strictly enforce r >= 0
theta0 = np.random.uniform(0, np.pi / 2, n_samples)  # Quadrant 1

# Class 1: Outer arc (r in [4, 8])
r1 = np.random.uniform(4, 8, n_samples) + np.random.normal(0, 0.15, n_samples)
r1 = np.clip(r1, 0, None)  # Strictly enforce r >= 0
theta1 = np.random.uniform(0, np.pi / 2, n_samples)  # Quadrant 1

# Convert to Cartesian coordinates (x, y)
x0, y0 = r0 * np.cos(theta0), r0 * np.sin(theta0)
x1, y1 = r1 * np.cos(theta1), r1 * np.sin(theta1)

# Enforce strict non-negative Cartesian coordinates
x0, y0 = np.clip(x0, 0, None), np.clip(y0, 0, None)
x1, y1 = np.clip(x1, 0, None), np.clip(y1, 0, None)

# --- 2. Train Linear SVM in Cartesian Space (x, y) ---
X_cartesian = np.vstack([
    np.column_stack((x0, y0)),
    np.column_stack((x1, y1))
])
y = np.array([0] * n_samples + [1] * n_samples)

# Linear kernel fits a straight line: w1*x + w2*y + b = 0
svm = SVC(kernel='linear', C=1)
svm.fit(X_cartesian, y)

# --- 3. Plotting ---
plt.figure(figsize=(7, 7))

# Scatter plot of points with larger marker size (s=80)
plt.scatter(x0, y0, color='red', marker='o', alpha=0.6, s=80, label='Class 0')
plt.scatter(x1, y1, color='blue', marker='^', alpha=0.6, s=80, label='Class 1')

# Grid generation for the first quadrant
x_max = max(x0.max(), x1.max()) + 0.5
y_max = max(y0.max(), y1.max()) + 0.5

xx, yy = np.meshgrid(
    np.linspace(0, x_max, 300),
    np.linspace(0, y_max, 300)
)

# Evaluate decision function on (x, y) grid
Z = svm.decision_function(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)

# Draw the straight line decision boundary (0) and margin lines (-1, 1) with thicker lines
plt.contour(xx, yy, Z, levels=[-1, 0, 1], colors='k', linestyles=['--', '-', '--'], linewidths=[1.5, 2.5, 1.5])

# Scaled text parameters for publication sizing
plt.title("Linear SVM Boundary Example", fontsize=18, fontweight='bold', pad=16)
plt.xlabel("x", fontsize=18)
plt.ylabel("y", fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.legend(loc='upper right', fontsize=16)

plt.xlim(left=0, right=x_max)
plt.ylim(bottom=0, top=y_max)
plt.gca().set_aspect('equal')
plt.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()

plt.savefig("SVM_Linear_Boundary_Example.png", dpi=300, bbox_inches="tight")
plt.savefig("SVM_Linear_Boundary_Example.pdf", bbox_inches="tight")

plt.show()
plt.close()