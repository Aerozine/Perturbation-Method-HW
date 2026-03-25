import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
#plt.style.use("seaborn-v0_8")  
plt.rcParams.update({
    "font.family": "serif",
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

colours = {
    "numerical": "blue",
    "composite": "orange",
    "outer":     "green",
    "inner":     "red",
}
def solve_numerical(eps, n_init=300):
    # eps*y'' = 2(x^2+1)*y' - 1 , y(0)=y(1)=0
    def ODE(x, y):
        return np.vstack((y[1], (2*(x**2 + 1)*y[1] - 1) / eps))
    def bc(ya, yb):
        return np.array([ya[0], yb[0]])
    x0 = np.linspace(0, 1, n_init)
    y0 = np.zeros((2, x0.size))
    sol = solve_bvp(ODE, bc, x0, y0, tol=1e-10)
    return sol

def y_outer(x, eps):
    y0 = 0.5 * np.arctan(x)
    y1 = 1 / (8*(1 + x**2)**2) - 1/8
    return y0 + eps*y1

def y_inner(x, eps):
    xx  = (1 - x) / eps
    e   = np.exp(-4 * xx)
    y0  = -np.pi/8 * (e - 1)
    y1  = (3/32 + (-np.pi/8)*xx + 2*(-np.pi/8)*xx**2)*e - 3/32 - xx/4
    return y0 + eps*y1

def y_match(x,eps):
    return np.pi/8 + (x-1)/4 - 3*eps/32

def y_composite(x,eps):
    return y_outer(x,eps)+y_inner(x,eps)-y_match(x,eps)

if __name__ == "__main__":
    eps_arr = [0.2,0.01]
    for eps in eps_arr:
        # compute
        sol = solve_numerical(eps)
        x_plot = np.linspace(0, 1, 1000)
        y_num = sol.sol(x_plot)[0]
        y_out = y_outer(x_plot, eps)
        y_in  = y_inner(x_plot, eps)
        y_comp = y_composite(x_plot, eps)
        #plot
        plt.figure(figsize=(10, 6))
        plt.plot(x_plot, y_num,
                label="Numerical (BVP)",
                color=colours["numerical"],
                linestyle="--",
                linewidth=2,
                 alpha=0.7
                 )
        plt.plot(x_plot, y_comp,
                label="Composite",
                color=colours["composite"],
                linewidth=2,
                 alpha=0.7
                 )
        plt.plot(x_plot, y_out,
                label="Outer",
                color=colours["outer"],
                linestyle=":",
                 alpha=0.7)
        plt.plot(x_plot, y_in,
                label="Inner",
                color=colours["inner"],
                linestyle="-.",alpha=0.7)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"Comparison of Numerical and Asymptotic Solutions (ε = {eps})")
        plt.xlim(0, 1)
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"figures/plot_{eps}.svg")
        plt.show()
    
