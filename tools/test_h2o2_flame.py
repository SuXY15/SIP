import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
import scipy.integrate

gas = ct.Solution("h2o2.yaml")

# 0D auto-ignition
def sim_0D(Tin = 1000):
    gas.TPX = Tin, ct.one_atm, "H2:2.0, O2:1.0"

    r = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([r])
    states = ct.SolutionArray(gas, extra=['t'])
    t = 0
    while t < 1e-2:
        t = sim.step()
        states.append(r.thermo.state, t=t)
    idt = states.t[np.diff(states.T).argmax()]
    plt.plot(states.t, states.T)
    plt.xlim([0, 2*idt])

    print("sim_0D, thermal_conductivity =", gas.thermal_conductivity)

# 0D with global reaction
def global_0D(Tin = 1000):
    gamma = 1.2
    Ru = 692.25
    T0 = 300
    Ea = 27 * Ru * T0
    Q = 36 * Ru * T0
    A = 2e7
    class ODE:
        def __init__(self, gas):
            self.P = gas.P

        def __call__(self, t, y):
            T = y[0]
            Y = y[1]
            rho = self.P / T / Ru
            omega = A * Y * np.exp(-Ea / Ru / T)
            dY = - rho * omega
            dT = - Q / Ru * dY * (gamma-1)
            return [dT, dY]
    solver = scipy.integrate.ode(ODE(gas))
    solver.set_integrator('vode', method='bdf')
    solver.set_initial_value([Tin, 1.0], 0)
    t_end = 1e-2
    states = [[0, Tin, 1.0]]
    while solver.t < t_end:
        solver.integrate(solver.t + 1e-6)
        states.append([solver.t, solver.y[0], solver.y[1]])
    states = np.array(states)
    plt.plot(states[:,0], states[:,1], '--')


# 1D flame
def sim_1D():
    gas.TPX = 300, ct.one_atm, "H2:2.0, O2:1.0"

    f = ct.FreeFlame(gas, width=0.03)
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    f.transport_model = "Multi"
    f.solve(loglevel=0, auto=True)

    print("S_L =", f.velocity[0])
    print("rho0 =", f.density[0])
    print("Tmax =", f.T[-1])

    fig, axs = plt.subplots(3, 1, figsize=(4,6))
    fig.subplots_adjust(top=0.95, right=0.95, left=0.2, bottom=0.1, wspace=0.35)
    axs[0].plot(f.grid, f.T)
    axs[0].set_ylabel("Temperature")
    axs[1].plot(f.grid, f.density)
    axs[1].set_ylabel("Density")
    axs[2].plot(f.grid, f.velocity)
    axs[2].set_ylabel("velocity")

    fig, axs = plt.subplots(3, 1, figsize=(4,6))
    fig.subplots_adjust(top=0.95, right=0.95, left=0.2, bottom=0.1, wspace=0.35)
    axs[0].plot(f.grid, f.cp / f.cv)
    axs[0].set_ylabel("gamma")
    species_diff = [np.dot(f.Y[:,i], f.mix_diff_coeffs[:,i]) for i in range(len(f.grid))]
    axs[1].plot(f.grid, f.thermal_conductivity)
    axs[1].set_ylabel("f.thermal_conductivity")
    axs[2].plot(f.grid, f.thermal_conductivity / f.density / f.cp, label="thermal alpha")
    axs[2].plot(f.grid, species_diff, label="species D")
    axs[2].set_ylabel("diffusivity")
    axs[2].legend()


sim_0D(1000)
global_0D(1000)

# sim_0D(990)
# global_0D(990)

sim_1D()

plt.show()