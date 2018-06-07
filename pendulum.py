from numpy import sin, cos
from numpy import pi
from numpy import deg2rad
from numpy import arange
import argparse

parser = argparse.ArgumentParser(
    description="Pendulum simulation. Pass length of pendulum and starting angle. Optionally pass method of approximation and mass of the point, and step size.")

parser.add_argument('--length', help="Length of the pendulum", type=float, required=True)
parser.add_argument('--angle', help="Starting angle", type=float, required=True)
parser.add_argument('--method', choices=['euler', 'midpoint', 'mk4'], default='euler',
                    help="Function for approximation")
parser.add_argument('--step', default=0.05, type=float, help="Step size for approximation")
parser.add_argument('--mass', default=1, type=float, help="Mass of the point")
args = parser.parse_args()

length = args.length
alpha = args.angle
func = args.method
dt = args.step
mass = args.mass

print(vars(args))


class Pendulum:
    h0 = 3
    g = 9.81

    def __init__(self, l, alpha, f, dt, m) -> None:
        super().__init__()
        self.f = f
        self.l = l
        self.dt = dt
        self.m = m
        self.alpha = deg2rad(alpha) + pi
        self.omega = 0
        self.x = sin(self.alpha) * self.l
        self.y = cos(self.alpha) * self.l

    def solve_euler(self) -> None:
        self.omega = self.omega + (Pendulum.g / self.l) * sin(self.alpha) * self.dt
        self.alpha = self.alpha + self.omega * self.dt

    def solve_midpoint(self) -> None:
        k1_omega = (Pendulum.g / self.l) * sin(self.alpha)
        k1_alpha = self.omega
        t1_omega = self.omega + k1_omega * self.dt * 0.5
        t1_alpha = self.alpha + k1_alpha * self.dt * 0.5
        k2_omega = (Pendulum.g / self.l) * sin(t1_alpha)
        k2_alpha = t1_omega
        self.omega = self.omega + k2_omega * dt
        self.alpha = self.alpha + k2_alpha * dt

    def solve_mk4(self) -> None:  # TODO
        k1_omega = (Pendulum.g / self.l) * sin(self.alpha)
        k1_alpha = self.omega
        t1_omega = self.omega + k1_omega * self.dt * 0.5
        t1_alpha = self.alpha + k1_alpha * self.dt * 0.5
        k2_omega = (Pendulum.g / self.l) * sin(t1_alpha)
        k2_alpha = t1_omega
        t2_omega = self.omega + k2_omega * self.dt * 0.5
        t2_alpha = self.alpha + k2_alpha * self.dt * 0.5
        k3_omega = (Pendulum.g / self.l) * sin(t2_alpha)
        k3_alpha = t2_omega
        t3_omega = self.omega + k3_omega * self.dt
        t3_alpha = self.alpha + k3_alpha * self.dt
        k4_omega = (Pendulum.g / self.l) * sin(t3_alpha)
        k4_alpha = t3_omega
        self.omega = self.omega + (1 / 6) * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega) * self.dt
        self.alpha = self.alpha + (1 / 6) * (k1_alpha + 2 * k2_alpha + 2 * k3_alpha + k4_alpha) * self.dt

    def do_step(self) -> None:
        self.f(self)
        self.x = sin(self.alpha) * self.l
        self.y = cos(self.alpha) * self.l

    def calc_potential_energy(self) -> float:
        h = self.l + self.y + Pendulum.h0
        return self.m * Pendulum.g * h

    def calc_kinectic_energy(self) -> float:
        return self.m * (self.omega * self.omega * self.l * self.l) / 2


if (func == "euler"):
    func = Pendulum.solve_euler
elif (func == "midpoint"):
    func = Pendulum.solve_midpoint
else:
    func = Pendulum.solve_mk4
pendulum = Pendulum(length, alpha, func, dt, mass)

import matplotlib.pyplot as plt
import matplotlib.animation as animation

t = arange(0.0, 20, 0.05)
fig = plt.figure(figsize=(5, 5))
pendulum_ax = fig.add_subplot(111, autoscale_on=False, xlim=(-pendulum.l - 1, pendulum.l + 1),
                              ylim=(-pendulum.l - 1, pendulum.l + 1))
pendulum_ax.grid()
plt.title('Pendulum %s' % args.method)
line, = pendulum_ax.plot([], [], 'o-', lw=2)
time_template = '%f s'
time_text = pendulum_ax.text(0.05, 0.9, '', transform=pendulum_ax.transAxes)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    pendulum.do_step()
    thisx = [pendulum.x, 0]
    thisy = [pendulum.y, 0]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i * pendulum.dt))

    return line, time_text


anim = []

anim.append(animation.FuncAnimation(fig, animate, frames=1000,
                                    interval=pendulum.dt, blit=True, init_func=init))

# ENERGIA ANIMACJA

const = pendulum.calc_kinectic_energy() + pendulum.calc_potential_energy()

t_energy = arange(0.0, 20, 0.05)
fig_energy = plt.figure(2)
ax_kinectical_energy = fig_energy.add_subplot(111, autoscale_on=True, xlim=(-const - 1, const + 1),
                                              ylim=(-const - 1, const + 1))
ax_kinectical_energy.grid()
plt.title('Pendulum energy')
line_kinectical_energy, = ax_kinectical_energy.plot([], [], 'ro')
energy_txt = ax_kinectical_energy.text(0.05, 0.9, '', transform=ax_kinectical_energy.transAxes)


def init_energy():
    line_kinectical_energy.set_data([], [])
    energy_txt.set_text('')
    return line_kinectical_energy, energy_txt


def animate_energy(i):
    pendulum.do_step()
    thisx = [pendulum.x, pendulum.x]
    potential_energy = pendulum.calc_potential_energy()
    kinectic_energy = pendulum.calc_kinectic_energy()
    thisy = [potential_energy, kinectic_energy]

    line_kinectical_energy.set_data(thisx, thisy)

    sum = (potential_energy + kinectic_energy)
    error = (const - sum)

    energy_txt.set_text('blad = %f, suma = %f' % (error, sum))

    return line_kinectical_energy, energy_txt


anim.append(animation.FuncAnimation(fig_energy, animate_energy, 1000,
                                    interval=pendulum.dt, blit=True, init_func=init_energy))

# # ani.save('double_pendulum.mp4', fps=15)
plt.show()
