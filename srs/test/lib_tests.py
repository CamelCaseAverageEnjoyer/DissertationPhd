"""Проверка на работоспособность установленных библиотек. Среди них:
- kiamastro
    Нет FKIAMToolbox, не может ничего запустить
- nasapy
    Есть класс Nasa, что-то делает
- celmech
    Большая библиотека. Обязательно изучить, когда займёшься уравнениями
- poliastro
    Атмосфера: https://docs.poliastro.space/en/stable/examples/Atmospheric%20models.html
- astropy
"""
from matplotlib import pyplot as plt
import numpy as np

from astropy import units as u
from poliastro.earth.atmosphere import COESA62, COESA76

# We build the atmospheric instances
coesa62 = COESA62()
coesa76 = COESA76()

# Collect all atmospheric models and define their plotting properties
atm_models = {
    coesa62: ["--r", "r", "Coesa 1962"],
    coesa76: ["-b", "b", "Coesa 1976"],
}


fig, axs = plt.subplots(1, 3, figsize=(12, 5))
fig.suptitle("State variables against altitude", fontweight="bold")

# Complete altitude range and initialization of state variables sets
alt_span = np.linspace(0, 1000, 1001) * u.km
T_span = np.array([]) * u.K
p_span = np.array([]) * u.Pa
rho_span = np.array([]) * u.kg / u.m**3

# We solve for each property at given altitude
for alt in alt_span:
    T, p, rho = coesa76.properties(alt)
    T_span = np.append(T_span, T)
    p_span = np.append(p_span, p.to(u.Pa))
    rho_span = np.append(rho_span, rho)

# Temperature plot
axs[0].set_title("Temperature")
axs[0].set_xlabel("T [K]")
axs[0].set_ylabel("Altitude [K]")
axs[0].plot(T_span, alt_span)

# Pressure plot
axs[1].set_title("Pressure")
axs[1].set_xlabel("p [Pa]")
axs[1].plot(p_span, alt_span)
axs[1].set_xscale("log")

# Density plot
axs[2].set_title("Density")
axs[2].set_xlabel(r"$\rho$ [kg/m3]")
axs[2].plot(rho_span, alt_span)
axs[2].set_xscale("log")

plt.show()
