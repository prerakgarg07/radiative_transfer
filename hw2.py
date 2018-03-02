import numpy as np
import astropy.units as u
from matplotlib import pyplot as plt

# Assumed Values
n = 1*u.cm**-3
D = 100*u.pc
I0 = 1e-9 *u.erg*u.s**-1*u.cm**-2*u.sr**-1*u.Hz**-1
Snu = 0.2e-9 *u.erg*u.s**-1*u.cm**-2*u.sr**-1*u.Hz**-1
sigma = 2e-25*u.m**2
del_s = 0.1*u.pc
fwhm = 2
freq = np.arange(10, 20, 0.1)

# Part a
def find_tau(sigma, s):
    return (n*sigma*s).decompose()


# Part b
def Inu(I_0,tau,Snu):
    return I_0*np.exp(-tau) + Snu*(1-np.exp(-tau))


def get_Inu(tau,I0):
    I = []
    del_tau = tau/100
    I.append(I0.value)
    for i in range(100):
        I.append(Inu(I[i]*u.erg*u.s**-1*u.cm**-2*u.sr**-1*u.Hz**-1,del_tau,Snu).value)
    return I[100]

# Part c
def sigma_nu_gauss(max):
    return max*np.exp((-4*np.log(2)*(freq-15)**2)/fwhm)*u.m**2


#Part d
# i.
tau = find_tau(sigma_nu_gauss(2e-4), del_s)
print tau
Intensity = []
for tau_nu in tau:
    Intensity.append(get_Inu(tau_nu,I0))
plt.plot(Intensity)
plt.xlabel("Frequency")
plt.ylabel("Intensity")
plt.show()


# ii
tau = find_tau(sigma_nu_gauss(2e-20), del_s)
Intensity = []
I0 = 0 * u.erg*u.s**-1*u.cm**-2*u.sr**-1*u.Hz**-1
for tau_nu in tau:
    Intensity.append(get_Inu(tau_nu,I0))
plt.plot(Intensity)
plt.xlabel("Frequency")
plt.ylabel("Intensity")
plt.show()

#iii
tau = find_tau(sigma_nu_gauss(2e-22), del_s)
Intensity = []
print tau
I0 = 0.05e-9 * u.erg*u.s**-1*u.cm**-2*u.sr**-1*u.Hz**-1
for tau_nu in tau:
    Intensity.append(get_Inu(tau_nu,I0))
plt.plot(Intensity)
plt.xlabel("Frequency")
plt.ylabel("Intensity")
plt.show()

#iv
tau = find_tau(sigma_nu_gauss(2e-22), del_s)
Intensity = []
I0 = 1e-9 * u.erg*u.s**-1*u.cm**-2*u.sr**-1*u.Hz**-1
for tau_nu in tau:
    Intensity.append(get_Inu(tau_nu,I0))
plt.plot(Intensity)
plt.xlabel("Frequency")
plt.ylabel("Intensity")
plt.show()