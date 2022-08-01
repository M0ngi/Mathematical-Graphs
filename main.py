# Some graphs from this PDF: https://www.scirp.org/pdf/AM20110500013_23450279.pdf
# Page 3, 4

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

PI = np.pi           # Pi constant to use
POW = np.power       # Power function to use
LN = np.log          # Log function to use


def configFigure(ymin, ymax, step=10):
  """
    Configure the figure display.
  """
  plt.ylim(ymin, ymax) # Set Y display limit from -20 to 80
  plt.xticks(x) # Set X axis values (0, 0.4, 0.8 ...)
  plt.yticks(np.arange(ymin, ymax+1, step)) # Set Y axis values
  plt.xlabel('Q') # Set X axis label
  plt.grid()
  plt.legend()


def legendText(phi, delta):
  """
    Generates a legend text for the given phi & delta
  """
  return "Phi: " + str(phi) + " - Delta: " + str(delta)


def h(phi, z):
  return 1 + phi * np.sin(2 * PI * z)


def g(z, q, phi, delta):
  """
    G function.
    If delta == 0, returns G0 value.
  """
  h_value = h(phi, z)

  if delta != 0:
    top = 8*(q/PI + ( POW(h_value, 2) - 1 - phi/2 ) )
    bot = ( POW(delta, 4) - POW(h_value, 4) - POW(POW(delta, 2) - POW(h_value, 2), 2)/LN(delta/h_value) )
    return -top/bot
  
  # delta == 0, g0 value
  return 8*(q/PI + ( POW(h_value, 2) - 1 - phi/2 ) ) / POW(h_value, 4)
  

# Cette fonction calcule la force de frottement externe
def Dp(q, phi, delta):
  """
    Integral of G(z), from 0 to 1
    phi, delta, q: constants
  """
  return integrate.quad(lambda x: g(x, q, phi, delta), 0, 1)[0]


def f0(z, q, phi, delta):
  """
    Used to calculate F0(Q)
  """
  h_value = h(phi, z)
  if delta != 0:
    return PI*g(z, q, phi, delta)*(POW(h_value, 2) - (POW(delta, 2) - POW(h_value, 2))/(2*LN(delta/h_value)))
  
  # delta == 0
  return PI * POW(h_value, 2) * g(z, q, phi, delta)


# Cette fonction calcule la force de frottement interne
def F0(q, phi, delta):
  """
    Integral of f0(z) from 0 to 1
    phi, delta, q: constants
  """
  return integrate.quad(lambda x: f0(x, q, phi, delta), 0, 1)[0]


def fi(z, q, phi, delta):
  """
    Used to calculate Fi(Q)
  """
  h_value = h(phi, z)
  if delta != 0:
    return PI * g(z, q, phi, delta) * ( POW(delta, 2) - ( POW(delta, 2) - POW(h_value, 2) )/(2*LN(delta/h_value)) )
  
  # Delta == 0
  return 0


# Cette fonction calcule la perte de charge
def Fi(q, phi, delta):
  """
    Integral of fi(z) from 0 to 1
    phi, delta, q: constants
    
    Fi = 0 if delta = 0
  """
  return integrate.quad(lambda x: fi(x, q, phi, delta), 0, 1)[0]

x = np.arange(0, 6.1, 0.4) # X axis values

# Figure 1: Dp(Q)
plt.figure(1)

phi_v   = [0.2, 0.4 ]
delta_v = [0  , 0.44]

y = [Dp(q=i, phi=phi_v[0], delta=delta_v[0]) for i in x]
plt.plot(x, y, 'bv-', label=legendText(phi_v[0], delta_v[0]))

y = [Dp(q=i, phi=phi_v[0], delta=delta_v[1]) for i in x]
plt.plot(x, y, 'r^-', label=legendText(phi_v[0], delta_v[1]))

y = [Dp(q=i, phi=phi_v[1], delta=delta_v[0]) for i in x]
plt.plot(x, y, 'g1-', label=legendText(phi_v[1], delta_v[0]))

y = [Dp(q=i, phi=phi_v[1], delta=delta_v[1]) for i in x]
plt.plot(x, y, 'mx-', label=legendText(phi_v[1], delta_v[1]))

plt.ylabel('Dp')
configFigure(ymin=-20, ymax=80)

# Figure 2: F0(Q)
plt.figure(2)

phi_v   = [0.2, 0.4 ]
delta_v = [0  , 0.44]

y = [F0(q=i, phi=phi_v[0], delta=delta_v[0]) for i in x]
plt.plot(x, y, 'bv-', label=legendText(phi_v[0], delta_v[0]))

y = [F0(q=i, phi=phi_v[0], delta=delta_v[1]) for i in x]
plt.plot(x, y, 'r^-', label=legendText(phi_v[0], delta_v[1]))

y = [F0(q=i, phi=phi_v[1], delta=delta_v[0]) for i in x]
plt.plot(x, y, 'g1-', label=legendText(phi_v[1], delta_v[0]))

y = [F0(q=i, phi=phi_v[1], delta=delta_v[1]) for i in x]
plt.plot(x, y, 'mx-', label=legendText(phi_v[1], delta_v[1]))

plt.ylabel('F0')
configFigure(ymin=-10, ymax=40)

# Figure 3: Fi(Q)
plt.figure(3)

phi_v   = [0.2 , 0.4 ]
delta_v = [0.32, 0.44]

y = [Fi(q=i, phi=phi_v[0], delta=delta_v[0]) for i in x]
plt.plot(x, y, 'bv-', label=legendText(phi_v[0], delta_v[0]))

y = [Fi(q=i, phi=phi_v[0], delta=delta_v[1]) for i in x]
plt.plot(x, y, 'r^-', label=legendText(phi_v[0], delta_v[1]))

y = [Fi(q=i, phi=phi_v[1], delta=delta_v[0]) for i in x]
plt.plot(x, y, 'g1-', label=legendText(phi_v[1], delta_v[0]))

y = [Fi(q=i, phi=phi_v[1], delta=delta_v[1]) for i in x]
plt.plot(x, y, 'mx-', label=legendText(phi_v[1], delta_v[1]))

plt.ylabel('Fi')
configFigure(ymin=-20, ymax=30, step=5)

plt.show()
