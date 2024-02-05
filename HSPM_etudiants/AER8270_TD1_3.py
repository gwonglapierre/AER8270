import numpy as np
from matplotlib import pyplot as plt
import HSPM
import geometryGenerator as gg
import shutil, sys
import os


plt.rcParams['text.usetex'] = True


# Creation de la geometrie d'un profil NACA 6412



panels = gg.ReadPoints("naca6412.txt")
alpha = np.linspace(-20,20,9)
prob = HSPM.HSPM(panels,alpha,[0.,0.,0.])

prob.run()

CD = prob.CD
CL = prob.CL
CM = prob.CM

# Afficher les résultats

fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2,2)
fig.suptitle("Coefficients aérodynamiques d'un profil NACA 6412 par la méthode HSPM", fontsize=16)
fig.set_figwidth(8)
fig.set_figheight(10)

ax1.plot(alpha,CL)
ax1.set_title(r'$C_L$ vs. $\alpha$')
ax1.set_xlabel(r'$\alpha$')
ax1.set_ylabel(r'$C_L$')

ax3.plot(alpha,CD)
ax3.set_title(r'$C_D$ vs. $\alpha$')
ax3.set_xlabel(r'$\alpha$')
ax3.set_ylabel(r'$C_D$')

ax4.plot(alpha,CM)
ax4.set_title(r'$C_M$ vs. $\alpha$')
ax4.set_xlabel(r'$\alpha$')
ax4.set_ylabel(r'$C_M$')

ax2.axis('off')

ax1.grid()
ax3.grid()
ax4.grid()

plt.savefig('AER8270_TD1_3g.png',dpi=300)
plt.show()