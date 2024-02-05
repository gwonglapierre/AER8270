import numpy as np
from matplotlib import pyplot as plt
import HSPM
import geometryGenerator as gg

plt.rcParams['text.usetex'] = True


# Creation de la geometrie d'un profil NACA 4412

panels = gg.GenerateNACA4digit(4.,4.,12.,100)
alpha = np.linspace(-20,20,9)
prob = HSPM.HSPM(panels,alpha)

prob.run()

CD = prob.CD
CL = prob.CL
CM = prob.CM

print(CD)
print(CL)
print(CM)

# Afficher les résultats

fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2,2)

ax1.plot(alpha,CL)
ax1.set_title('$C_L$ vs. $\alpha$')
ax1.set_xlabel('$\alpha$')
ax1.set_ylabel('$C_L$')
ax3.plot(alpha,CD)
ax4.plot(alpha,CM)
ax1.grid()
plt.show()