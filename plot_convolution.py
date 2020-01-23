from matplotlib import pyplot as plt
import numpy as np

wavelengths = np.linspace(-0.2, 1.2, 500)

L = np.exp( np.e * wavelengths)
E = np.exp(-np.e * wavelengths)
R = L / E
S = np.zeros_like(wavelengths) ; S[(wavelengths >= 0.3) & (wavelengths <= 0.7)] = 1

fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(7,2.5), sharex=True, tight_layout=True)

axs[0,0].axis("off")

axs[0,1].plot(wavelengths, L, c='k')
axs[0,1].set_title("$L_w$")
axs[0,1].annotate(r"$\div$", xy=(0.5, 0.73), xycoords="figure fraction", size="xx-large", ha="center", va="center")

axs[0,2].plot(wavelengths, E, c='k')
axs[0,2].set_title("$E_d$")
axs[0,2].annotate(r"$=$", xy=(0.745, 0.73), xycoords="figure fraction", size="xx-large", ha="center", va="center")

axs[0,3].plot(wavelengths, R, c='k')
axs[0,3].set_title("$R_{rs}$")

axs[1,0].plot(wavelengths, S, c='k')
axs[1,0].fill_between(wavelengths, 0, S, facecolor="0.7")
axs[1,0].set_title("$S_B$")
axs[1,0].annotate("$ = 1$", xy=(0.63, 0.3), xycoords="axes fraction")

axs[1,1].plot(wavelengths, L * S, c='k')
axs[1,1].fill_between(wavelengths, 0, L * S, facecolor="0.7")
axs[1,1].set_title(r"$L_w \times S_B$")
axs[1,1].annotate(r"$ \approx 5.2$", xy=(0.63, 0.3), xycoords="axes fraction")
axs[1,1].annotate(r"$\div$", xy=(0.5, 0.24), xycoords="figure fraction", size="xx-large", ha="center", va="center")

axs[1,2].plot(wavelengths, E * S, c='k')
axs[1,2].fill_between(wavelengths, 0, E * S, facecolor="0.7")
axs[1,2].set_title(r"$E_d \times S_B$")
axs[1,2].annotate(r"$ \approx 0.3$", xy=(0.63, 0.3), xycoords="axes fraction")
axs[1,2].annotate(r"$\neq$", xy=(0.745, 0.24), xycoords="figure fraction", size="xx-large", ha="center", va="center")

axs[1,3].plot(wavelengths, R * S, c='k')
axs[1,3].fill_between(wavelengths, 0, R * S, facecolor="0.7")
axs[1,3].set_title(r"$R_{rs} \times S_B$")
axs[1,3].annotate(r"$ \approx 42$" + "\n" + r"$ \neq 15$", xy=(0.63, 0.3), xycoords="axes fraction")

for ax in axs.ravel():
    ax.tick_params(axis="both", left=False, labelleft=False, bottom=False, labelbottom=False)
    ax.set_xlabel("$\lambda$")

plt.savefig("convolution.pdf")
plt.show()
