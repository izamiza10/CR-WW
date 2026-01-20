from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
pp   = PdfPages('fig302.pdf')
tmp1 = plt.figure(1)
tmp1.set_size_inches(8.00,6.00)
plot = open('fig302-0.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
vale = [float(x[2]) for x in plot]
plt.hist( valx, vale, weights = valy, histtype='step', color='blue', label=r"SK I")
plot = open('fig302-1.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
vale = [float(x[2]) for x in plot]
plt.hist( valx, vale, weights = valy, histtype='step', color='red', label=r"SK II")
plt.xlim( 0.000e+00, 4.000e+02)
plt.ylim( 0.000e+00, 0.000e+00)
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,3))
plt.legend(frameon=False,loc='best')
plt.title(r"Reconnection rate as a function of CM energy for e$^+$e$^-$ $\to$ W$^+$W$^-$")
plt.xlabel(r"$E_{\mathrm{CM}}$ (GeV)")
plt.ylabel(r"Probability")
pp.savefig(tmp1,bbox_inches='tight')
plt.clf()
pp.close()
