import matplotlib as m
fig_size = (7,6)
font = 'Geneva'
params = {'backend': 'tkagg',
          'axes.labelsize': 10,
          'axes.titlesize': 10,
          'font.size': 10,
          'legend.fontsize': 10,
          'font.family': 'sans-serif',
          'font.sans-serif': font,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': False,
          'figure.figsize': fig_size,
          'figure.dpi': 80}
m.rcParams.update(params)
m.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]
m.rcParams['mathtext.fontset'] = 'custom'
m.rcParams['mathtext.rm'] = font
m.rcParams['mathtext.it'] = font+':italic'
m.rcParams['mathtext.bf'] = font+':bold'


import sparsetools as sp
import matplotlib.pyplot as plt
#import spectral as s

i = sp.profile('observed.nc')
o = sp.profile('synthetic_cycle1.nc')
m = sp.model('atmosout_cycle1.nc')

plt.close("all")

f = plt.figure(figsize=(7,5))
ax1 = plt.subplot2grid((5,6), (0,0), colspan=6, rowspan=3)
ax2 = plt.subplot2grid((5,6), (3,0), colspan=2, rowspan=2)
ax3 = plt.subplot2grid((5,6), (3,2), colspan=2, rowspan=2)
ax4 = plt.subplot2grid((5,6), (3,4), colspan=2, rowspan=2)


ss = 0
#wav = s.air2vac(i.wav)
ax1.plot(i.dat[0,0,0,:,ss],'.', color='black')
ax1.plot(o.dat[0,0,0,:,ss], color='orangered')
ax1.set_ylabel('Intensity')
ax1.set_xlabel(r'Wavelength [$\mathrm{\AA}$]')

if(ss == 0):
    ax1.set_ylim(0.05, 0.9)
elif(ss == 1 or ss == 2):
    ax1.set_ylim(-0.055, 0.005)
else:
    ax1.set_ylim(-0.03, 0.03)

dep = m.cmass[0,0,0]
    
ax2.plot(dep, m.temp.squeeze(), 'k-')
ax2.set_ylim(3000,12000)

ax3.plot(dep, m.vturb.squeeze()*1.e-5, 'k-')
ax3.set_ylim(0,10)

ax4.plot(dep, m.vlos.squeeze()*1.e-5, 'k-')
ax4.set_ylim(-7,7)

f.set_tight_layout(True)
plt.show()

#f.savefig('eb_vissers.pdf', dpi=300, format='pdf', bbox='standard')
