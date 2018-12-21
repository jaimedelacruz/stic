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

f = plt.figure(figsize=(7,7))
ax1 = plt.subplot2grid((7,6), (0,0), colspan=6, rowspan=3)
ax2 = plt.subplot2grid((7,6), (3,0), colspan=2, rowspan=2)
ax3 = plt.subplot2grid((7,6), (3,2), colspan=2, rowspan=2)
ax4 = plt.subplot2grid((7,6), (3,4), colspan=2, rowspan=2)
ax5 = plt.subplot2grid((7,6), (5,0), colspan=2, rowspan=2)
ax6 = plt.subplot2grid((7,6), (5,2), colspan=2, rowspan=2)
ax7 = plt.subplot2grid((7,6), (5,4), colspan=2, rowspan=2)

ss = 0
#wav = s.air2vac(i.wav)
ax1.plot(i.dat[0,0,0,:,ss],'.', color='black')
ax1.plot(o.dat[0,0,0,:,ss], color='orangered')
ax1.set_ylabel('Intensity')
ax1.set_xlabel(r'Wavelength [index]')

if(ss == 0):
    ax1.set_ylim(0.05, 0.9)
elif(ss == 1 or ss == 2):
    ax1.set_ylim(-0.055, 0.005)
else:
    ax1.set_ylim(-0.03, 0.03)

dep = m.ltau[0,0,0]
    
ax2.plot(dep, m.temp.squeeze()*1.e-3, 'k-')
ax2.set_ylim(3,12)
ax2.set_title('Temp [kK]')
ax2.set_xlabel(r'Log $\tau_{500}$')

ax3.plot(dep, m.vturb.squeeze()*1.e-5, 'k-')
ax3.set_ylim(0,10)
ax3.set_title('Vturb [km/s]')
ax3.set_xlabel(r'Log $\tau_{500}$')

ax4.plot(dep, m.vlos.squeeze()*1.e-5, 'k-')
ax4.set_ylim(-7,7)
ax4.set_title('Vlos [km/s]')
ax4.set_xlabel(r'Log $\tau_{500}$')

ax5.plot(dep, m.Bln.squeeze(), 'k-')
ax5.set_ylim(-1000,1000)
ax5.set_title('B_long [G]')
ax5.set_xlabel(r'Log $\tau_{500}$')

ax6.plot(dep, m.Bho.squeeze(), 'k-')
ax6.set_ylim(0,1000)
ax6.set_title('|B_hor| [G]')
ax6.set_xlabel(r'Log $\tau_{500}$')

ax7.plot(dep, m.azi.squeeze() * 180/3.1415926, 'k-')
ax7.set_ylim(0,180)
ax7.set_title('B_azi [deg]')
ax7.set_xlabel(r'Log $\tau_{500}$')



#f.set_tight_layout(True)
f.subplots_adjust(wspace=0.75, hspace=1.4, left=0.07, right=0.97, top=0.98, bottom=0.08)
f.savefig('fig_result.pdf', format='pdf', bbox_inches='tight')
#plt.ion()
plt.show()

