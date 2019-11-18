import subprocess, csv
import numpy as np
import matplotlib.pyplot as p
import astropy.constants as cs
from scipy import integrate


# Basic Functions
def convert_radec(RAh, RAm, RAs, Decd, Decm, Decs):
    RA = (RAh/24. + RAm/24./60. + RAs/24./60./60.)*360.
    Dec = Decd + Decm/60. + Decs/60./60.
    return RA, Dec

# Loading tables
def download_tepcat(save_dir='./'):

    import urllib2
    response = urllib2.urlopen('http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-csv.csv')
    allplanets = response.read()

    response = urllib2.urlopen('http://www.astro.keele.ac.uk/jkt/tepcat/observables.csv')
    observables = response.read()

    with open(save_dir+'allplanets-tepcat.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',)
        for line in allplanets.split('\n'):
            if not line.strip() == '':
                writer.writerow([ el.strip() for el in line.split(',') ])
    print('Downloaded system parameters, \"Well-studied transiting planets\"')
            
    with open(save_dir+'observables_tepcat.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',)
        for line in observables.split('\n'):
            if not line.strip() == '':
                writer.writerow([ el.strip() for el in line.split(',') ])
    print('Downloaded magnitudes, \"For planning observations\"')

# Planet Dictionary
def get_planet(name, od):
    try:
        ind = np.min(np.nonzero(od['System']==name)[0]) 
        planet = { var: arr[ind] for (var, arr) in od.items() }
    except TypeError: # wrong way around?
        name, od = od, name
        ind = np.min(np.nonzero(od['System']==name)[0])
        planet = { var: arr[ind] for (var, arr) in od.items() }
    return planet

def print_planet(name, od):
    
    planet = get_planet(name, od)
    keys = planet.keys()
    lkey = max([len(key) for key in keys])
    keys.remove('System')
    print '-'*len(planet['System'])
    print planet['System']
    print '-'*len(planet['System'])
    
    keyorder = ['M_b', 'R_b', 'g_b', 'Teq', 'TEQ', 'SH', 'Vmag', 'Kmag', 'a(AU)', 'Period', 'e', 
                    'M_A', 'R_A', 'g', 'Teff', '[Fe/H]', 'RA', 'Dec']
    for key in keyorder:
        if key in keys:
            print key.rjust(lkey), '\t', planet[key]
    for key in keys: 
        if not key in keyorder:
            print key.rjust(lkey), '\t', planet[key]

def read_exodata(csv_file='allplanets-tepcat.csv'):
    '''
    Read in csv file from TEPCat
    First read and output header as a list
    
    Then read each line in csv, splitting and putting in a dict
    using the header values as keys
    '''
    with open(csv_file, 'r') as g:
        reader = csv.reader(g, delimiter=',')
        header = [ name.strip() for name in reader.next() ]
        #units = reader.next()
        yield header
        for row in reader:
            if row[0].startswith('OGLE'): continue
            store = []
            for el in row:
                el = el.strip()
                try: el = float(el)
                except ValueError: pass #leave string
                store.append(el)
            row = store
        
            data_dict = dict(list(zip(header, row)))
            yield data_dict

def get_arrays(exo_data, hnames=['System', 'Teff', '[Fe/H]', 'M_A', 'R_A', 'loggA', 'rho_A', 'Period', \
                           'e', 'a(AU)', 'M_b', 'R_b', 'g_b', 'rho_b', 'Teq', \
                           'Discovery_reference', 'Recent_reference'], \
              need=['System']):
    '''
    Returns a dictionary of the requested parameter arrays, 
    dictionary contains all the parameters listed in hnames
    Only need to change hnames if changing the format of the csv 
    or using a different catalogue
    
    exo_data: read_exodata(csv_file) - list of dictionaries
    need: specifies parameters required, if a planet lacks a parameter in need 
                it is skipped and the name is placed in rejects
                
    Can put some other constraints on planet params here, e.g. Mass > 0
    but TEPCat is pretty solid. Some stuff like TEQ has -1
    '''
    rejects = []
    
    out_dict = {}
    for i in range(len(hnames)):
        out_dict[hnames[i]] = []
        
    for line in exo_data:
        # check that we have all hnames
        need_vals = [ line[hname] for hname in need ]
        #if np.any([ val == -1. for val in need_vals]): 
        #    rejects.append(line['System'])
        #    continue
        try:
            vals = [ line[hname] for hname in hnames ]
        except KeyError: # value missing
            continue

        for i in range(len(hnames)): out_dict[hnames[i]].append(vals[i])
            
    for i in range(len(hnames)):
        out_dict[hnames[i]] = np.array(out_dict[hnames[i]])
      
    #print 'Rejects:', rejects
    return out_dict

def get_tepcat_planets(bands_2mass=['J','H','K'], TESS=False):
    '''
    Combine TEPCatalogues of well studied planets and observables
    basically just to include V and K mags in dictionary
    
    Just use this to read in the data
    '''
    exo_data = read_exodata('allplanets-tepcat.csv')
    header = exo_data.next()
    od = get_arrays(exo_data, need=['System', 'M_b', 'R_b', 'Teq', 'a(AU)', 'R_A', 'Period'])
    
    exo_data2 = read_exodata('observables_tepcat.csv')
    header = exo_data2.next()
    od2 = get_arrays(exo_data2, hnames=['System', 'Vmag', 'Kmag', 'RAh', 'RAm', 'RAs', 'Decd', 'Decm', 'Decs'])
    
    od2_hnames = ['System', 'Vmag', 'Kmag', 'RA', 'Dec']
    RA, Dec = convert_radec(od2['RAh'], od2['RAm'], od2['RAs'], od2['Decd'], od2['Decm'], od2['Decs'])
    od2['RA'] = RA; od2['Dec'] = Dec
    
    for hname in od2_hnames[1:]: 
        od[hname] = np.ones(len(od['System']))*-1
    
    od_names = list(od['System'])
    od2_names = list(od2['System'])
    for i1, name in enumerate(od['System']):
        try: 
            i2 = od2_names.index(name)
            for hname in od2_hnames[1:]:
                od[hname][i1] = od2[hname][i2]
            
        except ValueError: 
            pass
        
    # Add info on TESS observations
    if TESS:
        with open('TESS_obs.dat', 'r') as g: out = g.read()
        TICs, Sectors = zip(*[ [int(lin.split('|')[ival]) for ival in [0,5]] for lin in out.split('\n')[1:-1]])
        TICs, Sectors = np.array(TICs), np.array(Sectors)
        _sectors = []; _TESS_days = []
        for i, planet in enumerate(od['System']):
            secs = Sectors[TICs==i]
            _sectors.append(secs)
            _TESS_days.append(27*len(secs))
        od['TESS_days'] = np.array(_TESS_days)
        od['TESS_sectors'] = np.array(_sectors)
    
    return od


# TESS
# or use website https://heasarc.gsfc.nasa.gov/cgi-bin/tess/webtess/wtv.py
def get_TESS_obs(od):
    tempfname = '_coords.dat'
    # Make a file of all RA, Decs
    RA, Dec = od['RA'], od['Dec']
    
    with open(tempfname, 'w') as g:
        for i,r,d in zip(np.arange(len(RA)),RA,Dec):
            g.write('{}\t{}\t{}\n'.format(i,r,d))
    
    pr = subprocess.Popen("python -m tess_stars2px -f /home/jacob/Exo_Data/{}".format(tempfname), 
                                                                   stdout=subprocess.PIPE, shell=True)
    out, err = pr.communicate()
    
    with open('TESS_obs.dat', 'w') as g:
        g.write(out)
    return out

def update_TESS():
    # Should update everytime new planets are added, otherwise they will be missed
    od = get_tepcat_planets(TESS=False)
    out = get_TESS_obs(od)

# Plotting 

from bokeh.io import output_notebook
output_notebook()
import matplotlib as mpl
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import ColumnDataSource, Range1d, LabelSet, Label, HoverTool, LinearColorMapper, ColorBar
from bokeh.palettes import magma

bokehpalette = magma(256)

def plot_html(od, xname='M_b', yname='R_b', zname=None, z_range=None, size=8, fname=None, 
              x_range=None, y_range=None, title=None, reverse_z=False, plot=False,
              keys=['M_b', 'R_b', 'Teq', 'Period', 'Vmag', 'Kmag', 'M_A', 'R_A', 'Teff', 'e']):
    
    
    if zname is not None:
        _zs = od[zname].copy()
        _zs[np.isnan(od[zname])] = 0.
        if z_range:
            _zs[od[zname]<z_range[0]] = z_range[0]
            _zs[od[zname]>z_range[1]] = z_range[1]
        cmap = mpl.cm.ScalarMappable(cmap='magma')
        if reverse_z: _zs = -_zs
        colors = cmap.to_rgba(_zs, bytes=True)
        colors = ["#{0:02x}{1:02x}{2:02x}".format(c[0],c[1],c[2]) for c in colors]
        od['colors'] = colors
    else:
        od['colors'] = ['blue']*len(od['System'])
        
    source = ColumnDataSource(data=od)
    tools = "pan,wheel_zoom,box_zoom,reset,save".split(',')
    
    keys = keys
    if xname in keys: keys.remove(xname)
    if yname in keys: keys.remove(yname)

    tooltips=[ ("Name", "@System"),(xname+':', "@"+xname),(yname+':', "@"+yname)]
    tooltips += [ (key+':', '@'+key) for key in keys]
    
    hover = HoverTool(tooltips=tooltips)
    tools.append(hover)
        
    plt = figure(tools=tools, title=title, x_range=x_range, 
                 y_range=y_range)

    plt.scatter(x=xname, y=yname, source=source, size=size, line_color='black', fill_color='colors')
    plt.xaxis.axis_label = xname
    plt.yaxis.axis_label = yname
    
    if zname is not None:
        mapper = LinearColorMapper(palette=bokehpalette, 
                                   low=np.nanmin(_zs), high=np.nanmax(_zs))
        color_bar = ColorBar(color_mapper=mapper, label_standoff=12, border_line_color=None,
                            location=(0,0))
        color_bar_label = Label(text=zname, angle=90*np.pi/180.,
                                text_color='black', text_font_style='normal',
                                text_font_size='12px')#, x=5,y=5), 
                                #x_units='screen', y_units='screen')
        
        plt.add_layout(color_bar, 'right')
        #plt.add_layout(color_bar_label, 'right')
        
    if fname: save(plt, fname)
    if plot: show(plt)

# Tables 

def csv_arrays(fname, od, keys, names=None, units=['-'], formats=['{}'], sort='System', reverse=True, finite=True, nomin1=True):
    if names is None: names=keys
    with open(fname, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        names_dict = {}
        for key, val in zip(keys, names): names_dict[key] = val
        writer.writerow(names_dict)
        units_dict = {}
        for key, val in zip(keys, units): units_dict[key] = val
        writer.writerow(units_dict)
        
        indexes = np.argsort(od[sort])
        if reverse: indexes = reversed(indexes)
            
        for i in indexes:
            name = od['System'][i]
            planet = get_planet(name, od)
            
            # check values make sense, otherwise skip
            vals = [val for val in planet.values() if not type(val) in [str, np.string_]]
            
            if finite: # ensure valus are finite
                if not np.all([np.isfinite(val) for val in vals if type(val) is not str]): continue 
                        
            if nomin1: # ensure values are not = -1
                if not np.all([val != -1 for val in vals if type(val) is not str]): continue 
                       
            row = {}
            for key, form in zip(keys, formats): row[key] = form.format(planet[key])
            writer.writerow(row)


# For transmission

def calc_eq_temp(a, Ts, Rs, A=0, f_redist=2.):
    '''Assume parameters in SI units'''
    # f_redist = 1 for no redistribution
    # f_redist = 2 for full day-night redistribution
    return np.sqrt(np.divide(Rs,a)/f_redist)*Ts*(2**0.25) * (1-A)**0.25

def calc_surface_g(M, R):
    '''Returns surface g in SI'''
    return cs.G.value*M/R**2

def calc_scale_h(T,g,mu=1.2):
    '''
    Calculate atmospheric scale height,
    data should be in SI units.
    '''
    return cs.k_B.value*T/g/(mu*cs.m_p.value)

def calc_transit_depth(sh, Rp, Rs, nsh=1):
    return 2*nsh*sh*Rp/(Rs**2)

def weins_peak(T):
    return cs.b_wien.value*T*1e6 # in micron

def calc_impact_parameter(a, i, Rs):
    return a / Rs * np.cos(i)

def calc_transit_duration(Rp, Rs, a, i, Per):
    b = calc_impact_parameter(a, i, Rs)
    return Per/np.pi * np.arcsin(((Rp+Rs)**2-(b*Rs)**2)**0.5 / a)

# For emission


def bb(l, T):
    # the energy per unit time (or the power) radiated per unit area of emitting surface 
    # in the normal direction per unit solid angle by a black body at temperature T at wavelength l
    return 2*cs.h.value*cs.c.value**2 / l**5 / (np.exp(cs.h.value*cs.c.value/l/cs.k_B.value/T)-1)

def bb_ph(l, T):
    # the PHOTONS per unit time radiated per unit area of emitting surface 
    # in the normal direction per unit solid angle by a black body at temperature T at wavelength l
    # Energy per photon is hc/lambda
    return bb(l, T) / (cs.h.value*cs.c.value / l)


def band_wave(band):
    if band == 'V':
        l1, l2 = (551-88)*10**-9, (551+88)*10**-9
    elif band == 'B':
        l1, l2 = (445-47)*10**-9, (445+47)*10**-9
    elif band == 'J':
        l1, l2 = (1.235-0.081)*10**-6, (1.235+0.081)*10**-6
    elif band == 'H':
        l1, l2 = (1.662-0.1255)*10**-6, (1.662+0.1255)*10**-6
    elif band == 'K':
        l1, l2 = (2.159-0.131)*10**-6, (2.159+0.131)*10**-6
    elif band == 'G102':
        l1, l2 = .800*10**-6, 1.150*10**-6
    elif band == 'G141':
        l1, l2 = 1.075*10**-6, 1.700*10**-6
    elif band == 'Spitzer1':
        l1, l2 = 3.1*10**-6, 4.1*10**-6
    else:
        l1, l2 = band
    return l1, l2


def bb_flux(T, R, band=(1e-6,2e-6)):
    # calculate flux in units of J/s in a given band
    l1, l2 = band_wave(band)
    I, error = integrate.quad(bb, l1, l2, args=T)
    F = I*np.pi*R**2 * np.pi #pi solidangle, piR^2 area
    return F

def bb_photon_flux(T, R, band=(1e-6,2e-6)):
    # calculate flux in units of J/s in a given band
    l1, l2 = band_wave(band)
    I, error = integrate.quad(bb_ph, l1, l2, args=T)
    F = I*np.pi*R**2 * np.pi #pi solidangle, piR^2 area
    return F


def eclipse_depth(Tp, Ts, Rp, Rs, band):
    try:
        Fp = np.array([ bb_flux(tp, rp, band) for tp, rp in zip(Tp, Rp)])
        Fs = np.array([ bb_flux(ts, rs, band) for ts, rs in zip(Ts, Rs)])
    except TypeError: # not iterable
        Fp = bb_flux(Tp, Rp, band)
        Fs = bb_flux(Ts, Rs, band) # J/s
    return Fp/Fs

import pysynphot as S
def eclipse_depth_phoenix(Tp, Ts, Rp, Rs, Ms, band):
    # WIP! Does not make sense
    gs = calc_surface_g(Ms, Rs) # surface gravity of the star

    sp = S.Icat('k93models', *stellar_params)
    phoenix_wvs, phoenix_fl = sp.wave, sp.flux*1e-7*100**2 *np.pi # in Angstrom, W/m2/A

    model_wv = np.array([np.mean(wbin) for wbin in wv_bins])
    binsizes = np.array([wbin[1]-wbin[0] for wbin in wv_bins]) # in micron

    spectrum_obj = S.ArraySpectrum(phoenix_wvs, phoenix_fl, waveunits='angstrom')
    def phoenix_spec(l):
        # Interpolate stellar spectrum properly, accounting for change in resolution
        # l in micron, has to be in Angstrom (pysynphot just breaks in micron for no reason)
        # return spec in W/m2/A
        l = np.array(l)*1e4
        filt = S.spectrum.ArraySpectralElement(phoenix_wvs, np.ones_like(phoenix_wvs))
        obs = S.observation.Observation(spectrum_obj, filt, binset=l, force='taper')
        flux = obs.binflux
        return flux # return flux in W/m2/A

    fs = phoenix_spec(model_wv) * binsizes*1e4 # should really integrate here
    fp = bb(model_wv*1e-6, Tp) * binsizes*1e-6 * np.pi # in W/m^2, need to integrate over solid angle (pi)
    fpfs = fp/fs * rp**2

    return fpfs
    

def reflected_depth(Rp,a,A):
    # Fp/Fs is just geometry in this case
    return (Rp/a/2)**2 *A

def emission_over_reflection(Tp, Ts, Rp, Rs, a, l, A=1):
    # calculate the emitted/relected flux for given parameters at a given wavelength
    # Note: albedo for reflected component may not be the same as that used 
    # to calculate eq temperature.
    try:
        Fp_ref = (Rp/a/2)**2 * A * np.array([ bb(l,ts) for ts in Ts ])
        Fp_em = np.array([bb(l,tp) for tp in Tp])
    except TypeError:
        Fp_ref = (Rp/a/2)**2 * A * bb(l,Ts)
        Fp_em = bb(l,Tp)
    return Fp_em/Fp_ref

def reflected_frac(Tp, Ts, Rp, Rs, a, l, A=1):
    # calculate the fraction of relected flux for given parameters at a given wavelength
    # Note: albedo for reflected component may not be the same as that used 
    # to calculate eq temperature.
    try:
        Fp_ref = (Rp/a/2)**2 * A * np.array([ bb(l,ts) for ts in Ts ])
        Fp_em = np.array([bb(l,tp) for tp in Tp])
    except TypeError:
        Fp_ref = (Rp/a/2)**2 * A * bb(l,Ts)
        Fp_em = bb(l,Tp)
    return Fp_ref/(Fp_ref+Fp_em)
