def DOS_grab(homepath, norm=True): #homepath is the file path for the DOSCARs
    '''
    Processes VASP DOSCAR format to make for easy plotting (must be LORBIT=True).

    Parameters:
    homepath (str): Path to the DOSCAR file.
    norm (bool): default(True) is to normalise energies to the Fermi level.
    lobster : say whether we're using a lobster file or not

    Returns:
    numpy array: DOS[0,0] = energies
                 DOS[0,1] = total up-spin DOS
                 DOS[0,2] = total down-spin DOS
                 DOS[j,0] = sum of jth site-projected up-spin
                 DOS[j,1] = sum of jth site-projected down-spin
                 DOS[j,k=2:12] = jth site-projected spin, even k = up, odd k = down
                                 order of k: $s,p_x,p_y,p_z,d_xy,d_yz$,d_z^2,d_xz,d_{x^2-y^2}$
                                 s orbitals: 2~k~3
                                 p orbitals: 4~k~9
                                 d orbitals: 10~k~19
    fermi (float): Fermi level in eV
    NEDOS (int): number of energy levels at which DOS was calculated
    '''
    import os
    import numpy as np

    with open(homepath,'r') as file:
        lines = file.readlines()
        centres = int(lines[0].split()[0])
        value_number = int(lines[5].split()[2])
        fermi = float(lines[5].split()[3])
        print('centres: %i\nNEDOS: %i\nFermi level \ eV: %f' % (centres,value_number,fermi))

        n = centres + 1 # number of centres you're grabbing
        '''m = len(lines[-1].split())+1 # +2 for sums at start of jth element, -1 as line contains energy.'''
        # m is number of spin-orbitals you want to grab, 2 for site totals, 18 for s,p,d orbitals
        '''DOS = [[[] for i in range(m)] for i in range(n)] # initialise empty list n x m'''
        DOS = [[[],[],[],[],[]]] # start with the list for the DOS[0] stuff

        if norm:
            for i in range(0,value_number):
                DOS[0][0].append(float(lines[6+i].split()[0]) - fermi) #energies, normalised Fermi
                DOS[0][1].append(float(lines[6+i].split()[1])) # total up spin
                DOS[0][2].append(-1*float(lines[6+i].split()[2])) # total down spin
                DOS[0][3].append(float(lines[6+i].split()[3]))  # integrated up
                DOS[0][4].append(-1*float(lines[6+i].split()[4])) # integrated down
        else:
            for i in range(0,value_number):
                DOS[0][0].append(float(lines[6+i].split()[0])) # raw energies 
                DOS[0][1].append(float(lines[6+i].split()[1])) # total up spin
                DOS[0][2].append(-1*float(lines[6+i].split()[2])) # total down spin
                DOS[0][3].append(float(lines[6+i].split()[3]))  # integrated up
                DOS[0][4].append(-1*float(lines[6+i].split()[4])) # integrated down

    # get the angular momentum-resolved data from the DOSCAR for the d-orbitals:
    # j = site index, i = energy index, k = orbital index (alternates up spin, down spin)
    for j in range(1,centres+1):
        m = len(lines[6+j*(value_number+1)].split()) # update the length to number of orbitals at site j
        DOS.append([ [] for k in range(m+1) ]) # add an extra list for each new site
        for i in range(0,value_number):
            DOS[j][0].append(sum(np.array(lines[6+i+j*(value_number+1)].\
                                         split()[1::2]).astype(float))) # sum of site-projected up
            DOS[j][1].append(sum(np.array(lines[6+i+j*(value_number+1)].\
                                         split()[2::2]).astype(float))) # sum of site-projected down
            for k in range(2,m+1):
                DOS[j][k].append(float(lines[6+i+j*(value_number+1)].split()[k-1]))

    for j in range(1,len(DOS)):
        for k in range(1,len(DOS[j]),2):
            DOS[j][k] = np.array(DOS[j][k])*-1 # make the down spins have negative DOS   
        for k in range(0,len(DOS[j]),2):
            DOS[j][k] = np.array(DOS[j][k])

    return(DOS,fermi,value_number)

def DOSplotter_1atom(DOS,atom=1,emph=None):
    """Plots (atom)'s d-orbital and total DOS from DOS_grab output
    Parameters:
    (array) : output from DOS_grab() function
    atom=(int) : index of atom you want to plot
    emph=(int) or None : emphasise a particular line, 0:xy,1:yz,2:z2,3:xz,4:x2-y2

    Returns:
    fig object - use fig.suptitle to add a title
    axs object - use axs[0 or 1].set to change things like xlim, ylim, titles
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from cycler import cycler

    fig, axs = plt.subplots(1,2, figsize=(16,14))
    fig.tight_layout(rect=[0, 0.03, 1, 0.92])
    plt.rc('font',size=16)

    labels = ["$s$","$p_x$","$p_y$","$p_z$","$xy$","$yz$","$z^2$",\
              "$xz$","$x^2-y^2$"] # orbital labels in VASP order
    default_cycler = (cycler(color=['green', 'green', 'r', 'r','k','k','m','m','b','b']))
    # this cycles through 5 different colours with the same colour for up and down spins
    plt.rc('lines', linewidth=1.0)
    plt.rc('axes', prop_cycle=default_cycler)

    u1,d1 = [],[]

    for i in range(10,20):
        if i % 2 == 0: # checks if odd so that I can plot two-by-two the up and down spins
            u1.append(axs[0].plot(DOS[atom][i],DOS[0][0],label=labels[i//2-1], linewidth=0.75))
            d1.append(axs[0].plot(DOS[atom][i+1],DOS[0][0], linewidth=0.75))

    axs[1].plot(DOS[0][1],DOS[0][0]) # plot the sums of the DOS
    axs[1].plot(DOS[0][2],DOS[0][0])

    for i in axs:
        i.axhline(0, linestyle='-', color = 'k', linewidth = 0.5)
        i.axvline(0,linestyle='-', color = 'k', linewidth = 0.5)
        i.set(xlabel="Density of states"\
                  ,xlim=(-1.5,1.5),ylim=(-8,7)\
                  )
        i.legend(loc='upper right')

    axs[1].set(xlim=(-10,10))
    axs[0].set(ylabel="Energy relative to Fermi level/ eV")

    # changes which d-orbitals you want to emphasise, if any
    if not emph:
        pass
    else:

        u1[emph][0].set(linewidth=2.0)
        d1[emph][0].set(linewidth=2.0)

    return(fig,axs)

def DOSplotter_2atom(DOS,atom1=1,atom2=2,emph=None):
    """Plots (atom)'s d-orbital and total DOS from DOS_grab output

    Parameters:
    (array) : output from DOS_grab() function
    atom1=(int) : index of first atom you want to plot
    atom2=(int) : index of second atom you want to plot
    emph=(int) or None : emphasise a particular line, 0:xy,1:yz,2:z2,3:xz,4:x2-y2

    Returns:
    fig object - use fig.suptitle to add a title
    axs object - use axs[0 or 1].set to change things like xlim, ylim, titles
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from cycler import cycler

    fig, axs = plt.subplots(1,2, figsize=(16,14))
    fig.tight_layout(rect=[0, 0.03, 1, 0.92])
    plt.rc('font',size=16)

    labels = ["$s$","$p_x$","$p_y$","$p_z$","$xy$","$yz$","$z^2$",\
              "$xz$","$x^2-y^2$"] # orbital labels in VASP order
    default_cycler = (cycler(color=['g', 'g', 'r', 'r','k','k','m','m','b','b']))
    # this cycles through 5 different colours with the same colour for up and down spins
    plt.rc('lines', linewidth=1.0)
    plt.rc('axes', prop_cycle=default_cycler)

    u1,d1,u2,d2 = [],[],[],[]

    for i in range(10,20):
        if i % 2 == 0: # checks if odd so that I can plot two-by-two the up and down spins
            u1.append(axs[0].plot(DOS[atom1][i],DOS[0][0],label=labels[i//2-1], linewidth=0.75))
            d1.append(axs[0].plot(DOS[atom1][i+1],DOS[0][0], linewidth=0.75))
            u2.append(axs[1].plot(DOS[atom2][i],DOS[0][0],label=labels[i//2-1], linewidth=0.75))
            d2.append(axs[1].plot(DOS[atom2][i+1],DOS[0][0], linewidth=0.75))

    for i in axs:
        i.axhline(0, linestyle='-', color = 'k', linewidth = 0.5)
        i.axvline(0,linestyle='-', color = 'k', linewidth = 0.5)
        i.set(xlabel="Density of states"\
                  ,xlim=(-1.5,1.5),ylim=(-8,7)\
                  )
        i.legend(loc='upper right')

    axs[0].set(ylabel="Energy relative to Fermi level/ eV")

    # changes which d-orbitals you want to emphasise, if any
    if not emph:
        pass
    else:

        u1[emph][0].set(linewidth=2.0)
        d1[emph][0].set(linewidth=2.0)
        u2[emph][0].set(linewidth=2.0)
        d2[emph][0].set(linewidth=2.0)

    return(fig,axs)

def DOSplotter_Natom(DOS,*args):
    """Plots (atom)'s d-orbital and total DOS from DOS_grab output

    Parameters:
    (array) : output from DOS_grab() function
    *args   : atoms for whcih you want to see the DOS 

    Returns:
    fig object - use fig.suptitle to add a title
    axs object - use axs[0 or 1].set to change things like xlim, ylim, titles
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from cycler import cycler

    fig, axs = plt.subplots(1,len(args), figsize=(16,14))
    fig.tight_layout(rect=[0, 0.03, 1, 0.92])
    plt.rc('font',size=16)

    labels = ["$s$","$p_x$","$p_y$","$p_z$","$xy$","$yz$","$z^2$",\
              "$xz$","$x^2-y^2$"] # orbital labels in VASP order
    default_cycler = (cycler(color=['g', 'g', 'r', 'r','k','k','m','m','b','b']))
    # this cycles through 5 different colours with the same colour for up and down spins
    plt.rc('lines', linewidth=1.0)
    plt.rc('axes', prop_cycle=default_cycler)

    uplots = []
    dplots = []
    for j, val in enumerate(args):
        for i in range(10,20):
            if i % 2 == 0: # checks if odd so that I can plot two-by-two the up and down spins
                uplots.append(axs[j].plot(DOS[val][i],DOS[0][0],label=labels[i//2-1], linewidth=0.75))
                dplots.append(axs[j].plot(DOS[val][i+1],DOS[0][0], linewidth=0.75))

    for i in axs:
        i.axhline(0, linestyle='-', color = 'k', linewidth = 0.5)
        i.axvline(0,linestyle='-', color = 'k', linewidth = 0.5)
        i.set(xlabel="Density of states"\
                  ,xlim=(-1.5,1.5),ylim=(-8,7)\
                  )
        i.legend(loc='upper right')

    axs[0].set(ylabel="Energy relative to Fermi level/ eV")

    return(fig,axs)

def DOSplotter_ax(DOS,axs,atom,emph=False):
    """Plots (atom)'s d-orbital and total DOS from DOS_grab output

    Parameters:
    DOS (array) : output from DOS_grab() function
    axs (plt axis obj)    : axs to which you want to plot
    atom (int)  : atom for which you want to see the DOS 
    emph (bool or int) : False for no highlighting, 0,1,2,3,4 for highlighting d-orbitals (vasp order)

    Returns:
    (list) uplots,dplots
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from cycler import cycler

    labels = ["$s$","$p_x$","$p_y$","$p_z$","$xy$","$yz$","$z^2$",\
              "$xz$","$x^2-y^2$"] # orbital labels in VASP order
    default_cycler = (cycler(color=['g', 'g', 'r', 'r','k','k','m','m','b','b']))
    # this cycles through 5 different colours with the same colour for up and down spins
    plt.rc('lines', linewidth=1.0)
    plt.rc('axes', prop_cycle=default_cycler)

    uplots = []
    dplots = []
    for i in range(10,20):
        if i % 2 == 0: # checks if odd so that I can plot two-by-two the up and down spins
            uplots.append(axs.plot(DOS[atom][i],DOS[0][0],label=labels[i//2-1], linewidth=0.75))
            dplots.append(axs.plot(DOS[atom][i+1],DOS[0][0], linewidth=0.75))

    axs.axhline(0, linestyle='-', color = 'k', linewidth = 0.5)
    axs.axvline(0,linestyle='-', color = 'k', linewidth = 0.5)
    axs.set(xlabel="Density of states"\
              ,xlim=(-1.5,1.5),ylim=(-8,7)\
              )
    axs.legend(loc='upper right')

    axs.set(ylabel="Energy relative to Fermi level/ eV")

    if not emph:
        pass
    else:
        uplots[emph][0].set(linewidth=2.0)
        dplots[emph][0].set(linewidth=2.0)

    return uplots,dplots

def DOSplotter_lobster(DOS,axs,atom,orbitals,emph=False):
    """Plots (atom)'s d-orbital and total DOS from DOS_grab output

    Parameters:
    DOS (array) : output from DOS_grab() function
    axs (plt axis obj)    : axs to which you want to plot
    atom (int)  : atom for which you want to see the DOS 
    orbitals (tuple) : orbital indices range to plot
    emph (bool or int) : False for no highlighting, 0,1,2,3,4 for highlighting d-orbitals (vasp order)

    Returns:
    (list) uplots,dplots
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from cycler import cycler

    labels = ["$s$","$p_x$","$p_y$","$p_z$","$xy$","$yz$","$z^2$",\
              "$xz$","$x^2-y^2$"] # orbital labels in VASP order
    default_cycler = (cycler(color=['g', 'g', 'r', 'r','k','k','m','m','b','b']))
    # this cycles through 5 different colours with the same colour for up and down spins
    plt.rc('lines', linewidth=1.0)
    plt.rc('axes', prop_cycle=default_cycler)

    uplots = []
    dplots = []
    for i in range(*orbitals):
        if i % 2 == 0: # checks if odd so that I can plot two-by-two the up and down spins
            uplots.append(axs.plot(DOS[atom][i],DOS[0][0],label=labels[i//2-1], linewidth=0.75))
            dplots.append(axs.plot(DOS[atom][i+1],DOS[0][0], linewidth=0.75))

    axs.axhline(0, linestyle='-', color = 'k', linewidth = 0.5)
    axs.axvline(0,linestyle='-', color = 'k', linewidth = 0.5)
    axs.set(xlabel="Density of states"\
              ,xlim=(-1.5,1.5),ylim=(-8,7)\
              )
    axs.legend(loc='upper right')

    axs.set(ylabel="Energy relative to Fermi level/ eV")

    if not emph:
        pass
    else:
        uplots[emph][0].set(linewidth=2.0)
        dplots[emph][0].set(linewidth=2.0)

    return uplots,dplots

def Eig_grab(homepath,fermi,nlin):
    """
    Extracts data from VASP EIGENVAL file to plot band-structures
    
    Parameters:
    (str) Path to EIGENVAL file
    (float) Fermi level, found by searching OUTCAR for fermi (2nd result) from good single point calc.
    (int) Number of lines calculated through k-space (inspect KPOINTS for this)
    
    Returns:
    (array,float) Eigenvalue data, format:
                  Eig[0][0] = kpoints [a,b,c]
                  Eig[i][0] = ith kpoint spin-up eigenvalues [band1,band2....nband]
                  Eig[i][1] = ith kpoint spin-down eigenvalues [band1,band2....nband]   
    (int) Number of kpoints
    (int) Number of points per line through k-space    
    (int) Number of lines thorugh k-space
    """
    import os
    import numpy as np
   
    with open(homepath,'r') as file:
        elines = file.readlines()        
    nkpt = int(elines[5].split()[1])
    nband = int(elines[5].split()[2])
    point_per_lin = int(nkpt/nlin)
    print('Number of kpoints: %i\nNumber of bands: %i\nPoints per line: %i\n' \
          % (nkpt,nband,point_per_lin))
    
    EIG = [[[],[]] for i in range(nband)]
    KPT = []
    for i in range(0, nkpt):
        KPT.append(list(map(float, elines[7+i*(nband+2)].split()[:3]))) # append coords of k-points as float
    for j in range(0, nband):    
        for i in range(0,nkpt):
            try:
                EIG[j][0].append(float(elines[8+i*(nband+2)+j].split()[1])-fermi) # spin-up eigenvalues
                EIG[j][1].append(float(elines[8+i*(nband+2)+j].split()[2])-fermi) # spin-down eigenvalues
            except:
                print('Problems with %i,%i' % j,i)

    EIG = np.array(EIG)
    KPT = np.array(KPT)
    return EIG,KPT,point_per_lin,nlin

def simple_bands(EIG,KPT,ppl,nlin):
    """Plots band structure given data from Eig_grab
    Parameters:
    EIG (array) : output from Eig_grab
    KPT (array) : kpoint info, output from Eig_grab
    ppl (int)   : number of kpoints per line (find in KPOINTS file
    nlin (int)  : number of lines through kspace (find in KPOINTS file))

    Returns:
    fig, axs objects
    """
    import matplotlib.pyplot as plt
    import numpy as np
    ubands = []
    dbands = []

    grid = plt.GridSpec(4,2,wspace=0.4, hspace=0.5)

    fig = plt.figure(figsize=(14,10))
    ax = [fig.add_subplot(grid[:2,:])]
    ax.append(fig.add_subplot(grid[2:,:]))
    plt.rc('font',size=12)

    ks = np.arange(0,np.shape(EIG)[2])

    for i in range(len(EIG)): # loop over bands
        ubands.append(ax[0].plot(ks,EIG[i,0],'b')) # plot band i, spin up
        dbands.append(ax[1].plot(ks,EIG[i,1],'g')) # plot band i, spin down)

    ax[0].axhline(0, linestyle='-', color = 'k', linewidth = 0.5)
    ax[1].axhline(0, linestyle='-', color = 'k', linewidth = 0.5)

    line_dividers = np.arange(0,1+np.shape(EIG)[2],ppl)
    special_ks    = np.append(KPT[::ppl],[KPT[-1]],axis=0).astype(int)

    for i in range(nlin+1):
        ax[0].axvline(line_dividers[i],linestyle='-', color = 'k', linewidth = 0.5)
        ax[1].axvline(line_dividers[i],linestyle='-', color = 'k', linewidth = 0.5)       

    ax[0].set(ylim=[-1.5,1.5],ylabel='Energy relative to Fermi level / eV',\
               xticks=line_dividers,xticklabels=special_ks,title='Majority Spin Channel')
    ax[1].set(ylim=[-1.5,1.5],xlabel='k-space',ylabel='Energy relative to Fermi level / eV',\
               xticks=line_dividers,xticklabels=special_ks,title='Minority Spin Channel')

    return fig,ax

class COHP_grab:
    """
    Splits up COHPCAR into 4 arrays: E, averaged, integrated and bond-decomposed

    parameters:
    homepath (str) : path to the COHPCAR file

    attributes:
    E : np array of the Energy scale
    av : 2D array of the averaged COHPS for up and down spin
    integ : 2D array for the integrated COHPs
    cohp : bond-decomposed array, format is cohp[spin][bond][value at E]
    cohp_axplot( ax to plot to, bond to plot, labels - show axis labels or not - False when plotting in loop)
    
    """
    def __init__(self,homepath):
        
        import os
        import numpy as np

        with open(homepath,'r') as file:
            lines = file.readlines()

            bonds = int(lines[1].split()[0]) - 1
            print('Number of bonds: %i' % bonds) # number of different atomic bonds in the calc
            value_number = int(lines[1].split()[2])
            print('Number of E values per centre: %i' % value_number) # number of different values per centre
            self.int_list = lines[3:3+bonds]
        
        self.E, self.av, self.integ = [], [[],[]], [[],[]]
        self.cohp = [[[] for i in range(bonds)],[[] for i in range(bonds)]]
        self.icohp = [[[] for i in range(bonds)],[[] for i in range(bonds)]]
        for i in range(0,value_number):
            self.E.append(float(lines[bonds+3+i].split()[0]))
            self.av[0].append(float(lines[bonds+3+i].split()[1]))
            self.av[1].append(float(lines[bonds+3+i].split()[2*bonds+3]))
            self.integ[0].append(float(lines[bonds+3+i].split()[2]))
            self.integ[1].append(float(lines[bonds+3+i].split()[2*bonds+4]))
            for j in range(0,2*bonds,2):
                self.cohp[0][j//2].append(float(lines[bonds+3+i].split()[(j+3)]))
                self.icohp[0][j//2].append(float(lines[bonds+3+i].split()[(j+4)]))
                self.cohp[1][j//2].append(float(lines[bonds+3+i].split()[2*bonds+5+j]))
                self.icohp[1][j//2].append(float(lines[bonds+3+i].split()[2*bonds+6+j]))
            
        self.E = np.array(self.E)
        self.av = np.array(self.av)
        self.integ = np.array(self.integ)
        self.cohp = -1*np.array(self.cohp)
        self.icohp = -1*np.array(self.icohp)

    def cohp_axplot(self,axs,bond,labels=True):
        """
        plots a cohp to an axis object of your choice (helper function for plotting on to subplots)
        Parameters:
        axs : axs slice sufficiently large to plot all the bonds you want to 
        """
        import matplotlib.pyplot as plt
        
        u = axs.plot(self.cohp[0][bond],self.E,'r-',label='Up spin')
        d = axs.plot(self.cohp[1][bond],self.E,'b-',label='Down spin')
        axs.axhline(0,color='k',linewidth=0.5)
        axs.axvline(0,color='k',linewidth=0.5)
        
        inter = self.int_list[bond].split('(')
        titel = inter[0].split(':')[1]+'\n'+inter[1][:-13]+u'\u212B'
        axs.set(title=titel)
        
        if labels:
            axs.set(xlabel = "pCOHP",xlim=(-1.5,1.5),ylim=(-10,10))
            axs.set(ylabel='Energy / eV')
        
        return(u,d)


def swap(plotter,xlim=None,ylim=None):
    """
    Swaps the axes and associated objects of a pymatgen Dosplotter() object
    """
    fig = plotter.get_plot().gcf()
    axs = fig.get_axes()[0]
    nxlabel = axs.get_ylabel()
    nylabel = axs.get_xlabel()
    line_list = axs.get_lines()
    for lines in line_list:
        try:
            iter(lines)
        except:
            lines = [lines]

        for line in lines:
            xdata, ydata = line.get_xdata(), line.get_ydata()
            line.set_xdata(ydata)
            line.set_ydata(xdata)
    if not xlim:
        nxlim = axs.get_ylim()
    else:
        nxlim = xlim
    if not ylim:
        nylim = axs.get_xlim()
    else:
        nylim = ylim
    axs.set(xlim = nxlim, ylim = nylim, xlabel = nxlabel, ylabel = nylabel)
    axs.autoscale_view()
    return fig
