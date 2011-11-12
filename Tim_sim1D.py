#Chad, this may be kind of a lot to throw at you, but I'm running a
#simulation in which I have many possible variations, and I want to
#flexibly specify them and maintain backward compatibility.
#You'll want to comment out line 287, which adds comments in the Finder
#specifying which parameters were used to run the model.

# Generate a piece of applescript that will take a POSIX file name
# and add spotlight tags to this. 
# We will run this script using the osascript command in Mac OS X.
# taken from http://snipplr.com/view.php?codeview&id=5015
def gen_add_comment_script(file,new_tags):
    temp_file = '"' + file + '"'
    osa_comment =  '"' + new_tags + '"'
    osascript = """osascript<<END
    tell application "Finder"
    set filePath to POSIX file %s
    set fileComment to %s
    set the_File to filePath as alias
    set comment of the_File to fileComment
    end tell 
END""" % (temp_file,  osa_comment) 
    return osascript

class nylon12():
    def get_k_type(self):
        return self.k_type
    def set_k_type(self,k_type):
        self.k_type = k_type
    k_type = property(get_k_type,set_k_type)

    def get_k_const(self):
        return self.k_const
    def set_k_const(self,k):
        self.k_const = k
    k_const = property(get_k_const,set_k_const)
    
    def get_c_type(self):
        return self.c_type
    def set_c_type(self,c_type):
        self.c_type = c_type
    c_type = property(get_c_type,set_c_type)

    def get_c_const(self):
        return self.c_const
    def set_c_const(self,cp):
        self.c_const = c_const
    c_const = property(get_c_const,set_c_const)
    
    def get_rho(self):
        return self.rho
    def set_rho(self,rho):
        self.rho = rho
    rho = property(get_rho,set_rho)
    
    def k(self,T=350.):
        from numpy import array,polyval
        if 'con' in self.k_type:
            return  self.k_const
        elif 'lin' in self.k_type:
            p=array([ 0.00018587,  0.00110113])
            return polyval(p,T)#1.94e-4*T - 8.68e-4 #W/mK

    def c(self,T=350.):
        from numpy import array,polyval
        if 'con' in self.c_type:
            return  self.c_const
        elif 'lin' in self.c_type:
            p=[7.9806,-685.27]
            return polyval(p,T)
    
    def alpha(self,T=350.):
        return self.k(T)/self.rho/self.c(T) #m2/s - thermal diffusivity

    def __init__(self,k_type='linear',c_type='linear'):
        self.k_const = 0.072852 #W/mk
        self.k_type = k_type
        self.rho = 490. #kg/m3
        self.c_const = 1800. #kJ/kg/K
        self.c_type = c_type
            
def Fo_T(alpha,dt,dx):
    return alpha*dt/dx**2 #Fourier number (dimensionless time)

def sim1D(**kwargs):
    """
    mse = sim1D(**kwargs)
    1-dimensional simulation of proportional controlled heat flux into a semi-infinite medium.
    Returns the mean sqared error between the surface solution and 'data' if supplied.
    Available keyword arguments are shown with default values

    ---solution parameters---
    dt = 1.0                   [s] time increment
    time = arange(0.,1800.,dt) [s] time domain for comparison, if comparison data is available
    tmax = 100.                [s] if no time domain is given, this is the maximum time of the simulation
    xmax = 0.2                 [m] maximum depth at 'infinity'
    dx = 1.016e-4              [m] spatial domain differential element size .005
    x = arange(0.,xmax,dx)     [m]
    data = []                  [K] temperature data for comparison, if comparison data is available
                                   if temperature data is available, it will be plotted at the end.
    k_type = 'linear'          [W/mK] linear model for thermal conductivity
             'constant'        [W/mK] 0.072852 unless specific with a value for k

    ---material parameters---
    k = 0.072852               [W/mk] thermal conductivity, if k_type = 'constant'
    cp = 1800                  [kJ/kgK] specific heat
    rho = 490                  [kg/m3] density
    Kp = 0.04                  [W/m2/K] 0.2, proportional control coefficient

    ---problem parameters---
    mtas = 13                  mean-time-average samples to use in calculating mean temperature
    T_initial = 330.           [K] initial temperature
    T_set = 413.               [K] set-point temperature
    T_offset = 3.0             [K] offset to add to set-point to remove steady-state error

    ---output parameters---
    filename='last sim'      filename for saving

"""
    import matplotlib.pyplot as plt
    from   matplotlib import rc
    import numpy as np
    import os
    import progressbar as pb

    #Settings to make the plots appropriate for inclusion in TeX generated publications
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text',usetex=True)
    FONTSIZE = 10
    FIGSIZE = (3.5,3.5)
    FIGDIM = ([0.15,0.1,0.8,0.85])

#Proportional control coefficient
    if 'Kp' in kwargs:
        Kp = kwargs['Kp']
    else:
        Kp = .04

#number of time samples
    if 'mtas' in kwargs:
        moving_time_average_samples = mtas
    else:
        moving_time_average_samples = 15

#surface flux
    if 'qs_nom' in kwargs:
        qs_nom = kwargs['qs_nom']
    else:
        qs_nom = 600. #500. #585. #W

#material properties
    if 'k_type' in kwargs:
        m=nylon12(kwargs['k_type'])      #instantiate m - material
        if 'const' in kwargs['k_type']:
            if 'k' in kwargs:
                m.k_const = kwargs['k']
                print 'k found\n'
    else:
        m = nylon12('linear')
        print 'using default linear thermal conductivity.\n'
        
#specific heat
    if 'c_type' in kwargs:
        m.c_type = kwargs['c_type']
        if 'const' in kwargs['c_type']:
            if 'c' in kwargs:
                m.c_const = kwargs['c']
                print 'constant c found'
    else:
        print 'using default linear specific heat'
        
#density
    if 'rho' in kwargs:
        m.rho = kwargs['rho']

#spatial domain
    if 'xmax' in kwargs:
        xmax = kwargs['xmax']
    else:
        xmax = 0.02 #[m] depth of powder to consider
    if 'dx' in kwargs:
        dx = kwargs['dx']
    else:
        dx = 1.016e-4
    if 'x' in kwargs:
        x = np.asarray(kwargs['x'])
    else:
        x = np.arange(0,xmax,dx)

#Temperatures
    if 'T_initial' in kwargs:
        T_initial = kwargs['T_initial']
    else:
        T_initial = 300
        
    if 'T_offset' in kwargs:
        T_offset = kwargs['T_offset']
    else:
        T_offset = 3
        
    if 'T_set' in kwargs:
        T_set = kwargs['T_set']
    else:
        T_set = 470

#time domain
    if 'time' in kwargs: #set up time variable
        time = kwargs['time']
        dt = time[1] - time[0]
        if 'data' in kwargs:
            data = kwargs['data']
            Compare = True
        else:
            Compare = False
    else: #use default
        dt = dx**2/(5*m.alpha(T_set)) #stability criterion Fo<=1/2
        if 'tmax' in kwargs:
            tmax = float(kwargs['tmax'])
        else:
            tmax = 100.
        time = np.arange(0.,tmax+dt,dt)
        Compare = False
    tmax = max(time)
    num_time_steps = len(time)

#initialize the working variables
    T     = np.ones((num_time_steps,len(x)))*T_initial
    qs    = np.zeros(num_time_steps)
    err   = np.zeros(num_time_steps)
    u     = np.zeros(num_time_steps)

#loop through the time and space domains
    inf = len(x)-1
    print "Solving ...\n"
    pbar=pb.ProgressBar().start()
    for i in range(1,num_time_steps): #time step
        dt = time[i] - time[i-1]
        #constant flux boundary condition
        err[i] = T_set + T_offset - np.mean(T[range(max(0,i-moving_time_average_samples),i),0])
        u[i] = err[i] * Kp
        qs[i] = max(min(1.,u[i]) * qs_nom,-10)
        T[i,0] = 2*Fo_T(m.alpha(T[i-1,0]),dt,dx)*(T[i-1,1] + qs[i]*dx/m.k(T[i-1,1])) + (1 - 2*Fo_T(m.alpha(T[i-1,0]),dt,dx)) * T[i-1,0]

        #adiabatic far wall boundary condition
        T[i,inf] = 2*Fo_T(m.alpha(T[i-1,inf-1]),dt,dx) * T[i-1,inf-1] + (1 - 2*Fo_T(m.alpha(T[i-1,inf]),dt,dx)) * T[i-1,inf]

        #internal nodes heat equation
        for j in range(1,len(x)-1):
            T[i,j] = Fo_T(m.alpha(T[i-1,j]),dt,dx) * (T[i-1,j-1] + T[i-1,j+1]) + (1 - 2*Fo_T(m.alpha(T[i-1,j]),dt,dx)) * T[i-1,j]
        pbar.update(100.*float(i)/float(num_time_steps))
    pbar.finish()

#plot the results
    print "Plotting ...\n"
    fig = plt.figure(1,figsize=FIGSIZE)
    ax = fig.add_axes(FIGDIM)
    plotlabel = 'dx=%1.2e, Fo=%1.2e' %(dx,Fo_T(m.alpha(T_set),dt,dx))
    line = ax.plot(time,T[:,0],label=plotlabel)
    if(Compare):
        line2 = ax.plot(time,data,label='Reference')
    xtext = ax.set_xlabel('Time (s)',fontsize=FONTSIZE,family='sans-serif')
    ytext = ax.set_ylabel('Surface Temperature (K)',fontsize=FONTSIZE,family='sans-serif')
    for label in ax.get_xticklabels():
        label.set_family('sans-serif')

    if 'filename' in kwargs:
        filename = kwargs['filename']
    else:
        filename = 'last_sim'

    np.savez(filename,T=T,time=time,qs=qs)

    figfilename = filename+'.pdf'
    plt.savefig(figfilename,format='pdf')

    comment_info = "qs_nom = %.0f\nT_set = %1.1f\nKp = %1.3f\nT_initial = %1.3f\nT_set = %1.1f\nT_offset = %1.1f\ndx = %1.3e\ndt=%1.3e" % (qs_nom,
    T_set,
    Kp,
    T_initial,
    T_set,
    T_offset,
    dx,
    dt)
    
    os.system(gen_add_comment_script(figfilename,comment_info))
    try:
        rmse = np.sqrt( np.mean( (T[:,0]-data)**2 ) )
        return rmse
    except:
        return -1.
    
