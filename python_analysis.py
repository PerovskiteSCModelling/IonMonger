'''
A script to load a solution file from IonMonger and extract the data for
further analysis and plotting. The package that transfers the data from MATLAB,
`scipy.io`, formats structures as dictionaries in Python. However, during the
conversion, some variables are placed in a seemingly arbitrary number of nested
containers. For this reason, when extracting variables, it may be necessary to
add additional indexing. For example, `params['dE'][0][0][0]` to obtain the
ETL doping density.
'''
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

# Load the solution file
mat = scipy.io.loadmat('Data/simulation.mat')
'''
Parameters extracted from solution file:
----------------------------------------
    x:  vector
        the gridpoints in the perovskite (nm)
    xE: vector
        the gridpoints in the ETL (nm)
    xH: vector
        the gridpoints in the HTL (nm)
        
        
    params :    array
                Contains all of the fields of the params structure.
    
    P :     array
            Anion vacancy density in the perovskite (m-3). Arranged as P[t,x].
    phi :   array
            Electric potential in the perovskite (V). Arranged as phi[t,x].
    n :     array
            Electron density in the perovskite (m-3). Arranged as n[t,x].
    p :     array
            Hole density in the perovskite (m-3). Arranged as p[t,x].
    phiE :  array
            Electric potential in the ETL (V). Arranged as phiE[t,x].
    nE :    array
            Electron density in the ETL (m-3). Arranged as nE[t,x].
    phiH :  array
            Electric potential in the HTL (V). Arranged as phiH[t,x].
    pH :    array
            Hole density in the HTL (m-3). Arranged as pH[t,x].
            
    time :  vector
            Time grid (s)
    V :     vector
            Applied voltage (V)
    J :     vector
            Current density (mAcm-2)
    Jl :    vector
            Current density lost to ETL-perovskite surface recombination (mAcm-2)
    Jr :    vector
            Current density lost to HTL-perovskite surface recombination (mAcm-2)
            
'''

# ============= Extract data =============

# Vectors
vectors = mat['sol']['vectors'][0][0][0]

x = vectors['x'][0] # gridpoints in perovskite layer
x = np.resize(x,[x.shape[0],]) # reshape into vector
xE = vectors['xE'][0] # gridpoints in ETL
xE = np.resize(xE,[xE.shape[0],]) # reshape into vector
xH = vectors['xH'][0] # gridpoints in HTL
xH = np.resize(xH,[xH.shape[0],]) # reshape into vector

# Parameters
params = mat['sol']['params'][0][0][0]

# Distributions
dstrbns = mat['sol']['dstrbns'][0][0]

P = dstrbns['P'][0][0]
phi = dstrbns['phi'][0][0]
n = dstrbns['n'][0][0]
p = dstrbns['p'][0][0]
phiE = dstrbns['phiE'][0][0]
nE = dstrbns['nE'][0][0]
phiH = dstrbns['phiH'][0][0]
pH = dstrbns['pH'][0][0]

# Time
time = mat['sol']['time'][0][0][0]

# Applied voltage
V = mat['sol']['V'][0][0]
V = np.resize(V,[V.shape[0],]) # reshape into vector

# Currents
J = mat['sol']['J'][0][0]
J = np.resize(J,[J.shape[0],]) # reshape into vector
Jl = mat['sol']['Jl'][0][0]
Jl = np.resize(Jl,[Jl.shape[0],]) # reshape into vector

# ============= Plotting =============

'''
An example for plotting electric potential during the first section of the
protocol.
'''
ind = 180 # select the timestep(s) to plot
plt.plot(x, np.transpose(phi[ind,:]),'r')
plt.plot(xE, np.transpose(phiE[ind,:]),'r')
plt.plot(xH, np.transpose(phiH[ind,:]),'r')
plt.xlabel(r'$x$ (nm)')
plt.ylabel(r'electric potential, $\phi$ (V)')
plt.title(f'electric potential at t = {round(time[ind],2)}s')
plt.grid()
plt.show()

'''
An example for plotting a J-V curve.
'''
plt.plot(V[100:200],J[100:200],'-b',label='reverse')
plt.plot(V[200:300],J[200:300],'--b',label='forward')
plt.ylim(0,10)
plt.xlim(0,1.2)
plt.xlabel(r'applied voltage (V)')
plt.ylabel(r'current density (mAcm$^{-2}$)')
plt.legend()
plt.show()

