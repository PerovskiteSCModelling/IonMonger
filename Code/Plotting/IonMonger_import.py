'''
Functions to load a solution file from IonMonger and extract the data for
further analysis and plotting. The package that transfers the data from MATLAB,
`scipy.io`, formats structures as dictionaries in Python. However, during the
conversion, some variables are placed in a seemingly arbitrary number of nested
containers. For this reason, when extracting variables, it may be necessary to
add additional indexing. For example, `params['dE'][0][0][0]` to obtain the
ETL doping density. Note that impedance spectroscopy simulations must be
unpacked using the `unpack_IS` function, instead of `unpack_sol`. Ensure that
Python's directory is set to the folder containing the solution file.
'''
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import os

#%% Define the unpacking functions

def unpack_sol(filename):
    '''
    A function to load a MATLAB file containing a solution structure produced
    by IonMonger into the Python environment and unpack the key variables for
    analysis/plotting.
    
    Inputs
    ------
    filename : string
        The name of the MATLAB file containing the solution structure. Must 
        contain the entire path from the working directory to the file.

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
    
    # Load the solution file
    mat = scipy.io.loadmat(filename)

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

    # ============= Analysis =============
    '''
    Here can be written any analysis to be performed on the data, including
    plotting. Examples for plotting JV curves and carrier densities are 
    included.
    '''
    plt.figure(figsize=(8, 5), dpi=200)
    plt.plot(V[100:200],J[100:200],'-b',label='reverse')
    plt.plot(V[200:300],J[200:300],'--b',label='forward')
    plt.ylim(-2,25)
    plt.xlim(0,1.2)
    plt.xlabel(r'applied voltage (V)')
    plt.ylabel(r'current density (mAcm$^{-2}$)')
    plt.legend()
    plt.show()
    
    plt.figure(figsize=(5, 3), dpi=200)
    plt.semilogy(x,P[150,:],'-g',label='anion vacancies')
    plt.semilogy(x,n[150,:],'-b',label='electrons')
    plt.semilogy(x,p[150,:],'-r',label='holes')
    plt.semilogy(xE,nE[150,:],'-b',label='')
    plt.semilogy(xH,pH[150,:],'-r',label='')
    plt.xlabel(r'$x$ (nm)')
    plt.ylabel(r'number density (m$^{-3}$)')
    plt.legend()
    plt.show()
    
    
def unpack_IS(filename):
    '''
    A function to load a MATLAB file containing a solution structure produced
    by IonMonger for an impedance spectroscopy simulation into the Python
    environment and unpack the key variables for analysis/plotting. This
    function is compatible only with the reduced output solution of IonMonger.
    
    Inputs
    ------
    filename : string
        The name of the MATLAB file containing the solution structure. Must 
        contain the entire path from the working directory to the file.

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
            
    J :     vector
            Current density at DC voltage (mAcm-2)
    Jl :    vector
            Current density at DC voltage lost to ETL-perovskite surface
            recombination (mAcm-2)
    Jr :    vector
            Current density at DC voltage lost to HTL-perovskite surface
            recombination (mAcm-2)
            
    f :     vector
            Sample frequencies (Hz)
    X :     vector
            Imaginary part of the impedance at each sample frequency (Ohm cm2)
    R :     vector
            Real part of the impedance at each sample frequency (Ohm cm2)
                
    '''
    
    # Load the solution file
    mat = scipy.io.loadmat(filename)

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
    
    # Steady state distributions at DC voltage
    dstrbns = mat['sol']['dstrbns'][0][0]
    
    P = dstrbns['P'][0][0][0]
    phi = dstrbns['phi'][0][0][0]
    n = dstrbns['n'][0][0][0]
    p = dstrbns['p'][0][0][0]
    phiE = dstrbns['phiE'][0][0][0]
    nE = dstrbns['nE'][0][0][0]
    phiH = dstrbns['phiH'][0][0][0]
    pH = dstrbns['pH'][0][0][0]
    
    # Currents
    J = mat['sol']['J'][0][0]
    J = np.resize(J,[J.shape[0],]) # reshape into vector
    Jl = mat['sol']['Jl'][0][0]
    Jl = np.resize(Jl,[Jl.shape[0],]) # reshape into vector
    
    # Impedance variables
    f = mat['sol']['freqs'][0][0][0]
    R = mat['sol']['R'][0][0]
    R = np.resize(R,[R.shape[0],]) # reshape into vector
    X = mat['sol']['X'][0][0]
    X = np.resize(X,[X.shape[0],]) # reshape into vector
    
    # ============= Analysis =============
    '''
    Here can be written any analysis to be performed on the data, including
    plotting. Examples for plotting Nyquist and Bode plots are included.
    '''
    plt.figure(figsize=(5, 3), dpi=200)
    plt.plot(R,X,'-or',markersize=2,linewidth=1)
    plt.gca().invert_yaxis()
    plt.xlabel(r'$R$ ($\Omega cm^2$)')
    plt.ylabel(r'$X$ ($\Omega cm^2$)')
    plt.gca().set_aspect('equal')
    plt.show()
    
    fig, ax = plt.subplots(2, 1,figsize=(5, 3), dpi=200)
    ax[0].semilogx(f,X,'-g')
    ax[0].invert_yaxis()
    ax[0].set_xticklabels([])
    ax[1].semilogx(f,R,'-g')
    plt.xlabel(r'frequency (Hz)')
    ax[0].set_ylabel(r'$X$ ($\Omega cm^2$)')
    ax[1].set_ylabel(r'$R$ ($\Omega cm^2$)')
    plt.show()
    
# ======== Examples to execute the above functions ========
# unpack_sol('Data/simulation.mat')
unpack_IS('Data/simulation.mat')
