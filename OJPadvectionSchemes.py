# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py 

# The numpy package for numerical functions and pi
import numpy as np

def FTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTCS, Courant number c"
    "for nt time-steps"
    
    nx = len(phiOld)
    ## new time-step array for phi
    phi = phiOld.copy()

    ## FTCS for each time-step
    for it in range(nt):
        ## Loop through all space using remainder after division (%)
        ## to cope with periodic boundary conditions
        
        for j in range(nx):
            phi[j] = phiOld[j] - 0.5*c*\
                     (phiOld[(j+1)%nx] - phiOld[(j-1)%nx])
                     
            
        
        ## update arrays for next time-step
        phiOld = phi.copy()
        

    return phi

    
def FTBS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTBS, Courant number c"
    "for nt time-steps"
    nx=len(phiOld)
    
    phi = phiOld.copy()
    ## FTBS for each time step
    for it in range(nt):
        for j in range(nx):
            phi[j] = phiOld[j] - c*(phiOld[j] - phiOld[(j-1)%nx])
            
            
        phiOld = phi.copy()
            
    return phi

    


        

def CTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using CTCS, Courant number c"
    "for nt time-steps, using FTCS to create second time step"
    nx=len(phiOld)
    phi = np.zeros((nt,nx))  ##making matrix for n time steps and j spacial 
                             ##steps
    
    phi[0,:] = phiOld.copy()  ##initial conditions for time step t=0
    phi[1,:] = FTCS(phi[0,:], c, 1)  ##initial conditions for t=1 using FTCS  

    ##CTCS for each time step
    for n in range(1, nt-1):
        for j in range(nx):
            phi[n+1,j] = phi[n-1,j] - c*(phi[n,(j+1)%nx] - phi[n, (j-1)%nx])
            
    return phi
    



def LW(phiOld, c, nt):
    "Linear advection of profile in phiOld using Lax-Wendroff, Courant number c"
    "for nt time-steps"
    nx=len(phiOld)
    ##initialsing phi and phi at plus half and minus half  for spacial steps
    ##calculation
    phi = phiOld.copy()
    phijplus = phiOld.copy()
    phijminus = phiOld.copy()
   
    for n in range(nt):
        
        ##loop for each spacial step, doing half plus and half minus steps
        ##then computing for next time step
        for j in range(nx):
            phijplus[j] = 0.5*(1+c)*phiOld[j] +0.5*(1-c)*phiOld[(j+1)%nx]
            phijminus[j]= 0.5*(1+c)*phiOld[(j-1)%nx] +0.5*(1-c)*phiOld[j]
            phi[j]=phiOld[j]-c*(phijplus[j]-phijminus[j])
              
       
        phiOld = phi.copy()
    
    return phi
    
