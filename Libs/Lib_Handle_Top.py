import numpy as np
from numba import jit

def Calc_MatricesTop_3D(nx,nz,nu,om):
    jqmax = int(nx/2)
    MatricesTop = np.zeros((2*jqmax+1,2*jqmax+1,3,3),dtype = np.complex128)
    for jx in range(-jqmax,jqmax+1):
        for jz in range(-jqmax,jqmax+1):
            qx = 2.*np.pi*jx/(nx-1)
            qz = 2.*np.pi*jz/(nz-1)
            phi = np.arctan2(qz,qx)
            q = qx*np.cos(phi) + qz*np.sin(phi)
            Z = (1j*om*nu)**0.5 #rho ==1
            k0 =(1j*om/nu)**0.5; xi = (q**2 + k0**2)**0.5/k0
            R = np.array([[np.cos(phi) ,0,np.sin(phi)],\
                          [0           ,1,0          ],\
                          [-np.sin(phi),0,np.cos(phi)]])
            R_inv = np.linalg.inv(R)    
            sq = np.sign(q)
            Fq = np.zeros((3,3),dtype = np.complex128)
            if jx == 0 and jz == 0: 
                Fq[0,0] = 1./(Z+1./3.)
                Fq[2,2] = 1./(Z+1./3.)
            else:     
                E = ((q**2 - 2*k0*q*sq*xi + k0**2*xi**2)*Z**2)/\
                   ((q - k0*sq*xi)*(-q + k0*sq*xi)) +\
                   (1 - (k0**2*xi*Z)/(q**2 - k0*q*sq*xi))*\
                     (1./3 + (sq*(-q**2 + k0**2*xi**2)*Z)/(k0*(-q + k0*sq*xi)))
                Fq[0,0] = 1./E*\
                    (1-(k0**2*xi*Z)/(q**2-k0*q*sq*xi))
                Fq[0,1] = 1./E*\
                   1j*(q**2 - 2*k0*q*sq*xi + k0**2*xi**2)*Z /(k0*(-q + k0*sq*xi))
                Fq[0,2] = 0.
                Fq[1,0] = 1./E*\
                    ( 1j*k0*Z)/(q-k0*sq*xi)
                Fq[1,1] = 1./E*\
                    (1./3. + (sq*(-q**2 + k0**2*xi**2)*Z)/(k0*(-q + k0*sq*xi)))
                Fq[1,2] = 0.
                Fq[2,0] = 0
                Fq[2,1] = 0
                Fq[2,2] = 1./(1./3. + xi*Z)
            MatricesTop[jx,jz] = np.dot(R_inv,np.dot(Fq,R))
    return MatricesTop 

@jit(nopython = True)
def Calc_uT_Fourier_3D(h,nx,nz,MatricesTop):
    uxT = np.zeros((nx,nz),dtype = np.complex128)
    uyT = np.zeros((nx,nz),dtype = np.complex128)
    uzT = np.zeros((nx,nz),dtype = np.complex128)
    jqmax = int(nx/2)
    for jx in range(-jqmax,jqmax+1):
        for jz in range(-jqmax,jqmax+1):
            qx             = 2.*np.pi*jx/nx
            if nz == 1: qz = 0 
            else      : qz = 2.*np.pi*jz/nz
            Hxy = 0; Hyy = 0; Hzy = 0
            for x in range(nx):
                for z in range(nz):
                    Hxy  += 2 * (h[x,-1,z,7 ] - h[x,-1,z,14]) *\
                        np.exp(-1j*qx*x) * np.exp(-1j*qz*z)/(nx*nz)  
                    Hyy  += 2 * (h[x,-1,z,7 ] + h[x,-1,z,14] + \
                                  h[x,-1,z,11] + h[x,-1,z,17] + \
                                  h[x,-1,z,3 ]) *\
                        np.exp(-1j*qx*x) * np.exp(-1j*qz*z)/(nx *nz)
                    Hzy  += 2 * (h[x,-1,z,11] - h[x,-1,z,17]) * \
                        np.exp(-1j*qx*x) * np.exp(-1j*qz*z)/(nx *nz)
            Hq_vec = np.array([Hxy,Hyy,Hzy],dtype = np.complex128)
            uq_vec = np.dot(MatricesTop[jx,jz],Hq_vec)
            for x in range(nx):
                for z in range(nz):
                    uxT[x,z] += uq_vec[0] * np.exp(1j*qx*x)*np.exp(1j*qz*z)
                    uyT[x,z] += uq_vec[1] * np.exp(1j*qx*x)*np.exp(1j*qz*z)
                    uzT[x,z] += uq_vec[2] * np.exp(1j*qx*x)*np.exp(1j*qz*z)
    return uxT,uyT,uzT

@jit(nopython = True)
def Calc_uT_Local_3D(h,nx,nz,om,tauInvs):  # for comppleteness, never called
    uxT = np.zeros((nx,nz),dtype = np.complex128)
    uyT = np.zeros((nx,nz),dtype = np.complex128)
    uzT = np.zeros((nx,nz),dtype = np.complex128)
    nus = np.zeros((nx,nz),dtype = np.complex128)
    ZLiqs = np.zeros((nx,nz),dtype = np.complex128)
    cS = 3**-0.5; ZComp = (cS)**0.5; ZComp = np.complex128(ZComp)
    for x in range(nx):
        for z in range(nz):
            nus[  x,z] = (1./tauInvs[x,-1,z]-0.5)/3.
            ZLiqs[x,z] = (1j*om*nus[x,z])**0.5 # rho = 1
            uxT[x,z] = 2*(h[x,-1,z,7 ]-h[x,-1,z,14]             )/(ZLiqs[x,z] + 1./3.)
            uzT[x,z] = 2*(h[x,-1,z,11]-h[x,-1,z,17]             )/(ZLiqs[x,z] + 1./3.)
            uyT[x,z] = 2*(h[x,-1,z,7 ]+h[x,-1,z,14]+\
                          h[x,-1,z,11]+h[x,-1,z,17]+h[x,-1,z,3 ])/(ZComp + 1.)      
    return uxT,uyT,uzT

def Calc_MatricesTop_2D(nx,nu,om):
    jqmax_Top = int(nx/2)
    MatricesTop = np.zeros((2*jqmax_Top+1,2,2),dtype = np.complex128)
    for jx in range(-jqmax_Top ,jqmax_Top +1):
        qx = 2.*np.pi*jx/(nx-1)
        q = qx
        Z = (1j*om*nu)**0.5 #rho ==1
        k0 =(1j*om/nu)**0.5; xi = (q**2 + k0**2)**0.5/k0
        sq = np.sign(q)
        Fq = np.zeros((2,2),dtype = np.complex128)
        if jx == 0 : 
            Fq[0,0] = 1./(Z+1./3.)
            Fq[2,2] = 0
        else:     
            E = ((q**2 - 2*k0*q*sq*xi + k0**2*xi**2)*Z**2)/\
               ((q - k0*sq*xi)*(-q + k0*sq*xi)) +\
               (2 - (k0**2*xi*Z)/(q**2 - k0*q*sq*xi))*\
                 (1./3 + (sq*(-q**2 + k0**2*xi**2)*Z)/(k0*(-q + k0*sq*xi)))
            Fq[0,0] = 1./E*\
                (2-(k0**2*xi*Z)/(q**2-k0*q*sq*xi))
            Fq[0,1] = 1./E*\
               1j*(q**2 - 2*k0*q*sq*xi + k0**2*xi**2)*Z /(k0*(-q + k0*sq*xi))
            Fq[1,0] = 1./E*\
                ( 1j*k0*Z)/(q-k0*sq*xi)
            Fq[1,1] = 1./E*\
                (1./3. + (sq*(-q**2 + k0**2*xi**2)*Z)/(k0*(-q + k0*sq*xi)))
        MatricesTop[jx] = Fq
    return MatricesTop 

@jit(nopython = True)
def Calc_uT_Fourier_2D(h,nx,MatricesTop):
    uxT = np.zeros((nx),dtype = np.complex128)
    uyT = np.zeros((nx),dtype = np.complex128)
    jqmax_Top = int(nx/2)
    for jx in range(-jqmax_Top,jqmax_Top+1):
        qx = 2.*np.pi*jx/nx
        Hxy = 0; Hyy = 0
        for x in range(nx):
            Hxy+=2*(h[x,-1,5]-h[x,-1,6])          *np.exp(-1j*qx*x)/nx
            Hyy+=2*(h[x,-1,5]+h[x,-1,6]+h[x,-1,7])*np.exp(-1j*qx*x)/nx
        Hq_vec = np.array([Hxy,Hyy],dtype = np.complex128)
        uq_vec = np.dot(MatricesTop[jx],Hq_vec)
        for x in range(nx):
            uxT[x] += uq_vec[0] * np.exp(1j*qx*x)
            uyT[x] += uq_vec[1] * np.exp(1j*qx*x)
    return uxT,uyT

def Calc_uT_Local_2D(h,nx,om,tauInvs): # for comppleteness, never called
    uxT = np.zeros((nx),dtype = np.complex128)
    uyT = np.zeros((nx),dtype = np.complex128)
    nus   = np.zeros((nx),dtype = np.complex128)
    ZLiqs = np.zeros((nx),dtype = np.complex128)
    cS = 3**-0.5; ZComp = (cS)**0.5; ZComp = np.complex128(ZComp)
    for x in range(nx):
        nus[  x] = (1./tauInvs[x,-1]-0.5)/3.
        ZLiqs[x] = (1j*om*nus[x])**0.5 # rho = 1
        uxT[x] = 2*(h[x,-1,5 ]-h[x,-1,6 ]          )/(ZLiqs[x] + 1./3.)   
        uyT[x] = 2*(h[x,-1,5 ]+h[x,-1,6 ]+h[x,-1,2])/(ZComp + 2.)
    return uxT,uyT



