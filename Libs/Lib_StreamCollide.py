import numpy as np
from numba import jit
from Libs import Lib_Handle_Top as Handle_Top
from Libs import Lib_Relax      as Relax

@jit(nopython = True)
def FreqDLBMStep_SoftPt_3D(h,nx,ny,nz,\
        nd,cxs,cys,czs,wi,ibars,i_ups,i_notups,i_downs,i_notdowns,MatricesTop,\
        tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym):
    h_str      = np.zeros((nx,ny,nz,nd),dtype = np.complex128)
    Fx_on_Wall = np.zeros((nx,nz)      ,dtype = np.complex128)
    uxT,uyT,uzT = Handle_Top.Calc_uT_Fourier_3D(h,nx,nz,MatricesTop)
    for x in range(nx):
        for z in range(nz):
            for y in range(1,ny-1):
                for i in range(nd):
                    xmd = (nx+x-cxs[i])%nx 
                    ymd =     y-cys[i]
                    zmd = (nz+z-czs[i])%nz 
                    h_str[x,y,z,i] = h[xmd,ymd,zmd,i]
            for y in range(1):
                for i in i_notups: 
                    xmd = (nx+x-cxs[i])%nx 
                    ymd =     y-cys[i]
                    zmd = (nz+z-czs[i])%nz 
                    h_str[x,y,z,i] = h[xmd,ymd,zmd,i]
                for i in i_ups:  
                    ibar = ibars[i]
                    h_str[x,y,z,i] = h[x,y,z,ibar] - 6*wi[ibar]*cxs[ibar] 
                    Fx_on_Wall[x,z] += (h_str[x,y,z,i]+h[x,y,z,ibar])*cxs[i]                  
            for y in range(ny-1,ny):
                for i in i_notdowns: 
                    xmd = (nx+x-cxs[i])%nx 
                    ymd =     y-cys[i]
                    zmd = (nz+z-czs[i])%nz 
                    h_str[x,y,z,i] = h[xmd,ymd,zmd,i]
                for i in i_downs:  
                    ibar = ibars[i]
                    h_str[x,y,z,i] = h[x,y,z,ibar] \
                        - 6 * wi[ibar]*cxs[ibar] * uxT[x,z]\
                        - 6 * wi[ibar]*cys[ibar] * uyT[x,z]\
                        - 6 * wi[ibar]*czs[ibar] * uzT[x,z]
    for x in range(nx):
        for y in range(ny): 
            for z in range(nz): 
                h[x,y,z] = Relax.Relax_3D(h_str[x,y,z],cxs,cys,czs,wi,ibars,\
                tauInvs[x,y,z],tauInvs_Asym[x,y,z],\
                one_m_tauInvs_m_Iom[x,y,z],one_m_tauInvs_m_Iom_Asym[x,y,z])
    return h,Fx_on_Wall

@jit(nopython = True)
def FreqDLBMStep_OscBnd_3D(h,nx,ny,nz,nd,cxs,cys,czs,wi,ibars,MatricesTop,\
        nLstot,OutsideLBMDomains,i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
        tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym):
    h_str      = np.zeros((nx,ny,nz,nd),dtype = np.complex128)
    Fx_on_Wall = np.zeros((nx,nz)      ,dtype = np.complex128)
    FxLs       = np.zeros(nLstot       ,dtype = np.complex128) 
    FyLs       = np.zeros(nLstot       ,dtype = np.complex128) 
    FzLs       = np.zeros(nLstot       ,dtype = np.complex128) 
    uxT,uyT,uzT = Handle_Top.Calc_uT_Fourier_3D(h,nx,nz,MatricesTop)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if not OutsideLBMDomains[x,y,z]: 
                    for i in range(nd):
                        if i_BCs[x,y,z,i] == 1: #bulk
                            xmd = (nx+x-cxs[i])%nx
                            ymd =     y-cys[i]
                            zmd = (nz+z-czs[i])%nz
                            h_str[x,y,z,i] = h[xmd,ymd,zmd,i]
                        if i_BCs[x,y,z,i] == 2: #bottom
                            ibar = ibars[i] # i and ibar reversed compared to Krüger et al.
                            h_str[x,y,z,i] = h[x,y,z,ibar]-6*wi[ibar]*cxs[ibar]
                            Fx_on_Wall[x,z]+=(h_str[x,y,z,i]+h[x,y,z,ibar])*cxs[i]                  
                        if i_BCs[x,y,z,i] == 3: #top
                            ibar = ibars[i]
                            h_str[x,y,z,i] = h[x,y,z,ibar]\
                               -6*wi[ibar]*cxs[ibar]*uxT[x,z]\
                               -6*wi[ibar]*cys[ibar]*uyT[x,z]\
                               -6*wi[ibar]*czs[ibar]*uzT[x,z]
                        if i_BCs[x,y,z,i] == 4: #Surface, no small gap
                            ibar = ibars[i]  
                            iL = PoiLs[x,y,z,i]
                            q = qs[iL]
                            if q <= 0.5: 
                                xpd = (nx+x+cxs[i])%nx
                                ypd =     y+cys[i]
                                zpd = (nz+z+czs[i])%nz
                                h_str[x,y,z,i] = (2.*q)*h[x  ,y,    z,ibar]+\
                                              (1.-2.*q)*h[xpd,ypd,zpd,ibar]
                                h_str[x,y,z,i] -= 3.*2.*wi[ibar]*\
                                    (cxs[ibar]*uxLs[iL]+\
                                     cys[ibar]*uyLs[iL]+\
                                     czs[ibar]*uzLs[iL])
                            else: 
                                h_str[x,y,z,i] = 1./(2.*q) *h[x,y,z,ibar]+\
                                             (1.-1./(2.*q))*h[x,y,z,i]  
                                h_str[x,y,z,i] -= 3./q* wi[ibar]*\
                                    (cxs[ibar]*uxLs[iL]+\
                                     cys[ibar]*uyLs[iL]+\
                                     czs[ibar]*uzLs[iL])
                            FxLs[iL] = (h_str[x,y,z,i]+h[x,y,z,ibar])*cxs[i]                  
                            FyLs[iL] = (h_str[x,y,z,i]+h[x,y,z,ibar])*cys[i]                  
                            FzLs[iL] = (h_str[x,y,z,i]+h[x,y,z,ibar])*czs[i]
                        if i_BCs[x,y,z,i] == 5: # #Surface, small gap
                            ibar = ibars[i]
                            iL = PoiLs[x,y,z,i]
                            h_str[x,y,z,i] = h[x,y,z,ibar] -6*wi[ibar]*\
                                (cxs[ibar]*uxLs[iL]+\
                                 cys[ibar]*uyLs[iL]+\
                                 czs[ibar]*uzLs[iL])
                            FxLs[iL] = (h_str[x,y,z,i]+h[x,y,z,ibar])*cxs[i]                  
                            FyLs[iL] = (h_str[x,y,z,i]+h[x,y,z,ibar])*cys[i]                  
                            FzLs[iL] = (h_str[x,y,z,i]+h[x,y,z,ibar])*czs[i]                  
    for x in range(nx):
        for y in range(ny): 
            for z in range(nz): 
                if not OutsideLBMDomains[x,y,z]: 
                    h[x,y,z] = Relax.Relax_3D(h_str[x,y,z],cxs,cys,czs,wi,ibars,\
                    tauInvs[x,y,z],tauInvs_Asym[x,y,z],\
                    one_m_tauInvs_m_Iom[x,y,z],one_m_tauInvs_m_Iom_Asym[x,y,z])
    return h,Fx_on_Wall,FxLs,FyLs,FzLs

@jit(nopython = True)
def FreqDLBMStep_SoftPt_2D(h,nx,ny,\
        nd,cxs,cys,czs,wi,ibars,i_ups,i_notups,i_downs,i_notdowns,MatricesTop,\
        tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym):
    h_str      = np.zeros((nx,ny,nd),dtype = np.complex128)
    Fx_on_Wall = np.zeros((nx)      ,dtype = np.complex128)
    uxT,uyT = Handle_Top.Calc_uT_Fourier_2D(h,nx,MatricesTop)
    for x in range(nx):
        for y in range(1,ny-1):
            for i in range(nd):
                xmd = (nx+x-cxs[i])%nx 
                ymd =     y-cys[i]
                h_str[x,y,i] = h[xmd,ymd,i]
        for y in range(1):
            for i in i_notups: 
                xmd = (nx+x-cxs[i])%nx 
                ymd =     y-cys[i]
                h_str[x,y,i] = h[xmd,ymd,i]
            for i in i_ups:  
                ibar = ibars[i]
                h_str[x,y,i] = h[x,y,ibar] - 6*wi[ibar]*cxs[ibar] 
                Fx_on_Wall[x] += (h_str[x,y,i]+h[x,y,ibar])*cxs[i]                  
        for y in range(ny-1,ny):
            for i in i_notdowns: 
                xmd = (nx+x-cxs[i])%nx 
                ymd =     y-cys[i]
                h_str[x,y,i] = h[xmd,ymd,i]
            for i in i_downs:  
                ibar = ibars[i]
                h_str[x,y,i] = h[x,y,ibar] \
                    - 6 * wi[i]*cxs[ibar] * uxT[x]\
                    - 6 * wi[i]*cys[ibar] * uyT[x]
    for x in range(nx):
        for y in range(ny): 
            h[x,y] = Relax.Relax_2D(h_str[x,y],cxs,cys,czs,wi,ibars,\
            tauInvs[x,y],tauInvs_Asym[x,y],\
            one_m_tauInvs_m_Iom[x,y],one_m_tauInvs_m_Iom_Asym[x,y])
    return h,Fx_on_Wall

@jit(nopython = True)
def FreqDLBMStep_OscBnd_2D(h,nx,ny,nd,cxs,cys,czs,wi,ibars,MatricesTop,\
        nLstot,OutsideLBMDomains,i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
        tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym):
    h_str      = np.zeros((nx,ny,nd),dtype = np.complex128)
    Fx_on_Wall = np.zeros((nx)      ,dtype = np.complex128)
    FxLs       = np.zeros(nLstot       ,dtype = np.complex128) 
    FyLs       = np.zeros(nLstot       ,dtype = np.complex128) 
    uxT,uyT = Handle_Top.Calc_uT_Fourier_2D(h,nx,MatricesTop)
    for x in range(nx):
        for y in range(ny):
            if not OutsideLBMDomains[x,y]: 
                for i in range(nd):
                    if i_BCs[x,y,i] == 1: #bulk
                        xmd = (nx+x-cxs[i])%nx
                        ymd =     y-cys[i]
                        h_str[x,y,i] = h[xmd,ymd,i]
                    if i_BCs[x,y,i] == 2: #bottom
                        ibar = ibars[i] # i and ibar reversed compared to Krüger et al.
                        h_str[x,y,i] = h[x,y,ibar]-6*wi[ibar]*cxs[ibar]
                        Fx_on_Wall[x]+=(h_str[x,y,i]+h[x,y,ibar])*cxs[i]                  
                    if i_BCs[x,y,i] == 3: #top
                        ibar = ibars[i]
                        h_str[x,y,i] = h[x,y,ibar]\
                           -6*wi[ibar]*cxs[ibar]*uxT[x]\
                           -6*wi[ibar]*cys[ibar]*uyT[x]
                    if i_BCs[x,y,i] == 4: #Surface, no small gap
                        ibar = ibars[i]  
                        iL = PoiLs[x,y,i]
                        q = qs[iL]
                        if q <= 0.5: 
                            xpd = (nx+x+cxs[i])%nx
                            ypd =     y+cys[i]
                            h_str[x,y,i] = (2.*q)*h[x,y,ibar]+\
                                          (1.-2.*q)*h[xpd,ypd,ibar]
                            h_str[x,y,i] -= 3.*2.*wi[ibar]*\
                                (cxs[ibar]*uxLs[iL]+\
                                 cys[ibar]*uyLs[iL])
                        else: 
                            h_str[x,y,i] = 1./(2.*q) *h[x,y,ibar]+\
                                         (1.-1./(2.*q))*h[x,y,i]  
                            h_str[x,y,i] -= 3./q* wi[ibar]*\
                                (cxs[ibar]*uxLs[iL]+\
                                 cys[ibar]*uyLs[iL])
                        FxLs[iL] = (h_str[x,y,i]+h[x,y,ibar])*cxs[i]                  
                        FyLs[iL] = (h_str[x,y,i]+h[x,y,ibar])*cys[i]                  
                    if i_BCs[x,y,i] == 5: # #Surface, small gap
                        ibar = ibars[i]
                        iL = PoiLs[x,y,i]
                        h_str[x,y,i] = h[x,y,ibar] -6*wi[ibar]*\
                            (cxs[ibar]*uxLs[iL]+\
                             cys[ibar]*uyLs[iL]+\
                             czs[ibar]*uzLs[iL])
                        FxLs[iL] = (h_str[x,y,i]+h[x,y,ibar])*cxs[i]                  
                        FyLs[iL] = (h_str[x,y,i]+h[x,y,ibar])*cys[i]                  
    for x in range(nx):
        for y in range(ny): 
            if not OutsideLBMDomains[x,y]: 
                h[x,y] = Relax.Relax_2D(h_str[x,y],cxs,cys,czs,wi,ibars,\
                tauInvs[x,y],tauInvs_Asym[x,y],\
                one_m_tauInvs_m_Iom[x,y],one_m_tauInvs_m_Iom_Asym[x,y])
    return h,Fx_on_Wall,FxLs,FyLs

@jit(nopython = True)
def FreqDLBMStep_1D(h,ny,\
        nd,cxs,cys,wi,ibars,i_ups,i_notups,i_downs,i_notdowns,\
        tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym,om):
    h_str      = np.zeros((ny,nd),dtype = np.complex128)
    Fx_on_Wall = np.complex128(0)
    nu = (1./tauInvs[-1]-0.5)/3.
    ZLiq = (1j*om*nu)**0.5
    uxT = 2*(h[-1,5]-h[-1,6])/(ZLiq + 1./3.)   
    for y in range(1,ny-1):
        for i in range(nd):
            ymd = y-cys[i]
            h_str[y,i] = h[ymd,i]
    for y in range(1):
        for i in i_notups: 
            ymd = y-cys[i]
            h_str[y,i] = h[ymd,i]
        for i in i_ups:  
            ibar = ibars[i]
            h_str[y,i] = h[y,ibar] - 6*wi[ibar]*cxs[ibar] 
            Fx_on_Wall += (h_str[y,i]+h[y,ibar])*cxs[i]                  
    for y in range(ny-1,ny):
        for i in i_notdowns: 
            ymd =     y-cys[i]
            h_str[y,i] = h[ymd,i]
        for i in i_downs:  
            ibar = ibars[i]
            h_str[y,i] = h[y,ibar] - 6 *wi[i]*cxs[ibar] * uxT
    for y in range(ny): 
            h[y] = Relax.Relax_2D(h_str[y],cxs,cys,wi,ibars,\
            tauInvs[y],tauInvs_Asym[y],one_m_tauInvs_m_Iom[y],one_m_tauInvs_m_Iom_Asym[y])
    return h,Fx_on_Wall

# def FreqDLBMStep_OscBnd(h,nx,ny,nz,nd,cxs,cys,czs,wi,ibars,i_ups,i_notups,i_downs,i_notdowns,om,\
#         tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym,rhos,\
#         MatricesTop,nSph,UpdateMotionFac,OscBndLocked,OscBndLockedTo,n,\
#         OscBndAmps,SphPoss,SphRespPars,nLs,nLstot,iSs,xLs,yLs,zLs,OutsideLBMDomains,\
#         i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,iSiL_Lists,dimensions):
#     if dimensions == 2: 
#         h,Fx_on_Wall,FxLs,FyLs = \
#             FreqDLBMStep_OscBnd_3D(h,nx,ny,nz,nd,cxs,cys,wi,ibars,\
#                 om,MatricesTop,nLstot,OutsideLBMDomains,\
#                 i_BCs,PoiLs,qs,uxLs,uyLs,\
#                 tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym)
#     if dimensions == 3: 
#         h,Fx_on_Wall,FxLs,FyLs,FzLs = \
#             FreqDLBMStep_OscBnd_3D(h,nx,ny,nz,nd,cxs,cys,czs,wi,ibars,\
#                 om,MatricesTop,nLstot,OutsideLBMDomains,\
#                 i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
#                 tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym)
#     return h,Fx_on_Wall,FxLs,FyLs,FzLs            

@jit(nopython = True)
def FD_LBM_Step_Ref(h,ny,nd,cxs,cys,czs,wi,ibars,\
        i_ups,i_notups,i_downs,i_notdowns,ZBulk,om,\
        tauInv_Ref,tauInv_Asym_Ref,one_m_tauInv_m_Iom_Ref,one_m_tauInv_m_Iom_Asym_Ref):
    uxT = 2*(h[-1,7]-h[-1,14])/(ZBulk + 1./3.)   
    h_str = np.zeros((ny,nd),dtype = np.complex128)
    Fx_on_Wall= np.complex128(0)
    for y in range(1,ny-1): 
        for i in range(nd): h_str[y   ,i] = h[y-   cys[i],i] #bulk
    for i in i_notups     : h_str[0   ,i] = h[0   -cys[i],i] #bottom
    for i in i_notdowns   : h_str[ny-1,i] = h[ny-1-cys[i],i] #top
    for i in i_ups:                                          #bottom
        ibar = ibars[i]
        h_str[0,i] = h[0,ibar] - 6*wi[ibar]*cxs[ibar] 
        Fx_on_Wall += (h_str[0,i] + h[0,ibar]) * cxs[i]  
    for i in i_downs:                                        #top
        ibar = ibars[i]
        h_str[ny-1,i] = h[ny-1,ibar] - 6*wi[ibar]*cxs[ibar]*uxT
    for y in range(ny): 
        h[y] = Relax.Relax_3D(h_str[y],cxs,cys,czs,wi,ibars,\
            tauInv_Ref,tauInv_Asym_Ref,one_m_tauInv_m_Iom_Ref,one_m_tauInv_m_Iom_Asym_Ref)
    return h,Fx_on_Wall
