import numpy as np
from Libs import Lib_SetPos_3D as SetPos_3D
from numba import jit

def Calc_yCen_3D(SPs): 
    nx        = SPs['nx']   
    ny        = SPs['ny']    
    nz        = SPs['nz']    
    RSph      = SPs['RSph']
    ySphbyR   = SPs['ySphbyR']
    ySph      = ySphbyR*RSph
    if ySphbyR > 1: yCen = ySph
    else:     
        jmax = 100000 
        norm = 0. 
        yCen = 0; 
        for j in range(jmax):
            x =      - RSph + 2*RSph*np.random.random()
            y = ySph - RSph + 2*RSph*np.random.random()
            z =      - RSph + 2*RSph*np.random.random()
            if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,0,y,ySph,z,0,RSph) and y > 0:
                yCen += y
                norm += 1 
        if norm > 0: 
            yCen /= norm;
    SPs['yCen'] = yCen        
    return      

def Calc_FracVolSph_1D(SPs):  
    nx = int(SPs['nx'])
    ny = int(SPs['ny'])
    nz = int(SPs['nz'])
    FilmThickness = SPs['FilmThickness_nm'] / SPs['Dx0_nm']
    FracVolSph = np.zeros((nx,ny,nz),dtype=np.float64)
    for y in range(ny):
        if  y+1 <= FilmThickness : FracVolSph[0,y,0] = 1
        elif  y >= FilmThickness : FracVolSph[0,y,0] = 0
        elif  y <  FilmThickness and y+1 >= FilmThickness: FracVolSph[0,y,0] = FilmThickness - y
    return FracVolSph      

def Calc_FracVolSph_2D(SPs):  
    nx = int(SPs['nx'])
    ny = int(SPs['ny'])
    nz = int(SPs['nz'])
    FracVolSph = np.zeros((nx,ny,nz),dtype=np.float64)
    return FracVolSph      

def Calc_FracVolSph_3D(SPs):  
    nx = int(SPs['nx'])
    ny = int(SPs['ny'])
    nz = int(SPs['nz'])
    nSph       = int(SPs['nSph'])
    RSph       = SPs['RSph']
    SphPoss    = SPs['SphPoss']
    FracVolSph = np.zeros((nx,ny,nz),dtype=np.float64)
    xSphs = SphPoss[0]
    ySphs = SphPoss[1]
    zSphs = SphPoss[2]
    n_for_integration = 1000 
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for iS in range(nSph):
                    if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
                                                      y,ySphs[iS],\
                                                      z,zSphs[iS],RSph-2): 
                        FracVolSph[x,y,z] = 1
                    if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
                                                      y,ySphs[iS],z,zSphs[iS],RSph+2)\
                        and not SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
                                                               y,ySphs[iS],\
                                                               z,zSphs[iS],RSph-2): 
                        for j in range(n_for_integration):
                            x2 = np.float64(x+np.random.rand()-0.5)
                            y2 = np.float64(y+np.random.rand()-0.5)
                            z2 = np.float64(z+np.random.rand()-0.5)
                            if SetPos_3D.DisLTDMin3D(nx,ny,nz,x2,xSphs[iS],\
                                                              y2,ySphs[iS],\
                                                              z2,zSphs[iS],RSph): 
                                FracVolSph[x,y,z] += 1./n_for_integration
    return FracVolSph

def Set_RelaxPars(SPs):
    nx = int(SPs['nx'])
    ny = int(SPs['ny'])
    nz = int(SPs['nz'])
    tauInvBulk  = SPs['tauInvBulk']
    etaabsSph  = SPs['etaabsSph']
    tandelSph  = SPs['tandelSph']
    Lambda_TRT = SPs['Lambda_TRT']
    om         = SPs['om']
    rhoSph     = SPs['rhoSph']
    dimensions = SPs['dimensions']
    cosSph = 1/(1.+tandelSph**2)**0.5
    sinSph = tandelSph * cosSph
    etaSph = etaabsSph * (-1j*cosSph + sinSph)
    nuSph  = etaSph
    nuBulk = (1./(tauInvBulk)-0.5)/3.
    if dimensions == 1 : FracVolSph = Calc_FracVolSph_1D(SPs)
    if dimensions == 2 : FracVolSph = Calc_FracVolSph_2D(SPs)
    if dimensions == 3 :
        if SPs['ProblemType'] in ['StiffParticles','SFA']:
            FracVolSph = np.zeros((nx,ny,nz),dtype=np.float64)
        if SPs['ProblemType'] == 'SoftParticles' :
            FracVolSph = Calc_FracVolSph_3D(SPs)
    nus = (np.exp(  FracVolSph *np.log(nuSph**0.5 ) + \
                 (1-FracVolSph)*np.log(nuBulk**0.5)) )**2
    rhos = FracVolSph *rhoSph + (1-FracVolSph)
    tauInvs = np.complex128(1./(3.*nus+0.5))
    if Lambda_TRT != 0 :  
        tauInvs_Asym = (4.-2.*tauInvs)/(2.-tauInvs+4*Lambda_TRT*tauInvs)
    else : tauInvs_Asym = tauInvs    
    one_m_tauInvs_m_Iom      = 1 - tauInvs      - 1j*om*rhos + 5./2.*(om*rhos)**2
    one_m_tauInvs_m_Iom_Asym = 1 - tauInvs_Asym - 1j*om*rhos  + 5./2.*(om*rhos)**2
    return FracVolSph,tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym,rhos

def Calc_sigxy_sigyy_3D(hl,cxs,cys,czs,wi,tauInvl):
    drl = np.sum(hl    )
    uxl = np.sum(hl*cxs)
    uyl = np.sum(hl*cys)
    uzl = np.sum(hl*czs)
    cdotu = cxs*uxl + cys*uyl + czs*uzl
    heq   = wi * (drl + 3.*cdotu)
    sigxy = -(1.-0.5*tauInvl)*np.sum(cxs*cys*(hl-heq))
    sigyy = -(1.-0.5*tauInvl)*np.sum(cys*cys*(hl-heq))
    return sigxy,sigyy
    
def Calc_MotionPars_3D_UU(B, nx,ny,nz,nd,cxs,cys,czs,i_ups,i_notups,ibars,wi,h,\
                 tauInvs,FracVolSph,SPs,SphPoss):
    
    nSph      = int(SPs['nSph'])
    RSph      = SPs['RSph']
    etaabsSph = SPs['etaabsSph']
    tandelSph = SPs['tandelSph']
    
    xSphs = SphPoss[0]
    ySphs = SphPoss[1]
    zSphs = SphPoss[2]
    
    ux = np.zeros((nx,ny,nz) ,dtype = np.complex128)
    uy = np.zeros((nx,ny,nz) ,dtype = np.complex128)
    uz = np.zeros((nx,ny,nz) ,dtype = np.complex128)
    
    uxSphs    = np.zeros((nSph),dtype = np.complex128)
    uySphs    = np.zeros((nSph),dtype = np.complex128)
    uzSphs    = np.zeros((nSph),dtype = np.complex128)
    
    OmxSphsxR = np.zeros((nSph),dtype = np.complex128)
    OmySphsxR = np.zeros((nSph),dtype = np.complex128)
    OmzSphsxR = np.zeros((nSph),dtype = np.complex128)
    
    sigxySphs = np.zeros((nSph),dtype = np.complex128)
    TzbyAreabyRSphs = np.zeros((nSph),dtype = np.complex128)
 
    norms   = np.zeros((nSph),dtype = np.complex128)
    
    # could be further optimized...
    def normalize_nsph(arr, norms):
        for i in range(nSph):
            if norms[i] > 0 : arr[i] /= norms[i]
            else            : arr[i] = np.nan      
        return arr
    
    # no optimization over the spheres for now
    uxl = np.sum(h*cxs, axis=3)
    uyl = np.sum(h*cys, axis=3)
    uzl = np.sum(h*czs, axis=3)
    
    for i in range(nSph):
        uxSphs[i] = np.sum(uxl*FracVolSph, axis=(0,1,2), where=(B[..., i]==1))
        uySphs[i] = np.sum(uyl*FracVolSph, axis=(0,1,2), where=(B[..., i]==1))
        uzSphs[i] = np.sum(uzl*FracVolSph, axis=(0,1,2), where=(B[..., i]==1))
        
        norms[i] = np.sum(FracVolSph, axis=(0,1,2), where=(B[..., i]==1))
    
    # double checked these: they are correct
    ux = np.sum(h*cxs, axis=3)
    uy = np.sum(h*cys, axis=3)
    uz = np.sum(h*czs, axis=3)
    
    uxSphs = np.where(norms > 0, uxSphs/norms, np.nan)
    uySphs = np.where(norms > 0, uySphs/norms, np.nan)
    uzSphs = np.where(norms > 0, uzSphs/norms, np.nan)
    
    eijS = np.zeros((3,3), dtype = complex);
    norms  = np.zeros((nSph),dtype = np.complex128)
    
    B_ = np.copy(B)
    B_[:, 0, :, :] = 0
    
    xm = [(nx+x-1) % nx for x in range(nx)]
    xp = [(nx+x+1) % nx for x in range(nx)]
    
    ym = [max(y-1, 0) for y in range(ny)]
    yp = [min(y+1, ny-1) for y in range(ny)]
    
    zm = [(nz+z-1) % nz for z in range(nz)]
    zp = [(nz+z+1) % nz for z in range(nz)]
    
    for iS in range(nSph):
        eijS[0,0] += np.sum((ux[xp, ...] - ux[xm, ...]) * FracVolSph / 2, where=(B_[..., iS]==1))
        eijS[1,0] += np.sum((uy[xp, ...] - uy[xm, ...]) * FracVolSph / 2, where=(B_[..., iS]==1))
        eijS[2,0] += np.sum((uz[xp, ...] - uz[xm, ...]) * FracVolSph / 2, where=(B_[..., iS]==1))
        
        eijS[0,1] += np.sum((ux[:, yp, :] - ux[:, ym, :]) * FracVolSph / 2, where=(B_[..., iS]==1))
        eijS[1,1] += np.sum((uy[:, yp, :] - uy[:, ym, :]) * FracVolSph / 2, where=(B_[..., iS]==1))
        eijS[2,1] += np.sum((uz[:, yp, :] - uz[:, ym, :]) * FracVolSph / 2, where=(B_[..., iS]==1))
        
        eijS[0,2] += np.sum((ux[..., zp] - ux[..., zm]) * FracVolSph / 2, where=(B_[..., iS]==1))
        eijS[1,2] += np.sum((uy[..., zp] - uy[..., zm]) * FracVolSph / 2, where=(B_[..., iS]==1))
        eijS[2,2] += np.sum((uz[..., zp] - uz[..., zm]) * FracVolSph / 2, where=(B_[..., iS]==1))
        
        norms[iS] = np.sum(FracVolSph, axis=(0,1,2), where=(B_[..., iS]==1))
    
    # Static things
    OmxSphsxR[iS] = 0.5*(eijS[1,2] - eijS[2,1])*RSph
    OmySphsxR[iS] = 0.5*(eijS[2,0] - eijS[0,2])*RSph 
    OmzSphsxR[iS] = 0.5*(eijS[0,1] - eijS[1,0])*RSph
    
    OmxSphsxR = normalize_nsph(OmxSphsxR, norms)
    OmySphsxR = normalize_nsph(OmySphsxR, norms)
    OmzSphsxR = normalize_nsph(OmzSphsxR, norms)

    h_str      = np.zeros((nx,ny,nz,nd),dtype = np.complex128)

    # TODO: optimize this in the future
    for x in range(nx):
        for z in range(nz):
            for y in range(1):
                for i in i_notups: 
                    xmd = (nx+x-cxs[i]) % nx 
                    ymd =     y-cys[i]
                    zmd = (nz+z-czs[i]) % nz 
                    h_str[x,y,z,i] = h[xmd,ymd,zmd,i]
                
                for ibar in i_ups:
                    i = ibars[ibar]
                    h_str[x,y,z,ibar] = h[x,y,z,i] - 6*wi[i]*cxs[i] 

    norms   = np.zeros((nSph),dtype = np.complex128)
    for x in range(nx):
        for z in range(nz):
            sigxyl,sigyyl = Calc_sigxy_sigyy_3D(h_str[x,0,z],cxs,cys,czs,wi,tauInvs[x,0,z])
            for iS in range(nSph): 
                if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
                                               0,ySphs[iS],\
                                               z,zSphs[iS],RSph+2):
                    sigxySphs[iS]       +=               sigxyl*FracVolSph[x,0,z]
                    TzbyAreabyRSphs[iS] += (x-xSphs[iS])*sigyyl*FracVolSph[x,0,z]/RSph
                    norms[iS]           +=                      FracVolSph[x,0,z]
    
    
    # from here it's just summarizing the results
    sigxySphs = normalize_nsph(sigxySphs, norms)
    TzbyAreabyRSphs = normalize_nsph(TzbyAreabyRSphs, norms)
    
    uxSphm1     = np.nanmean(uxSphs)-1
    OmzSphxR    = np.nanmean(OmzSphsxR)
    sigxy       = np.nanmean(sigxySphs) 
    TzbyAreabyR = np.nanmean(TzbyAreabyRSphs)

    FxSphs_avg,TzSphs_avg = Calc_Force_Torque_3D(SPs,h,cxs,cys,czs,wi,tauInvs,FracVolSph)
    
    MotionPars = np.array([uxSphm1,OmzSphxR,sigxy,TzbyAreabyR,FxSphs_avg,TzSphs_avg])
    MotionParsTitles =   ['$u_{x,Sph}-1$','$\Omega_{z,Sph}$ $\\times$ $R_{Sph}$',\
                       '$\sigma_{xy}$','Torque/Area/$R_{Sph}$',\
                        '$F_{x,avg}$','$T_{z,avg}$']
    AuxPars = np.array([etaabsSph , tandelSph,np.nan,np.nan,np.nan,np.nan])
    AuxParsTitles =   ['etaabsSph','tandelSph','','','','']
    return MotionPars,MotionParsTitles,AuxPars,AuxParsTitles


# def Calc_MotionPars_3D(nx,ny,nz,nd,cxs,cys,czs,i_ups,i_notups,ibars,wi,h,\
#                  tauInvs,FracVolSph,SPs,SphPoss):
#     nSph      = int(SPs['nSph'])
#     RSph      = SPs['RSph']
#     etaabsSph = SPs['etaabsSph']
#     tandelSph = SPs['tandelSph']
    
#     xSphs = SphPoss[0]
#     ySphs = SphPoss[1]
#     zSphs = SphPoss[2]
#     ux = np.zeros((nx,ny,nz) ,dtype = np.complex128)
#     uy = np.zeros((nx,ny,nz) ,dtype = np.complex128)
#     uz = np.zeros((nx,ny,nz) ,dtype = np.complex128)
#     uxSphs    = np.zeros((nSph),dtype = np.complex128)
#     uySphs    = np.zeros((nSph),dtype = np.complex128)
#     uzSphs    = np.zeros((nSph),dtype = np.complex128)
#     OmxSphsxR = np.zeros((nSph),dtype = np.complex128)
#     OmySphsxR = np.zeros((nSph),dtype = np.complex128)
#     OmzSphsxR = np.zeros((nSph),dtype = np.complex128)
#     sigxySphs = np.zeros((nSph),dtype = np.complex128)
#     TzbyAreabyRSphs = np.zeros((nSph),dtype = np.complex128)
 
#     norms   = np.zeros((nSph),dtype = np.complex128)
#     for x in range(nx):
#         for y in range(ny):
#             for z in range(nz):
#                 ux[x,y,z] = np.sum(h[x,y,z]*cxs)
#                 uy[x,y,z] = np.sum(h[x,y,z]*cys)
#                 uz[x,y,z] = np.sum(h[x,y,z]*czs)
#                 for iS in range(nSph): 
#                     if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
#                                                       y,ySphs[iS],\
#                                                       z,zSphs[iS],RSph+2):
#                         uxl = np.sum(h[x,y,z]*cxs)
#                         uyl = np.sum(h[x,y,z]*cys)
#                         uzl = np.sum(h[x,y,z]*czs)
#                         uxSphs[iS] += uxl*FracVolSph[x,y,z]
#                         uySphs[iS] += uyl*FracVolSph[x,y,z]
#                         uzSphs[iS] += uzl*FracVolSph[x,y,z]
#                         norms[iS]  +=     FracVolSph[x,y,z]
#     for iS in range(nSph): 
#         if norms[iS] > 0: 
#             uxSphs[iS] /= norms[iS]
#             uySphs[iS] /= norms[iS]
#             uzSphs[iS] /= norms[iS]
#         else: 
#             uxSphs[iS] = np.nan
#             uySphs[iS] = np.nan
#             uzSphs[iS] = np.nan
    
#     eijS = np.zeros((3,3), dtype = complex);     
#     norms  = np.zeros((nSph),dtype = np.complex128)
#     for x in range (nx):
#         for y in range (1,ny-1):
#             for z in range (nz):
#                 for iS in range(nSph): 
#                     if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
#                                                       y,ySphs[iS],\
#                                                       z,zSphs[iS],RSph+2):
#                         xm = (nx+x-1) % nx; 
#                         xp = (nx+x+1) % nx; 
#                         eijS[0,0] += (ux[xp,y,z] - ux[xm,y,z])*FracVolSph[x,y,z]/2 
#                         eijS[1,0] += (uy[xp,y,z] - uy[xm,y,z])*FracVolSph[x,y,z]/2  
#                         eijS[2,0] += (uz[xp,y,z] - uz[xm,y,z])*FracVolSph[x,y,z]/2  
    
#                         ym = y-1
#                         yp = y+1
#                         eijS[0,1] += (ux[x,yp,z] - ux[x,ym,z])*FracVolSph[x,y,z]/2.  
#                         eijS[1,1] += (uy[x,yp,z] - uy[x,ym,z])*FracVolSph[x,y,z]/2.  
#                         eijS[2,1] += (uz[x,yp,z] - uz[x,ym,z])*FracVolSph[x,y,z]/2.  
    
#                         zm = (nz+z-1) % nz; 
#                         zp = (nz+z+1) % nz; 
#                         eijS[0,2] += (ux[x,y,zp] - ux[x,y,zm])*FracVolSph[x,y,z]/2.  
#                         eijS[1,2] += (uy[x,y,zp] - uy[x,y,zm])*FracVolSph[x,y,z]/2.  
#                         eijS[2,2] += (uz[x,y,zp] - uz[x,y,zm])*FracVolSph[x,y,z]/2.  
#                         norms[iS]  +=     FracVolSph[x,y,z]
#     OmxSphsxR[iS] = 0.5*(eijS[1,2] - eijS[2,1])*RSph 
#     OmySphsxR[iS] = 0.5*(eijS[2,0] - eijS[0,2])*RSph 
#     OmzSphsxR[iS] = 0.5*(eijS[0,1] - eijS[1,0])*RSph 
#     for iS in range(nSph): 
#         if norms[iS] > 0: 
#             OmxSphsxR[iS] /= norms[iS]
#             OmySphsxR[iS] /= norms[iS]
#             OmzSphsxR[iS] /= norms[iS]
#         else: 
#             OmxSphsxR[iS] = np.nan
#             OmySphsxR[iS] = np.nan
#             OmzSphsxR[iS] = np.nan

#     h_str      = np.zeros((nx,ny,nz,nd),dtype = np.complex128)

#     for x in range(nx):
#         for z in range(nz):
#             for y in range(1):
#                 for i in i_notups: 
#                     xmd = (nx+x-cxs[i]) % nx 
#                     ymd =     y-cys[i]
#                     zmd = (nz+z-czs[i]) % nz 
#                     h_str[x,y,z,i] = h[xmd,ymd,zmd,i]
#                 for ibar in i_ups:  
#                     i = ibars[ibar]
#                     h_str[x,y,z,ibar] = h[x,y,z,i] - 6*wi[i]*cxs[i] 

#     norms   = np.zeros((nSph),dtype = np.complex128)
#     for x in range(nx):
#         for z in range(nz):
#             sigxyl,sigyyl = Calc_sigxy_sigyy_3D(h_str[x,0,z],cxs,cys,czs,wi,tauInvs[x,0,z])
#             for iS in range(nSph): 
#                 if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
#                                                   0,ySphs[iS],\
#                                                   z,zSphs[iS],RSph+2):
#                     sigxySphs[iS]       +=               sigxyl*FracVolSph[x,0,z]
#                     TzbyAreabyRSphs[iS] += (x-xSphs[iS])*sigyyl*FracVolSph[x,0,z]/RSph
#                     norms[iS]           +=                      FracVolSph[x,0,z]
#     for iS in range(nSph): 
#         if norms[iS] > 0: 
#             sigxySphs[iS]       /= norms[iS]
#             TzbyAreabyRSphs[iS] /= norms[iS]
#         else: 
#             sigxySphs[iS]       = np.nan
#             TzbyAreabyRSphs[iS] = np.nan
#     uxSphm1     = np.nanmean(uxSphs)-1 
#     OmzSphxR    = np.nanmean(OmzSphsxR) 
#     sigxy       = np.nanmean(sigxySphs) 
#     TzbyAreabyR = np.nanmean(TzbyAreabyRSphs) 
    
#     FxSphs_avg,TzSphs_avg = Calc_Force_Torque_3D(SPs,h,cxs,cys,czs,wi,tauInvs,FracVolSph)
    
#     MotionPars = np.array([uxSphm1,OmzSphxR,sigxy,TzbyAreabyR,FxSphs_avg,TzSphs_avg])
#     MotionParsTitles = ['$u_{\mathrm{x,Sph}}-1$',\
#                         '$\Omega_{\mathrm{z,Sph}}$ $\\times$ $R_{\mathrm{Sph}}$',\
#                         '$\sigma_{\mathrm{xy}}$','Torque/Area/$R_{\mathrm{Sph}}$',\
#                         '$F_{\mathrm{x,avg}}$','$T_{\mathrm{z,avg}}$']
#     AuxPars = np.array([etaabsSph , tandelSph,np.nan,np.nan,np.nan,np.nan])
#     AuxParsTitles =   ['etaabsSph','tandelSph','','','','']
#     return MotionPars,MotionParsTitles,AuxPars,AuxParsTitles

def Calc_dr_ux_uy_uz_SoftPt_3D(h,cxs,cys,czs,nx,ny,nz):
    dr = np.sum(h    , axis=3, dtype=np.complex128)
    ux = np.sum(h*cxs, axis=3, dtype=np.complex128)
    uy = np.sum(h*cys, axis=3, dtype=np.complex128)
    uz = np.sum(h*czs, axis=3, dtype=np.complex128)
    return dr,ux,uy,uz   

def Calc_Force_Torque_3D(SPs,h,cxs,cys,czs,wi,tauInvs,FracVolSph):
    def Calc_sig(hl,cxs,cys,czs,wi,tauInvl):
        sigl = np.zeros((3,3),dtype = np.complex128)
        cs = np.array([cxs,cys,czs])
        drl = np.sum(hl    )
        uxl = np.sum(hl*cxs)
        uyl = np.sum(hl*cys)
        uzl = np.sum(hl*czs)
        cdotu = cxs*uxl + cys*uyl + czs*uzl
        heq   = wi * (drl + 3.*cdotu)
        for xyz in range(3):
            for xyzp in range(3):
                sigl[xyz,xyzp] = -(1.-0.5*tauInvl)*np.sum(cs[xyz]*cs[xyzp]*(hl-heq))
        return sigl
    
    Calc_yCen_3D(SPs)
    nx = int(SPs['nx'])
    ny = int(SPs['ny'])
    nz = int(SPs['nz'])
    nSph    = int(SPs['nSph'])
    RSph    = SPs['RSph']
    SphPoss = SPs['SphPoss']
    ycen = SPs['yCen'] 
    xSphs = SphPoss[0]
    ySphs = SphPoss[1]
    zSphs = SphPoss[2]
    sig = np.zeros((nx,ny,nz,3,3),dtype = np.complex128)
    Fxs  = np.zeros((nx,ny,nz),dtype = np.complex128)
    Fys  = np.zeros((nx,ny,nz),dtype = np.complex128)
    Fzs  = np.zeros((nx,ny,nz),dtype = np.complex128)
    FxSphs = np.zeros((nSph),dtype = np.complex128)
    TzSphs = np.zeros((nSph),dtype = np.complex128)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                sig[x,y,z] = Calc_sig(h[x,y,z],cxs,cys,czs,wi,tauInvs[x,y,z])
    for x in range(nx):
        for y in range(1,ny-1):
            for z in range(nz):
                xm = (nx+x-1) % nx; 
                xp = (nx+x+1) % nx;
                for jdir in range(3):
                    Fxs[x,y,z] += (sig[xp,y,z,0,jdir] - sig[xm,y,z,0,jdir])/2 
                for jdir in range(3):
                    Fxs[x,y,z] += (sig[xp,y,z,jdir,0] - sig[xm,y,z,jdir,0])/2 

                ym = y-1
                yp = y+1
                for jdir in range(3):
                    Fys[x,y,z] += (sig[x,yp,z,1,jdir] - sig[x,ym,z,1,jdir])/2 
                for jdir in range(3):
                    Fys[x,y,z] += (sig[x,yp,z,jdir,1] - sig[x,ym,z,jdir,1])/2 

                zm = (nz+z-1) % nz; 
                zp = (nz+z+1) % nz; 
                for jdir in range(3):
                    Fzs[x,y,z] += (sig[x,y,zp,2,jdir] - sig[x,y,zm,2,jdir])/2 
                for jdir in range(3):
                    Fzs[x,y,z] += (sig[x,y,zp,jdir,2] - sig[x,y,zm,jdir,2])/2 

                for iS in range(nSph):
                    if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
                                                      y,ySphs[iS],\
                                                      z,zSphs[iS],RSph+1): 
                        FxSphs[iS] += Fxs[x,y,z]*FracVolSph[x,y,z]
                        TzSphs[iS] += ((Fys[x,y,z]*(x-xSphs[iS]))-\
                                       (Fxs[x,y,z]*(y-ycen)     ))*FracVolSph[x,y,z]
    FxSphs_avg = np.average(FxSphs)
    TzSphs_avg = np.average(TzSphs)
    return FxSphs_avg,TzSphs_avg

@jit(nopython=True, cache=True)
def compute_B(nx, ny, nz, nSph, RSph, xSphs, ySphs, zSphs):
    B = np.zeros((nx,ny,nz,nSph), dtype=np.int16)
    for x in range(nx): 
        for y in range(ny):
            for z in range(nz):
                for iS in range(nSph): 
                    if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
                                                y,ySphs[iS],\
                                                z,zSphs[iS],RSph+2):
                        B[x,y,z,iS] = 1
    return B
