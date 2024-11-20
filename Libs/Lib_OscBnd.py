import numpy as np
from Libs import Lib_SetPos_3D as SetPos_3D
from Libs import Lib_General   as General
from Libs import Lib_Plots_from_Main as Plots_from_Main 
from Libs import Lib_Plots_for_GUI  as Plots

def Calc_Domains_3D(SPs):
    nx = SPs['nx']   
    ny = SPs['ny']    
    nz = SPs['nz']    
    nSph = SPs['nSph']
    SphPoss = SPs['SphPoss']
    RSph = SPs['RSph']
    OutsideLBMDomains = np.zeros((nx,ny,nz     ),dtype = np.int64)    
    InParticles       = np.zeros((nx,ny,nz,nSph),dtype = np.int64)    
    xSphs = SphPoss[0]
    ySphs = SphPoss[1]
    zSphs = SphPoss[2]
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for iS in range(nSph):
                    if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,xSphs[iS],\
                                                      y,ySphs[iS],\
                                                      z,zSphs[iS],RSph): 
                        OutsideLBMDomains[x,y,z]    = 1
                        InParticles[      x,y,z,iS] = 1
    return OutsideLBMDomains,InParticles

def Calc_qs_xLs_yLs_zLs_3D(nLstot,iSs,xGridLs,yGridLs,zGridLs,iGridLs,SphPoss,RSph,SPs): 
    nd,cxs,cys,czs,ibars,wi,i_ups,i_notups,i_downs,i_notdowns = General.ReadStencil(SPs['dimensions'])
    nx = SPs['nx']    
    ny = SPs['ny']    
    nz = SPs['nz']
    xSphs = SphPoss[0]
    ySphs = SphPoss[1]
    zSphs = SphPoss[2]
    qs  = np.zeros((nLstot),dtype = np.float64)
    xLs = np.zeros((nLstot),dtype = np.float64)
    yLs = np.zeros((nLstot),dtype = np.float64)
    zLs = np.zeros((nLstot),dtype = np.float64)
    for iL in range(nLstot):
        iS =iSs[iL] 
        x = xGridLs[iL]
        y = yGridLs[iL]
        z = zGridLs[iL]
        i = iGridLs[iL]
        xmd = x-cxs[i]
        ymd = y-cys[i]
        zmd = z-czs[i]
        qs[iL] = 1
        lobound = 0.; upbound = 1.
        for iteratations in range(12):
            threshold = (lobound + upbound) / 2.
            xt = (nx+x + threshold * (xmd-x))%nx
            yt =     y + threshold * (ymd-y)
            zt = (nz+z + threshold * (zmd-z))%nz
            Outside = SetPos_3D.DisLTDMin3D(nx,ny,nz,xt,xSphs[iS],
                                                     yt,ySphs[iS],
                                                     zt,zSphs[iS],RSph) 
            if Outside: upbound = threshold
            else      : lobound = threshold
        qs[iL] = threshold         
        xLs[iL] = x - qs[iL]*cxs[i]
        yLs[iL] = y - qs[iL]*cys[i]
        zLs[iL] = z - qs[iL]*czs[i]
    return qs,xLs,yLs,zLs

def Set_BoundaryPars_3D(OutsideLBMDomains,InParticles,SPs):
    SphPoss = SPs['SphPoss']
    nx = SPs['nx']    
    ny = SPs['ny']    
    nz = SPs['nz']    
    nd,cxs,cys,czs,ibars,wi,i_ups,i_notups,i_downs,i_notdowns = General.ReadStencil(SPs['dimensions'])
    nSph  = SPs['nSph']
    RSph  = SPs['RSph']
    i_BCs = np.zeros((nx,ny,nz,nd),dtype = np.int64)
    PoiLs = np.zeros((nx,ny,nz,nd),dtype = np.int64) 
    nLs   = np.zeros((nSph)       ,dtype = np.int64) 
    nLstot= np.int64(0)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if not OutsideLBMDomains[x,y,z]:
                    for i in range(nd):
                        xmd = (nx+x-cxs[i])%nx 
                        ymd =     y-cys[i]
                        zmd = (nz+z-czs[i])%nz 
                        xpd = (nx+x+cxs[i])%nx 
                        ypd =     y+cys[i]
                        zpd = (nz+z+czs[i])%nz 
                        if ymd <  0  : i_BCs[x,y,z,i] = 2 # Bottom
                        if ymd > ny-1: i_BCs[x,y,z,i] = 3 # Top
                        if ymd >= 0 and ymd <= ny-1:
                            if not OutsideLBMDomains[xmd,ymd,zmd]:i_BCs[x,y,z,i]=1 # Bulk
                            else:
                                for iS in range(nSph):
                                    if InParticles[xmd,ymd,zmd,iS]: 
                                        if ypd <= ny-1: 
                                            if not OutsideLBMDomains[xpd,ypd,zpd]: 
                                                i_BCs[x,y,z,i]=4 # Surface, no small gap
                                            else: i_BCs[x,y,z,i]=5 # Surface, small gap 
                                        else: i_BCs[x,y,z,i] = 5 # Surface, small gap       
                                        nLs[iS] += 1
    nLstot = 5 * np.sum(nLs)
    xGridLs = np.zeros((nLstot),dtype = np.int64)
    yGridLs = np.zeros((nLstot),dtype = np.int64)
    zGridLs = np.zeros((nLstot),dtype = np.int64)
    iGridLs = np.zeros((nLstot),dtype = np.int64)
    iSs     = np.zeros((nLstot),dtype = np.int64)
    count = 0
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if not OutsideLBMDomains[x,y,z]:
                    for i in range(nd):
                        xmd = (nx+x-cxs[i])%nx 
                        ymd =     y-cys[i]
                        zmd = (nz+z-czs[i])%nz 
                        if ymd < ny:
                            if OutsideLBMDomains[xmd,ymd,zmd]: 
                                for iS in range(nSph):
                                    if InParticles[xmd,ymd,zmd,iS]: 
                                        PoiLs[x,y,z,i] = count
                                        iSs[    count] = iS
                                        xGridLs[count] = x
                                        yGridLs[count] = y
                                        zGridLs[count] = z
                                        iGridLs[count] = i
                                        count += 1
    nLstot = count
    xGridLs = xGridLs[:nLstot]
    yGridLs = yGridLs[:nLstot]
    zGridLs = zGridLs[:nLstot]
    iGridLs = iGridLs[:nLstot]
    iSs     = iSs[    :nLstot]
                                    
    qs,xLs,yLs,zLs = Calc_qs_xLs_yLs_zLs_3D(nLstot,iSs,
        xGridLs,yGridLs,zGridLs,iGridLs,SphPoss,RSph,SPs)   
    iSiL_Lists = np.zeros((nSph,int(np.max(nLs))),dtype = int)
    icounts = np.zeros(nSph,dtype = int)
    for iL in range(np.sum(nLs)):
        iS = iSs[iL]
        iSiL_Lists[iS,icounts[iS]] = iL
        icounts[iS] += 1
      
    if SPs['Do_from_GUI'] :  
        Plots.Plot_LinkProps_3D( xLs,yLs,zLs,yLs*SPs['Dx_nm'],' ',SPs)    
    else :  
        Plots_from_Main.Plot_LinkProps_3D(xLs,yLs,zLs,yLs*SPs['Dx_nm'],xLs,xLs,' ','','',SPs)    
    
    return i_BCs,nLs,nLstot,iSs,PoiLs,xGridLs,yGridLs,zGridLs,iGridLs,\
        xLs,yLs,zLs,qs,iSiL_Lists

def Ini_Motion_3D(nSph,nLstot,SPs):
    uxSphs  = np.ones( nSph,dtype = np.complex128)
    uySphs  = np.zeros(nSph,dtype = np.complex128)
    uzSphs  = np.zeros(nSph,dtype = np.complex128)
    OmxSphs = np.zeros(nSph,dtype = np.complex128)
    OmySphs = np.zeros(nSph,dtype = np.complex128)
    OmzSphs = np.zeros(nSph,dtype = np.complex128)
    OscBndAmps = np.array([uxSphs,uySphs,uzSphs,OmxSphs,OmySphs,OmzSphs])
    uxLs    = np.ones( (nLstot),dtype = np.complex128)
    uyLs    = np.zeros((nLstot),dtype = np.complex128)
    uzLs    = np.zeros((nLstot),dtype = np.complex128)
    return OscBndAmps,uxLs,uyLs,uzLs

def Calc_SphRespPars_3D(SPs,OscBndPars): 
    nx        = SPs['nx']   
    ny        = SPs['ny']    
    nz        = SPs['nz']    
    RSph      = SPs['RSph']
    ySphbyR   = SPs['ySphbyR']
    etaabsSph = SPs['etaabsSph']
    tandelSph = SPs['tandelSph']
    rhoSph    = SPs['rhoSph']
    om        = SPs['om']
    ySph      = ySphbyR*RSph
    if ySphbyR > 1: 
        yCen = ySph
        MSph = 4./3.*np.pi*RSph**3*rhoSph
        IxxSph = 2./5.*MSph*RSph**2
        IyySph = 2./5.*MSph*RSph**2
        IzzSph = 2./5.*MSph*RSph**2
    else:     
        jmax = 100000 
        norm = 0. 
        yCen = 0; IxxSph = 0; IyySph = 0; IzzSph = 0  
        for j in range(jmax):
            x =      - RSph + 2*RSph*np.random.random()
            y = ySph - RSph + 2*RSph*np.random.random()
            z =      - RSph + 2*RSph*np.random.random()
            if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,0,y,ySph,z,0,RSph) and y > 0:
                yCen += y
                norm += 1 
        if norm > 0: 
            yCen /= norm;  MSph = norm/jmax*(2*RSph)**3*rhoSph
        else:        yCen = np.nan; MSph  = np.nan
        for j in range(jmax):
            x =      - RSph + 2*RSph*np.random.random()
            y = ySph - RSph + 2*RSph*np.random.random()
            z =      - RSph + 2*RSph*np.random.random()
            if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,0,y,ySph,z,0,RSph) and y > 0:
                r2Cenx = (z**2+(y-yCen)**2)**0.5   
                r2Ceny = (x**2+ z**2      )**0.5   
                r2Cenz = (x**2+(y-yCen)**2)**0.5   
                IxxSph += rhoSph * r2Cenx**2*(2*RSph)**2/jmax
                IyySph += rhoSph * r2Ceny**2*(2*RSph)**2/jmax
                IzzSph += rhoSph * r2Cenz**2*(2*RSph)**2/jmax
                
    cos = 1/(1.+tandelSph**2)**0.5
    sin = cos * tandelSph
    etaSph = etaabsSph * (cos - 1j*sin)*1./6.
    GSph   = etaSph    * 1j*om
    if ySph >= RSph : RCont = 0
    else            : RCont = RSph * (1-(ySph/RSph)**2.)**0.5 
    kappaShear = 2. * RCont * GSph
    kappaBend  = 2. * kappaShear * RCont**2/RSph**2

    xiContShear = kappaShear/(1j*om)
    xiContVertl = xiContShear
    xiContBendg = 2. * xiContShear * RCont**2/RSph**2
    xiContTwist = xiContBendg
    xiLiqTrans = 6. * np.pi * RSph
    xiLiqRotat = 8. * np.pi * RSph
    SphRespPars = np.array([yCen,MSph,IxxSph,IyySph,IzzSph,\
        RCont,kappaShear,kappaBend,xiContShear,xiContVertl,xiContBendg,xiContTwist,\
        xiLiqTrans,xiLiqRotat,etaabsSph,tandelSph],dtype=np.complex128)
    OscBndPars['SphRespPars'] = SphRespPars
    return 



"""
def Calc_SphRespPars_3D(SPs): 
    nx        = SPs['nx']   
    ny        = SPs['ny']    
    nz        = SPs['nz']    
    RSph      = SPs['RSph']
    ySphbyR   = SPs['ySphbyR']
    etaabsSph = SPs['etaabsSph']
    tandelSph = SPs['tandelSph']
    rhoSph    = SPs['rhoSph']
    om        = SPs['om']
    ySph      = ySphbyR*RSph
    if ySphbyR > 1: 
        yCen = ySph
        MSph = 4./3.*np.pi*RSph**3*rhoSph
        IxxSph = 2./5.*MSph*RSph**2
        IyySph = 2./5.*MSph*RSph**2
        IzzSph = 2./5.*MSph*RSph**2
    else:     
        jmax = 100000 
        norm = 0. 
        yCen = 0; IxxSph = 0; IyySph = 0; IzzSph = 0  
        for j in range(jmax):
            x =      - RSph + 2*RSph*np.random.random()
            y = ySph - RSph + 2*RSph*np.random.random()
            z =      - RSph + 2*RSph*np.random.random()
            if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,0,y,ySph,z,0,RSph) and y > 0:
                yCen += y
                norm += 1 
        if norm > 0: 
            yCen /= norm;  MSph = norm/jmax*(2*RSph)**3*rhoSph
        else:        yCen = np.nan; MSph  = np.nan
        for j in range(jmax):
            x =      - RSph + 2*RSph*np.random.random()
            y = ySph - RSph + 2*RSph*np.random.random()
            z =      - RSph + 2*RSph*np.random.random()
            if SetPos_3D.DisLTDMin3D(nx,ny,nz,x,0,y,ySph,z,0,RSph) and y > 0:
                r2Cenx = (z**2+(y-yCen)**2)**0.5   
                r2Ceny = (x**2+ z**2      )**0.5   
                r2Cenz = (x**2+(y-yCen)**2)**0.5   
                IxxSph += rhoSph * r2Cenx**2*(2*RSph)**2/jmax
                IyySph += rhoSph * r2Ceny**2*(2*RSph)**2/jmax
                IzzSph += rhoSph * r2Cenz**2*(2*RSph)**2/jmax
                
    cos = 1/(1.+tandelSph**2)**0.5
    sin = cos * tandelSph
    etaSph = etaabsSph * (cos - 1j*sin)*1./6.
    GSph   = etaSph    * 1j*om
    if ySph >= RSph : RCont = 0
    else            : RCont = RSph * (1-(ySph/RSph)**2.)**0.5 
    kappaShear = 2. * RCont * GSph
    kappaBend  = 2. * kappaShear * RCont**2/RSph**2

    xiContShear = kappaShear/(1j*om)
    xiContVertl = xiContShear
    xiContBendg = 2. * xiContShear * RCont**2/RSph**2
    xiContTwist = xiContBendg
    xiLiqTrans = 6. * np.pi * RSph
    xiLiqRotat = 8. * np.pi * RSph
    SphRespPars = np.array([yCen,MSph,IxxSph,IyySph,IzzSph,\
        RCont,kappaShear,kappaBend,xiContShear,xiContVertl,xiContBendg,xiContTwist,\
        xiLiqTrans,xiLiqRotat,etaabsSph,tandelSph],dtype=np.complex128)
    return SphRespPars
"""
def Update_Motion_3D(nLs,nLstot,iSs,xLs,yLs,zLs,uxLs,uyLs,uzLs,\
        FxLs,FyLs,FzLs,nSph,om,\
        OscBndAmps,SphPoss,SphRespPars,\
        UpdateMotionFac,OscBndLocked,OscBndLockedTo,iSiL_Lists,n,SPs):
    xSphs = SphPoss[0]
    ySphs = SphPoss[1]
    zSphs = SphPoss[2]
    Fx = np.zeros((nSph),dtype = np.complex128)
    Fy = np.zeros((nSph),dtype = np.complex128)
    Fz = np.zeros((nSph),dtype = np.complex128)
    Tx = np.zeros((nSph),dtype = np.complex128)
    Ty = np.zeros((nSph),dtype = np.complex128)
    Tz = np.zeros((nSph),dtype = np.complex128)
    FxCont  = np.zeros((nSph),dtype = np.complex128)
    FyCont  = np.zeros((nSph),dtype = np.complex128)
    FzCont  = np.zeros((nSph),dtype = np.complex128)
    TxCont  = np.zeros((nSph),dtype = np.complex128)
    TyCont  = np.zeros((nSph),dtype = np.complex128)
    TzCont  = np.zeros((nSph),dtype = np.complex128)
    RCont       = SphRespPars[5]
    kappaShear  = SphRespPars[6]
    kappaBend   = SphRespPars[7]
    RContSI     = RCont * SPs['Dx_nm']*1e-9
    Mass_SI = 1e-24*SPs['Dx_nm']**3
    om_SI = 2.*np.pi*SPs['n']*SPs['f0_SI']; 
    Dt_SI = SPs['om']/om_SI
    kappaShearSI = kappaShear * Mass_SI / Dt_SI**2
    kappaBendSI  = kappaBend * Mass_SI / Dt_SI**2
    xiContShear = SphRespPars[8]
    xiContVertl = SphRespPars[9]
    xiContBendg = SphRespPars[10]
    xiContTwist = SphRespPars[11]
    xiLiqTrans  = SphRespPars[12]
    xiLiqRotat  = SphRespPars[13]
    etaabsSph   = SphRespPars[14]
    tandelSph   = SphRespPars[15]
    if not OscBndLocked : 
        uxSphs  = OscBndAmps[0]
        uySphs  = OscBndAmps[1]
        uzSphs  = OscBndAmps[2]
        OmxSphs = OscBndAmps[3]
        OmySphs = OscBndAmps[4]
        OmzSphs = OscBndAmps[5]
        yCen   = SphRespPars[0]
        MSph   = SphRespPars[1]
        IxxSph = SphRespPars[2]
        IyySph = SphRespPars[3]
        IzzSph = SphRespPars[4]
        
        for iS in range(nSph):
            Fx_Liq = 0; Fy_Liq = 0; Fz_Liq = 0 
            Tx_Liq = 0; Ty_Liq = 0; Tz_Liq = 0
            for iL in iSiL_Lists[iS,:nLs[iS]]:
                iL = int(iL)
                Fx_Liq += -FxLs[iL]
                Fy_Liq += -FyLs[iL]
                Fz_Liq += -FzLs[iL]
                Tx_Liq += ((yLs[iL]-ySphs[iS])*(-FzLs[iL])-(zLs[iL]-zSphs[iS])*(-FyLs[iL]))   
                Ty_Liq += ((zLs[iL]-zSphs[iS])*(-FxLs[iL])-(xLs[iL]-xSphs[iS])*(-FzLs[iL]))   
                Tz_Liq += ((xLs[iL]-xSphs[iS])*(-FyLs[iL])-(yLs[iL]-ySphs[iS])*(-FxLs[iL]))
            ux_at_ResSurf = uxSphs[iS] - OmzSphs[iS]*yCen
            uy_at_ResSurf = uySphs[iS]
            uz_at_ResSurf = uzSphs[iS] + OmxSphs[iS]*yCen      
            FxCont[iS] =  xiContShear * (1.-ux_at_ResSurf)
            FyCont[iS] =  xiContVertl * (0.-uy_at_ResSurf)
            FzCont[iS] =  xiContShear * (0.-uz_at_ResSurf)
            TxCont[iS] = -xiContBendg * OmxSphs[iS]*yCen**2 + FzCont[iS]*yCen     
            TyCont[iS] = -xiContTwist * OmySphs[iS]*yCen**2
            TzCont[iS] = -xiContBendg * OmzSphs[iS]*yCen**2 - FxCont[iS]*yCen     
            Fx_inertia = -1j*om*MSph   * uxSphs[ iS]
            Fy_inertia = -1j*om*MSph   * uySphs[ iS]
            Fz_inertia = -1j*om*MSph   * uzSphs[ iS]
            Tx_inertia = -1j*om*IxxSph * OmxSphs[iS]
            Ty_inertia = -1j*om*IyySph * OmySphs[iS]
            Tz_inertia = -1j*om*IzzSph * OmzSphs[iS]
            Fx[iS] = Fx_Liq + Fx_inertia + FxCont[iS]    
            Fy[iS] = Fy_Liq + Fy_inertia + FyCont[iS] 
            Fz[iS] = Fz_Liq + Fz_inertia + FzCont[iS] 
            Tx[iS] = Tx_Liq + Tx_inertia + TxCont[iS]
            Ty[iS] = Ty_Liq + Ty_inertia + TyCont[iS]
            Tz[iS] = Tz_Liq + Tz_inertia + TzCont[iS]
            uxSphs[iS]  += UpdateMotionFac/n * Fx[iS] / (xiLiqTrans + xiContShear)
            uySphs[iS]  += UpdateMotionFac/n * Fy[iS] / (xiLiqTrans + xiContVertl)
            uzSphs[iS]  += UpdateMotionFac/n * Fz[iS] / (xiLiqTrans + xiContShear)
            OmxSphs[iS] += UpdateMotionFac/n * Tx[iS] / (xiLiqRotat + xiContBendg)
            OmySphs[iS] += UpdateMotionFac/n * Ty[iS] / (xiLiqRotat + xiContTwist)
            OmzSphs[iS] += UpdateMotionFac/n * Tz[iS] / (xiLiqRotat + xiContBendg)
            for iL in iSiL_Lists[iS,:nLs[iS]]:
                iL = int(iL)
                uxLs[iL] = uxSphs[iS] + OmySphs[iS]*(zLs[iL]-zSphs[iS]) - \
                                        OmzSphs[iS]*(yLs[iL]-ySphs[iS]) 
                uyLs[iL] = uySphs[iS] + OmzSphs[iS]*(xLs[iL]-xSphs[iS]) - \
                                        OmxSphs[iS]*(zLs[iL]-zSphs[iS])
                uzLs[iL] = uzSphs[iS] + OmxSphs[iS]*(yLs[iL]-ySphs[iS]) - \
                                        OmySphs[iS]*(xLs[iL]-xSphs[iS])
        uxSphm1        = np.average(uxSphs)-1
        OmzSph         = np.average(OmzSphs) 
        if RCont > 0 :  
            sigContactpol  = np.average(Fx) / (np.pi*RCont**2)
            TzContbyAreaxR = np.average(Tz) / (np.pi*RCont**2)*(yCen)
        else : 
            sigContactpol  = np.nan
            TzContbyAreaxR = np.nan
    if OscBndLocked : 
        uySphs  = np.zeros(nSph,dtype = np.complex128)
        uzSphs  = np.zeros(nSph,dtype = np.complex128)
        OmxSphs = np.zeros(nSph,dtype = np.complex128)
        OmySphs = np.zeros(nSph,dtype = np.complex128)
        OmzSphs = np.zeros(nSph,dtype = np.complex128)
        for iS in range(nSph):
            Fx_Liq = 0; Fy_Liq = 0; Fz_Liq = 0 
            Tx_Liq = 0; Ty_Liq = 0; Tz_Liq = 0
            for iL in iSiL_Lists[iS,:nLs[iS]]:
                iL = int(iL)
                Fx_Liq += -FxLs[iL]
                Fy_Liq += -FyLs[iL]
                Fz_Liq += -FzLs[iL]
                Tx_Liq += ((yLs[iL]-ySphs[iS])*(-FzLs[iL])-(zLs[iL]-zSphs[iS])*(-FyLs[iL]))   
                Ty_Liq += ((zLs[iL]-zSphs[iS])*(-FxLs[iL])-(xLs[iL]-xSphs[iS])*(-FzLs[iL]))   
                Tz_Liq += ((xLs[iL]-xSphs[iS])*(-FyLs[iL])-(yLs[iL]-ySphs[iS])*(-FxLs[iL]))
            Fx[iS] = Fx_Liq
            Fy[iS] = Fy_Liq  
            Fz[iS] = Fz_Liq
            Tx[iS] = Tx_Liq
            Ty[iS] = Ty_Liq
            Tz[iS] = Tz_Liq
        OmzSph         = 0
        sigContactpol      = np.nan
        TzContbyAreaxR = np.nan
        if OscBndLockedTo == 'Zero':
            uxSphs  = np.zeros(nSph,dtype = np.complex128)
            uxSphm1 = 1
        if OscBndLockedTo == 'Substrate':
            uxSphs  = np.ones(nSph,dtype = np.complex128)
            uxSphm1 = 0

    OscBndAmps[0] = uxSphs
    OscBndAmps[1] = uySphs
    OscBndAmps[2] = uzSphs
    OscBndAmps[3] = OmxSphs
    OscBndAmps[4] = OmySphs
    OscBndAmps[5] = OmzSphs
    Fx_avg = np.average(Fx)
    Tz_avg = np.average(Tz)
    FxCont_tot = np.sum(FxCont)
    MotionPars = np.array([uxSphm1,\
                           OmzSph,\
                           sigContactpol,\
                           TzContbyAreaxR,\
                           Fx_avg,\
                           Tz_avg])
    MotionParsTitles = ['$u_{x,Sph}-1$',\
                        '$\Omega_{z,Sph}$',\
                        '$\sigma_{xy}$',\
                        'Torque/Area/$R_{Sph}$',\
                        '$F_{x,avg}$',\
                        '$T_{z,avg}$']
    AuxPars = np.array([etaabsSph,\
                        tandelSph,\
                        RContSI,\
                        kappaShearSI,\
                        kappaBendSI,\
                        np.nan])
    AuxParsTitles =    ['$|\eta_{abs,Sph}|$',\
                        '$tan(\delta)del_{Sph}$',\
                        '$R_{Cont,SI}$',\
                        '$\kappa_{Shear_SI}$',\
                        '$\kappa_{Bend,SI}$',\
                        '']
    return uxLs,uyLs,uzLs,FxCont_tot,OscBndAmps,MotionPars,MotionParsTitles,AuxPars,AuxParsTitles

def Collect_OscBndPars(iSs,nLs,nLstot,
        xLs,yLs,zLs,OutsideLBMDomains,InParticles,i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
        OscBndAmps,UpdateMotionFac,OscBndLocked,OscBndLockedTo,
        nSph,RSph,ySphbyR,rhoSph,delta0_nm,Dx_nm,iSiL_Lists,\
        xGridLs,yGridLs,zGridLs,iGridLs):
    OscBndPars = {}
    OscBndPars.update({'nSph'        : nSph}) 
    OscBndPars.update({'RSph'        : RSph}) 
    OscBndPars.update({'ySphbyR'     : ySphbyR}) 
    OscBndPars.update({'rhoSph'      : rhoSph}) 
    OscBndPars.update({'iSs'         : iSs}) 
    OscBndPars.update({'nLs'         : nLs}) 
    OscBndPars.update({'nLstot'      : nLstot}) 
    OscBndPars.update({'xLs'         : xLs}) 
    OscBndPars.update({'yLs'         : yLs}) 
    OscBndPars.update({'zLs'         : zLs}) 
    OscBndPars.update({'OutsideLBMDomains' : OutsideLBMDomains}) 
    OscBndPars.update({'InParticles' : InParticles}) 
    OscBndPars.update({'i_BCs'       : i_BCs}) 
    OscBndPars.update({'PoiLs'       : PoiLs}) 
    OscBndPars.update({'qs'          : qs}) 
    OscBndPars.update({'uxLs'        : uxLs}) 
    OscBndPars.update({'uyLs'        : uyLs}) 
    OscBndPars.update({'uzLs'        : uzLs}) 
    OscBndPars.update({'OscBndAmps'  : OscBndAmps}) 
    # OscBndPars.update({'SphRespPars' : SphRespPars}) 
    OscBndPars.update({'UpdateMotionFac': UpdateMotionFac}) 
    OscBndPars.update({'OscBndLocked': OscBndLocked}) 
    OscBndPars.update({'OscBndLockedTo': OscBndLockedTo}) 
    OscBndPars.update({'iSiL_Lists'  : iSiL_Lists}) 
    OscBndPars.update({'delta0_nm'   : delta0_nm}) 
    OscBndPars.update({'Dx_nm'       : Dx_nm}) 
    OscBndPars.update({'xGridLs'     : xGridLs}) 
    OscBndPars.update({'yGridLs'     : yGridLs}) 
    OscBndPars.update({'zGridLs'     : zGridLs}) 
    OscBndPars.update({'iGridLs'     : iGridLs}) 

    return OscBndPars


# def Collect_OscBndPars(iSs,nLs,nLstot,
#         xLs,yLs,zLs,OutsideLBMDomains,InParticles,i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
#         OscBndAmps,SphRespPars,UpdateMotionFac,OscBndLocked,OscBndLockedTo,
#         nSph,RSph,ySphbyR,rhoSph,delta0_nm,Dx_nm,iSiL_Lists):
#     OscBndPars = {}
#     OscBndPars.update({'nSph'        : nSph}) 
#     OscBndPars.update({'RSph'        : RSph}) 
#     OscBndPars.update({'ySphbyR'     : ySphbyR}) 
#     OscBndPars.update({'rhoSph'      : rhoSph}) 
#     OscBndPars.update({'iSs'         : iSs}) 
#     OscBndPars.update({'nLs'         : nLs}) 
#     OscBndPars.update({'nLstot'      : nLstot}) 
#     OscBndPars.update({'xLs'         : xLs}) 
#     OscBndPars.update({'yLs'         : yLs}) 
#     OscBndPars.update({'zLs'         : zLs}) 
#     OscBndPars.update({'OutsideLBMDomains' : OutsideLBMDomains}) 
#     OscBndPars.update({'InParticles' : InParticles}) 
#     OscBndPars.update({'i_BCs'       : i_BCs}) 
#     OscBndPars.update({'PoiLs'       : PoiLs}) 
#     OscBndPars.update({'qs'          : qs}) 
#     OscBndPars.update({'uxLs'        : uxLs}) 
#     OscBndPars.update({'uyLs'        : uyLs}) 
#     OscBndPars.update({'uzLs'        : uzLs}) 
#     OscBndPars.update({'OscBndAmps'  : OscBndAmps}) 
#     OscBndPars.update({'SphRespPars' : SphRespPars}) 
#     OscBndPars.update({'UpdateMotionFac': UpdateMotionFac}) 
#     OscBndPars.update({'OscBndLocked': OscBndLocked}) 
#     OscBndPars.update({'OscBndLockedTo': OscBndLockedTo}) 
#     OscBndPars.update({'iSiL_Lists'  : iSiL_Lists}) 
#     OscBndPars.update({'delta0_nm'   : delta0_nm}) 
#     OscBndPars.update({'Dx_nm'       : Dx_nm}) 
#     return OscBndPars
    
def Extract_OscBndPars(OscBndPars):
    nSph              = int(OscBndPars['nSph'])
    RSph              = OscBndPars['RSph'] 
    ySphbyR           = OscBndPars['ySphbyR'] 
    rhoSph            = OscBndPars['rhoSph']
    iSs               = OscBndPars['iSs']
    nLs               = OscBndPars['nLs']
    nLstot            = OscBndPars['nLstot']
    xLs               = OscBndPars['xLs']
    yLs               = OscBndPars['yLs']
    zLs               = OscBndPars['zLs']
    OutsideLBMDomains = OscBndPars['OutsideLBMDomains']
    InParticles       = OscBndPars['InParticles']
    i_BCs             = OscBndPars['i_BCs']
    PoiLs             = OscBndPars['PoiLs']
    qs                = OscBndPars['qs']
    uxLs              = OscBndPars['uxLs']
    uyLs              = OscBndPars['uyLs']
    uzLs              = OscBndPars['uzLs']
    OscBndAmps        = OscBndPars['OscBndAmps']
    SphRespPars       = OscBndPars['SphRespPars']
    UpdateMotionFac   = OscBndPars['UpdateMotionFac']
    OscBndLocked      = OscBndPars['OscBndLocked']
    OscBndLockedTo    = OscBndPars['OscBndLockedTo']

    iSiL_Lists    = OscBndPars['iSiL_Lists']

    return iSs,nLs,nLstot,xLs,yLs,zLs,OutsideLBMDomains,InParticles,i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
        OscBndAmps,SphRespPars,UpdateMotionFac,OscBndLocked,OscBndLockedTo,\
        nSph,RSph,ySphbyR,rhoSph,iSiL_Lists
        

def Setup_Boundaries_3D(SPs):
    nSph            = int(SPs['nSph'])
    UpdateMotionFac = SPs['UpdateMotionFac']
    RSph            = SPs['RSph']
    ySphbyR         = SPs['ySphbyR']
    rhoSph          = SPs['rhoSph']
    delta0_nm       = SPs['delta0_nm']
    Dx_nm           = SPs['Dx_nm']
    OscBndLocked    = SPs['OscBndLocked']
    OscBndLockedTo  = SPs['OscBndLockedTo']

    OutsideLBMDomains,InParticles = Calc_Domains_3D(SPs)

    i_BCs,nLs,nLstot,iSs,PoiLs,xGridLs,yGridLs,zGridLs,iGridLs,\
    xLs,yLs,zLs,qs,iSiL_Lists = \
        Set_BoundaryPars_3D(OutsideLBMDomains,InParticles,SPs)    

    OscBndAmps,uxLs,uyLs,uzLs = Ini_Motion_3D(nSph,nLstot,SPs)

    OscBndPars = Collect_OscBndPars(iSs,nLs,nLstot,
        xLs,yLs,zLs,OutsideLBMDomains,InParticles,i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
        OscBndAmps,UpdateMotionFac,OscBndLocked,OscBndLockedTo,
        nSph,RSph,ySphbyR,rhoSph,delta0_nm,Dx_nm,iSiL_Lists,\
        xGridLs,yGridLs,zGridLs,iGridLs)
    return OscBndPars

# def Setup_Boundaries_3D(SPs):
#     nSph            = int(SPs['nSph'])
#     UpdateMotionFac = SPs['UpdateMotionFac']
#     RSph            = SPs['RSph']
#     ySphbyR         = SPs['ySphbyR']
#     rhoSph          = SPs['rhoSph']
#     delta0_nm       = SPs['delta0_nm']
#     Dx_nm           = SPs['Dx_nm']
#     OscBndLocked    = SPs['OscBndLocked']
#     OscBndLockedTo  = SPs['OscBndLockedTo']

#     OutsideLBMDomains,InParticles = Calc_Domains_3D(SPs)

#     i_BCs,nLs,nLstot,iSs,PoiLs,xGridLs,yGridLs,zGridLs,iGridLs,\
#     xLs,yLs,zLs,qs,iSiL_Lists = \
#         Set_BoundaryPars_3D(OutsideLBMDomains,InParticles,SPs)

#     # SphRespPars = Calc_SphRespPars_3D(SPs)

#     OscBndAmps,uxLs,uyLs,uzLs = Ini_Motion_3D(nSph,nLstot,SPs)

#     OscBndPars = Collect_OscBndPars(iSs,nLs,nLstot,
#         xLs,yLs,zLs,OutsideLBMDomains,InParticles,i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
#         OscBndAmps,SphRespPars,UpdateMotionFac,OscBndLocked,OscBndLockedTo,
#         nSph,RSph,ySphbyR,rhoSph,delta0_nm,Dx_nm,iSiL_Lists)
#     return OscBndPars

def Calc_dr_ux_uy_uz_Rigid_StiffParticle(x,y,z,InParticles,OscBndAmps,SPs):
    nSph = int(SPs['nSph'])
    SphPoss = SPs['SphPoss']
    xSphs = SphPoss[0]
    ySphs = SphPoss[1]
    zSphs = SphPoss[2]
    uxSphs  = OscBndAmps[0]
    uySphs  = OscBndAmps[1]
    uzSphs  = OscBndAmps[2]
    OmxSphs = OscBndAmps[3]
    OmySphs = OscBndAmps[4]
    OmzSphs = OscBndAmps[5]
    drl = 0; uxl = 0; uyl = 0; uzl = 0
    for iS in range(nSph):
        if InParticles[x,y,z,iS]: 
            drl = 0
            uxl = uxSphs[iS] + OmySphs[iS]*(z-zSphs[iS]) - \
                               OmzSphs[iS]*(y-ySphs[iS])
            uyl = uySphs[iS] + OmzSphs[iS]*(x-xSphs[iS]) - \
                               OmxSphs[iS]*(z-zSphs[iS])
            uzl = uzSphs[iS] + OmxSphs[iS]*(y-ySphs[iS]) - \
                               OmySphs[iS]*(x-xSphs[iS])
    return drl,uxl,uyl,uzl

def Calc_dr_ux_uy_uz_OscBnd(h,cxs,cys,czs,OutsideLBMDomains,InParticles,OscBndAmps,SPs):
    nx = SPs['nx']
    ny = SPs['ny']
    nz = SPs['nz']
    dr = np.zeros((nx,ny,nz),dtype = np.complex128)
    ux = np.zeros((nx,ny,nz),dtype = np.complex128)
    uy = np.zeros((nx,ny,nz),dtype = np.complex128)
    uz = np.zeros((nx,ny,nz),dtype = np.complex128)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if not OutsideLBMDomains[x,y,z]: 
                    dr[x,y,z] = np.sum(h[x,y,z])
                    ux[x,y,z] = np.sum(h[x,y,z]*cxs)
                    uy[x,y,z] = np.sum(h[x,y,z]*cys)
                    uz[x,y,z] = np.sum(h[x,y,z]*czs)
                else: 
                    if SPs['ProblemType'] == 'StiffParticles':
                        dr[x,y,z],ux[x,y,z],uy[x,y,z],uz[x,y,z] = \
                            Calc_dr_ux_uy_uz_Rigid_StiffParticle(x,y,z,\
                                InParticles,OscBndAmps,SPs)
                    if SPs['ProblemType'] == 'Roughness':
                        dr[x,y,z] = uy[x,y,z] = uz[x,y,z] = 0
                        ux[x,y,z] = 1
                    if SPs['ProblemType'] == 'SFA':
                        dr[x,y,z] = ux[x,y,z] = uy[x,y,z] = uz[x,y,z] = 0
    return dr,ux,uy,uz                

def Calc_Domains_2D(SPs):
    nx = SPs['nx']   
    ny = SPs['ny']    
    RoughnessFourierComps = SPs['RoughnessFourierComps']
    OutsideLBMDomains = np.zeros((nx,ny),dtype = np.int64)    
    RoughSurface = np.zeros(nx)
    for x in range(nx):
        for y in range(ny):
            if y < RoughSurface[x]: OutsideLBMDomains[x,y]    = 1
    return OutsideLBMDomains

def Calc_qs_xLs_yLs_zLs_2D(SPs,nLstot,xGridLs,yGridLs,iGridLs): 
    nd,cxs,cys,czs,ibars,wi,i_ups,i_notups,i_downs,i_notdowns = General.ReadStencil(SPs['dimensions'])
    nx = SPs['nx']    
    RoughSurface = np.zeros(nx)
    qs  = np.zeros((nLstot),dtype = np.float64)
    xLs = np.zeros((nLstot),dtype = np.float64)
    yLs = np.zeros((nLstot),dtype = np.float64)
    for iL in range(nLstot):
        x = xGridLs[iL]
        y = yGridLs[iL]
        i = iGridLs[iL]
        xmd = x-cxs[i]
        ymd = y-cys[i]
        qs[iL] = 1
        lobound = 0.; upbound = 1.
        for iteratations in range(12):
            threshold = (lobound + upbound) / 2.
            xt = (nx+x + threshold * (xmd-x))%nx
            yt =     y + threshold * (ymd-y)
            Outside = yt < RoughSurface[xt]
            if Outside: upbound = threshold
            else      : lobound = threshold
        qs[iL] = threshold         
        xLs[iL] = x - qs[iL]*cxs[i]
        yLs[iL] = y - qs[iL]*cys[i]
    return qs,xLs,yLs

def Set_BoundaryPars_2D(OutsideLBMDomains,SPs):
    nx = SPs['nx']    
    ny = SPs['ny']    
    nd,cxs,cys,czs,ibars,wi,i_ups,i_notups,i_downs,i_notdowns = General.ReadStencil(SPs['dimensions'])
    i_BCs = np.zeros((nx,ny,nd),dtype = np.int64)
    PoiLs = np.zeros((nx,ny,nd),dtype = np.int64) 
    nLstot= np.int64(0)
    for x in range(nx):
        for y in range(ny):
            if not OutsideLBMDomains[x,y]:
                for i in range(nd):
                    xmd = (nx+x-cxs[i])%nx 
                    ymd =     y-cys[i]
                    xpd = (nx+x+cxs[i])%nx 
                    ypd =     y+cys[i]
                    if ymd <  0  : i_BCs[x,y,i] = 2 # Bottom
                    if ymd > ny-1: i_BCs[x,y,i] = 3 # Top
                    if ymd >= 0 and ymd <= ny-1:
                        if not OutsideLBMDomains[xmd,ymd]:i_BCs[x,y,i]=1 # Bulk
                        else:
                            if ypd <= ny-1: 
                                if not OutsideLBMDomains[xpd,ypd]: 
                                    i_BCs[x,y,i]=4 # Surface, no small gap
                                else: i_BCs[x,y,i]=5 # Surface, small gap 
                            else: i_BCs[x,y,i] = 5 # Surface, small gap       
                            nLstot += 1
    xGridLs = np.zeros((nLstot),dtype = np.int64)
    yGridLs = np.zeros((nLstot),dtype = np.int64)
    iGridLs = np.zeros((nLstot),dtype = np.int64)
    count = 0
    for x in range(nx):
        for y in range(ny):
            if not OutsideLBMDomains[x,y]:
                for i in range(nd):
                    xmd = (nx+x-cxs[i])%nx 
                    ymd =     y-cys[i]
                    if ymd < ny:
                        if OutsideLBMDomains[xmd,ymd]: 
                            PoiLs[x,y,i] = count
                            xGridLs[count] = x
                            yGridLs[count] = y
                            iGridLs[count] = i
                            count += 1
    xGridLs = xGridLs[:nLstot]
    yGridLs = yGridLs[:nLstot]
    iGridLs = iGridLs[:nLstot]
                                    
    qs,xLs,yLs = Calc_qs_xLs_yLs_zLs_2D(SPs,nLstot,xGridLs,yGridLs,iGridLs)   
    #Plots.Plot_LinkProps_2D(xLs,yLs,yLs*SPs['Dx_nm'],yLs,yLs,' ','','')
    return i_BCs,nLstot,PoiLs,xGridLs,yGridLs,iGridLs,xLs,yLs,qs

def Update_Motion_2D(nLstot,FxLs):
    FxLiq = -np.sum(FxLs)
    uxLs = np.ones( nLstot,dtype = np.complex128)
    uyLs = np.zeros(nLstot,dtype = np.complex128)
    MotionPars = np.array([np.nan]); MotionParsTitles = [''];     
    AuxPars    = np.array([np.nan]); AuxParsTitles    = ['']; 
    return uxLs,uyLs,FxLiq,MotionPars,MotionParsTitles,AuxPars,AuxParsTitles
