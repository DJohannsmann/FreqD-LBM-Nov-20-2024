import numpy as np

def Calc_Dfcbyn(SPs,Fx_on_Wall,FxCont,Do_Ref):
    Dx_SI = 1e-9*SPs['Dx_nm']; Mass_SI = 1e-24*SPs['Dx_nm']**3
    om_SI = 2.*np.pi*SPs['n']*SPs['f0_SI']; Dt_SI = SPs['om']/om_SI
    if Do_Ref : Stress_avg =  np.sum(Fx_on_Wall)
    else      : Stress_avg = (np.sum(Fx_on_Wall)+FxCont)/(SPs['nx']*SPs['nz'])
    Stress_SI   = Stress_avg*Mass_SI /(Dx_SI*Dt_SI**2)
    Velocity_SI = Dx_SI/Dt_SI
    ZLoad_SI    = Stress_SI/Velocity_SI
    Dfcbyn = SPs['f0_SI']*1j/(np.pi*SPs['Zq_SI'])*ZLoad_SI/SPs['n']
    return Dfcbyn

def Calc_tauInvBulk_ZBulk(SPs):
    cosBulk = 1/(1.+SPs['tandelBulk']**2)**0.5
    sinBulk = cosBulk * SPs['tandelBulk']
    etaBulk = (-1j*cosBulk + sinBulk)*1./6.
    nuBulk  = etaBulk
    SPs['tauInvBulk'] = 1./(3.*nuBulk+0.5)
    SPs['ZBulk']      = (1j*SPs['om']*etaBulk)**0.5

def Calc_etaabstandel(SPs):
    ncen = 1./2.*(np.max(SPs['ns'])+np.min(SPs['ns']))
    if  SPs['VEPars_Choice'] == 'etaabs_tandel':
        SPs['etaabsSph'] = SPs['etaabscenSphmPas'] * (SPs['n']/ncen)**SPs['betap_Sph' ] * 1./6.
        SPs['tandelSph'] = SPs['tandelcenSph']     * (SPs['n']/ncen)**SPs['betapp_Sph']
    if SPs['VEPars_Choice'] == 'from_J':
        eta_SI = 1e-3;
        omcen_SI = 2.*np.pi*ncen     *SPs['f0_SI']
        om_SI    = 2.*np.pi*SPs['n'] *SPs['f0_SI']
        Jp_cenSI = 1./(eta_SI * SPs['Jp_FacSph'] * omcen_SI)
        JppcenSI = 1./(eta_SI * SPs['Jpp_FacSph'] * omcen_SI)
        Jp_SI = Jp_cenSI * (SPs['n']/ncen)**SPs['betap_Sph']
        JppSI = JppcenSI * (SPs['n']/ncen)**SPs['betapp_Sph']
        GabsSphSI = 1./np.abs(Jp_SI+1j*JppSI)
        etaabsSphSI = GabsSphSI /(om_SI)
        SPs['etaabsSph']   = etaabsSphSI/1e-3*1./6.
        if Jp_SI > 1e-33 : SPs['tandelSph'] = JppSI/Jp_SI
        else             : SPs['tandelSph'] = 1e33
    if  SPs['VEPars_Choice'] == 'Maxwell':
        Rate_SI = SPs['MaxwellRelaxRate_MHz']*1e6
        tau_SI = 1./Rate_SI
        eta_SI = 1e-3;
        omcen_SI = 2.*np.pi*ncen     *SPs['f0_SI']
        om_SI    = 2.*np.pi*SPs['n'] *SPs['f0_SI']
        GinfSI = omcen_SI*SPs['etaabscenSphmPas']*1e-3
        GpSI  = om_SI**2 * tau_SI**2 / (1+om_SI**2*tau_SI**2) * GinfSI
        GppSI = om_SI    * tau_SI    / (1+om_SI**2*tau_SI**2) * GinfSI
        etapSI  = GppSI/om_SI + 1e-3*6
        etappSI = GpSI /om_SI 
        etaabsSI = np.abs(etapSI + 1j*etappSI)
        SPs['etaabsSph'] = etaabsSI / 1e-3 * 1./6.
        SPs['tandelSph'] = GppSI / GpSI
    return SPs

def ReadStencil(dimensions):
    if dimensions == 1 or dimensions == 2: 
        nd = 9; w0 = 4./9.; ws = 1./9.; wd = 1./36.
        cxs    = np.array([0,1,0,-1, 0,1,-1,-1, 1], dtype = np.int64)     
        cys    = np.array([0,0,1, 0,-1,1, 1,-1,-1], dtype = np.int64)  
        czs    = np.zeros(nd, dtype = np.int64)
        ibars  = np.array([0,3,4, 1, 2,7, 8, 5, 6], dtype = np.int64) 
        wi    = np.array([w0,ws,ws,ws,ws,wd,wd,wd,wd], dtype = np.float64)
        i_ups      = np.array([2,5,6])  # needed for bounce back at bottom
        i_notups   = np.array([0,1,3,4,7,8])
        i_downs    = np.array([4,7,8])  # needed for bounce back at top
        i_notdowns = np.array([0,1,2,3,5,6])
    if dimensions == 3 : 
        nd = 19; w0 = 1./3.; ws = 1./18.; wd = 1./36.
        cxs    = np.array([0,1,-1,0, 0,0, 0,1,-1 ,1 ,-1, 0, 0, 1,-1, 1,-1, 0, 0], dtype = np.int64)
        cys    = np.array([0,0,0, 1,-1,0, 0,1,-1 ,0 , 0, 1,-1,-1, 1, 0, 0, 1,-1], dtype = np.int64)
        czs    = np.array([0,0,0, 0, 0,1,-1,0, 0 ,1 ,-1, 1,-1, 0, 0,-1, 1,-1, 1], dtype = np.int64)
        ibars  = np.array([0,2, 1,4, 3,6, 5,8, 7 ,10, 9,12,11,14,13,16,15,18,17], dtype = np.int64)
        wi     = np.array([w0,ws,ws,ws,ws,ws,ws,wd,wd,wd,wd,wd,wd,wd,wd,wd,wd,wd,wd],dtype = np.float64)
        i_ups      = np.array([3,7,14,11,17])  # needed for bounce back at bottom
        i_notups   = np.array([0,1,2,4,5,6,8,9,10,12,13,15,16,18])
        i_downs    = np.array([4,8,13,12,18])  # needed for bounce back at top
        i_notdowns = np.array([0,1,2,3,5,6,7,9,10,11,14,15,16,17])
    return nd,cxs,cys,czs,ibars,wi,i_ups,i_notups,i_downs,i_notdowns
    
