import numpy as np
from Libs import Lib_General   as General
from Libs import Lib_RingIn    as RingIn
from Libs import Lib_Soft      as Soft
from Libs import Lib_OscBnd    as OscBnd
from Libs import Lib_SetPos_3D as SetPos_3D

def Handle_Geometry_FilmResonance(SPs): # FilmResonance option not debugged
    SPs['Width_nm']     = np.nan
    SPs['RSph']         = np.nan
    SPs['nx'] = SPs['nz'] = 1
    SPs['ny']           = int(np.round(SPs['FilmThickness']/SPs['Dx_nm'])) + 3
    SPs['SphPoss']      = [[np.nan],[np.nan],[np.nan]]
    SPs['nNodes']       = SPs['ny']
    SPs['CoverageTrue'] = np.nan

def Handle_Geometry_Roughness(SPs):  # Roughness option not debugged
    SPs['Width_nm']     = SPs['Roughn_Width_nm']
    SPs['RSph']         = np.nan
    nx_float            = SPs['Width_nm']/SPs['Dx_nm']
    SPs['nx']           = int(nx_float/2)*2+1  #odd
    SPs['nz']           = 1
    SPs['Roughn_VertScale'] = SPs['Roughn_VertScale_nm'] / SPs['Dx_nm']
    SPs['ny']           = int(np.round(2.*SPs['Roughn_VertScale'])+\
                              2.*SPs['Roughn_VertScale']*SPs['Gap2TopbyR'])
    SPs['SphPoss']      = [[np.nan],[np.nan],[np.nan]]
    SPs['nNodes']       = SPs['nx']*SPs['ny']*SPs['nz']
    SPs['CoverageTrue'] = np.nan

def Handle_Geometry_SFA(SPs): # SFA option not debugged
    SPs['Width_nm']     = SPs['Roughn_Width_nm']
    SPs['RSph']         = np.nan
    nx_float            = SPs['Width_nm']/SPs['Dx_nm']
    SPs['nx']           = int(nx_float/2)*2+1  #odd
    SPs['nz']           = 1
    SPs['Roughn_VertScale'] = SPs['Roughn_VertScale_nm'] / SPs['Dx_nm']
    SPs['ny']           = int(np.round(2.*SPs['Roughn_VertScale'])+\
                              2.*SPs['Roughn_VertScale']*SPs['Gap2TopbyR'])
    SPs['SphPoss']      = [[np.nan],[np.nan],[np.nan]]
    SPs['nNodes']       = SPs['nx']*SPs['ny']*SPs['nz']
    SPs['CoverageTrue'] = np.nan

def Handle_Geometry_Spheres(SPs):
    SPs['Width_nm']     = (SPs['nSph']*np.pi*SPs['RSph_nm']**2/SPs['CovTarget'])**0.5
    SPs['RSph']         = SPs['RSph_nm']/SPs['Dx_nm']
    nx_float            = SPs['Width_nm']/SPs['Dx_nm']
    SPs['nx']=SPs['nz'] = int(nx_float/2)*2+1  #odd
    SPs['ny']           = int(np.round(SPs['RSph']*SPs['ySphbyR']+\
                              SPs['RSph']*(1.+SPs['Gap2TopbyR'])))
    SPs['SphPoss']      = SetPos_3D.Set_SphPoss(SPs)
    SPs['nNodes']       = SPs['nx']*SPs['ny']*SPs['nz']
    SPs['CoverageTrue'] = SPs['nSph']*np.pi*SPs['RSph_nm']**2/(SPs['nx']*SPs['Dx_nm'])**2

def SingleSimulation(SPs):
    # ns=SPs['ns']
    # if SPs['ProblemType'] in ['SoftParticles','StiffParticles'] :  Handle_Geometry_Spheres(SPs)
    # if SPs['ProblemType'] == 'SFA'            : Handle_Geometry_SFA(SPs)
    # if SPs['ProblemType'] == 'Roughness'      : Handle_Geometry_Roughness(SPs)
    # if SPs['ProblemType'] == 'FilmResonance'  : Handle_Geometry_FilmResonance(SPs)

    # OscBndPars = OscBnd.Setup_Boundaries_3D(SPs)
    
    OscBndPars = SPs['OscBndPars']

    print('nx,ny,nz',SPs['nx'],SPs['ny'],SPs['nz'],'n',SPs['n'])
        # nu_for_om = 1./6.
        # SPs['delta'] = SPs['delta0_nm'] / SPs['Dx_nm'] / SPs['n']**0.5  
        # SPs['om']    = 2*nu_for_om / SPs['delta']**2
        
    General.Calc_tauInvBulk_ZBulk(SPs)
    General.Calc_etaabstandel(SPs)
    if SPs['ProblemType'] == 'StiffParticles' : OscBnd.Calc_SphRespPars_3D(SPs,OscBndPars)        
        
    FracVolSph,tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym,rhos = \
        Soft.Set_RelaxPars(SPs)              
        
    RingIn.RingIn(SPs,FracVolSph,OscBndPars,\
        tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym,rhos,Do_Ref = True)
    print('Dfcbyn_Ref' ,np.round(SPs['Dfcbyn_Ref' ],3),\
             'Dfratio_Ref',np.round(SPs['Dfratio_Ref'],3))
    RingIn.RingIn(SPs,FracVolSph,OscBndPars,\
        tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym,rhos,Do_Ref = False)     
