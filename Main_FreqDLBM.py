import numpy as np
import sys
import os
Base_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(Base_Path);

from Libs import Lib_SingleSim as Single_Sim
from Libs import Lib_IO        as IO
from Libs import Lib_OscBnd    as OscBnd

SPs = {}
SPs['ProblemType']  = 'StiffParticles' #'SoftParticles', 'StiffParticles',\
SPs['Dx_nm']     = 1  # important

#parameters which often are looped, "s" at the end of a name indicates an array
RSph_nms          = np.array([5])
ySphbyRs          = np.array([0.9])
CovTargets        = np.array([0.3,0.1])
ns                = np.array([7])
rhoSphs           = np.array([1])
etaabscenSphmPass = np.array([1e4]) # 1e4 corresponds to 1.9 GPa
tandelcenSphs     = np.array([0.1])

Jp_FacSphs   = np.array([10])
Jpp_FacSphs   = np.array([10])

Pars1 = etaabscenSphmPass;  SPs['Par1str'] = 'etaabscenSphmPas';
Pars2 = tandelcenSphs;      SPs['Par2str'] = 'tandelcenSph';
Pars3 = CovTargets;         SPs['Par3str'] = 'CovTarget'

SPs['folder'] = 'test'
SPs['fname0'] = 'test'
SPs['Do_from_GUI'] = False

if SPs['ProblemType'] in ['SoftParticles','StiffParticles','SFA','Roughness_3D'] : 
    SPs['dimensions'] = 3
if SPs['ProblemType'] in ['Roughness_2D']   : SPs['dimensions'] = 2
if SPs['ProblemType'] in ['FilmResonance']  : SPs['dimensions'] = 1
#'Roughness_2D','Roughness_3D','SFA','FilmResonance' not debugged yet

if SPs['ProblemType'] in ['StiffParticles','Roughness_2D','Roughness_3D','SFA'] :
    SPs['Do_OscBnd'] = True
    if SPs['ProblemType'] == 'StiffParticles':
        SPs['OscBndLocked']   = False
        SPs['OscBndLockedTo'] = 'Zero'
    if SPs['ProblemType'] in ['Roughness_2D','Roughness_3D']:
        SPs['OscBndLocked']   = True
        SPs['OscBndLockedTo'] = 'Substrate'
    if SPs['ProblemType'] == ['SFA']:
        SPs['OscBndLocked']   = True
        SPs['OscBndLockedTo'] = 'Zero'
else :
    SPs['Do_OscBnd']      = False
    SPs['OscBndLocked']   = False
    SPs['OscBndLockedTo'] = 'Zero'

SPs['VEPars_Choice']     = 'etaabs_tandel'  #'etaabs_tandel','from_J', 'Maxwell'

SPs['etaabsBulk'] = 1
SPs['tandelBulk'] = 1e33
SPs['f0_SI']      = 5e6
SPs['Zq_SI']      = 8.8e6
SPs['delta0_nm']  = 252 * SPs['etaabsBulk']**0.5
SPs['Gap_P2P']    = 1
SPs['rhoSph']     = 1
SPs['Gap2TopbyR'] = 1.5
SPs['betap_Sph']  = 0
SPs['betapp_Sph'] = 0
# only for SPs['VEPars_Choice' = 'from_J'
SPs['Jp_FacSph']  = Jp_FacSphs[0]
SPs['Jpp_FacSph'] = Jpp_FacSphs[0]
# only for SPs['VEPars_Choice' = 'Maxwell'
SPs['MaxwellRelaxRate_MHz'] = 1

# only for SPs['ProblemType'] = 'StiffParticles'
SPs['UpdateMotionFac']  = 0.02;

SPs['Do_SavePlots']       = False
SPs['Do_Plot_MotionPars'] = True
SPs['Do_Plot_RingIns']    = True
SPs['PrintIntervalFac']   = 1 # 1 : once per ring-in time

SPs['nSph']       = 3
SPs['navg']       = 1

# concerns when the ring-in in terminated and whether ring-in data are smoothed
SPs['TargetSlopeFitResults'] = 20
SPs['MaxtbytRI']             = 100
SPs['SigSmoothDfcbynsFac']   = 1e-2

# concerns the collision step
SPs['Lambda_TRT']            = 1./4. # 0: BGK collision
SPs['Do_UseQuadraticTerm']   = True
SPs['Do_Allow_rhoUneq1']     = False


# Roughness
SPs['Roughn_VertScale_nm']   = 3
SPs['Roughn_HoriScale_nm']   = 5
SPs['Roughn_Width_nm']       = 100
SPs['Single_Wave']           = True

# Filmresonance
SPs['FilmThickness_nm']     = 5

#  assing values to simulation parameters, which are not looped
SPs['RSph_nm']          = RSph_nms[0]
SPs['ySphbyR']          = ySphbyRs[0]
SPs['CovTarget']        = CovTargets[0]
SPs['n']                = ns[0]
SPs['rhoSph']           = rhoSphs[0]
SPs['etaabscenSphmPas'] = etaabscenSphmPass[0]
SPs['tandelcenSph']     = tandelcenSphs[0]

# for documentation in cfg-file
SPs['ns']                 = ns
SPs[SPs['Par1str'] + 's'] = Pars1
SPs[SPs['Par2str'] + 's'] = Pars2
SPs[SPs['Par3str'] + 's'] = Pars3
GridFactor = SPs['RSph_nm']/SPs['Dx_nm']

IO.Set_fname(SPs)
SPs['nPar1']=len(Pars1);SPs['nPar2']=len(Pars2);SPs['nPar3']=len(Pars3);SPs['novt']=len(ns);
for SPs['iavg'] in range(SPs['navg']):        
    for SPs['iPar1'],SPs[SPs['Par1str']] in enumerate(Pars1):
        for SPs['iPar2'],SPs[SPs['Par2str']] in enumerate(Pars2):
            for SPs['iPar3'],SPs[SPs['Par3str']] in enumerate(Pars3):
                if SPs['Par1str']=='RSph_nm':
                    SPs['Dx_nm'] =  SPs[SPs['Par1str']]/GridFactor
                elif SPs['Par2str']=='RSph_nm':
                    SPs['Dx_nm']     =  SPs[SPs['Par2str']]/GridFactor
                elif SPs['Par3str']=='RSph_nm':
                    SPs['Dx_nm'] =  SPs[SPs['Par3str']]/GridFactor
                else: 
                    SPs['Dx_nm'] =  SPs['RSph_nm']/GridFactor
                if SPs['ProblemType'] in ['SoftParticles','StiffParticles'] :  Single_Sim.Handle_Geometry_Spheres(SPs)
                if SPs['ProblemType'] == 'SFA'            : Single_Sim.Handle_Geometry_SFA(SPs)
                if SPs['ProblemType'] == 'Roughness'      : Single_Sim.Handle_Geometry_Roughness(SPs)
                if SPs['ProblemType'] == 'FilmResonance'  : Single_Sim.Handle_Geometry_FilmResonance(SPs)

                SPs['OscBndPars'] = OscBnd.Setup_Boundaries_3D(SPs)
                nu_for_om = 1./6.
                SPs['delta'] = SPs['delta0_nm'] / SPs['Dx_nm'] / SPs['n']**0.5  
                SPs['om']    = 2*nu_for_om / SPs['delta']**2
                for SPs['iovt'],SPs['n'] in enumerate(ns):
                    print('******************')
                    print('******************')
                    print('Dx_nm = ', SPs['Dx_nm'], 'n = ', SPs['n'], 'iavg = ', SPs['iavg']+1, SPs['Par1str'], ' = ', SPs[SPs['Par1str']], SPs['Par2str'], ' = ', SPs[SPs['Par2str']], SPs['Par3str'], ' = ', SPs[SPs['Par3str']])
                    print('******************')
                    print('******************')                     
                    Single_Sim.SingleSimulation(SPs)
