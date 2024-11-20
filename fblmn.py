# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 23:17:45 2024

@author: IlyaReviakine
"""
import sys
import os
import time
import numpy as np
from Libs import Lib_SingleSim as Single_Sim_3D
from Libs import Lib_OscBnd    as OscBnd

print('sys.argv',sys.argv)
SPs={}

if len(sys.argv) == 2:
    try:
        print(len(sys.argv))
        for arg in sys.argv:
            print(arg)
    except:
        print("Command line argument error.")
        input("Press enter to exit.")        
        sys.exit(1)
        
    try:        
        SPs = np.load(sys.argv[1],allow_pickle='TRUE').item()
    except:
        print("In fblmn,py No idea what happened here. File doesn't exist")
        input("Press enter to exit.")
        sys.exit(1)
        
        
    pfolder = SPs['folder'] + '\\tmpplot'
    if not os.path.exists(pfolder): os.mkdir(pfolder)        
    
    print(SPs['nPar1'], SPs['nPar2'], SPs['nPar3'])
    print(SPs['Par1str'], SPs['Par2str'], SPs['Par3str'])
    
    Pars1 = SPs[SPs['Par1str'] + 's']
    Pars2 = SPs[SPs['Par2str'] + 's']
    Pars3 = SPs[SPs['Par3str'] + 's']
    
    GridFactor = SPs['RSph_nm']/SPs['Dx_nm']
    ns=SPs['ns']        
    print('New Simulation Parameters:')
    before = time.time()
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
                    if SPs['ProblemType'] in ['SoftParticles','StiffParticles'] :  Single_Sim_3D.Handle_Geometry_Spheres(SPs)
                    if SPs['ProblemType'] == 'SFA'            : Single_Sim_3D.Handle_Geometry_SFA(SPs)
                    if SPs['ProblemType'] == 'Roughness'      : Single_Sim_3D.Handle_Geometry_Roughness(SPs)
                    if SPs['ProblemType'] == 'FilmResonance'  : Single_Sim_3D.Handle_Geometry_FilmResonance(SPs)

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
                        Single_Sim_3D.SingleSimulation(SPs)
                        
                        
    after = time.time()
    print('Done!')
    sp, sf = os.path.split(sys.argv[1])
    stponame = pfolder +'\\' + sf[:-12]+'_STOP_'
    np.save(stponame, SPs, allow_pickle=True)
    print(sys.argv[1])
    print(str(after - before)+' Press enter to exit.')
    # try:
    #     if os.path.isfile(sys.argv[1]): os.remove(sys.argv[1])
    # except:
    #     print('error')
    input()
        

    

