import numpy as np
import datetime 
import os
import configparser

Parameters={}
Parameters['etaabs_tandel'] = ['\u03C1, g/cm\u00B3', 'Radius, nm', 'Truncation', 'Number', '|\u03B7|', 'tan(\u03B4)', '\u03B2\'', '\u03B2\"', 'Coverage'] 
Parameters['from_J'] = ['\u03C1, g/cm\u00B3', 'Radius, nm', 'Truncation', 'Number', 'Jp_FacSphs', 'Jpp_FacSphs', '\u03B2\'', '\u03B2\"', 'Coverage'] 
Parameters['Maxwell'] = ['\u03C1, g/cm\u00B3', 'Radius, nm', 'Truncation', 'Number', '|\u03B7|', '\u03C4', 'Coverage'] 

SPsLoopMap={}
SPsLoopMap['etaabs_tandel'] =['rhoSph','RSph_nm','ySphbyR','nSph','etaabscenSphmPas','tandelcenSph','betap_Sph','betapp_Sph','CovTarget']
SPsLoopMap['from_J']=['rhoSph','RSph_nm','ySphbyR','nSph','Jp_FacSph','Jpp_FacSph','betap_Sph','betapp_Sph','CovTarget']
SPsLoopMap['Maxwell'] =['rhoSph','RSph_nm','ySphbyR','nSph','etaabscenSphmPas','tau','CovTarget']

LoopStartMap={}
LoopStartMap['etaabs_tandel']=[1.35, 10.0, 0.0, 3, 1e4, 0.1, 0.0, 0.0, 0.1]
LoopStartMap['from_J']=[1.35, 10.0, 0.0, 3, 10.0, 10.0, 0.0, 0.0, 0.1]
LoopStartMap['Maxwell'] =[1.35, 10.0, 0.0, 3, 1e4, 6e-7, 0.1]

LoopStepMap={}
LoopStepMap['etaabs_tandel']=['1','1','1','1','1','1','1','1','2']
LoopStepMap['from_J']=['1','1','1','1','1','1','1','1','2']
LoopStepMap['Maxwell']=['1','1','1','1','1','1','2']

Limits={}
Limits['etaabs_tandel']=[-2.0, 0.0, -2.0, 2.0]
Limits['from_J']=[-2.0, 0.0, -1.0, 1.0]
Limits['Maxwell']=[0.0, 0.0, 0.0, 0.0]

Lambda_TRT=[]
Lambda_TRT=['0.0', '1./12.', '1./6.', '3./16.', '1./4.']
Lambda_TRTvals=[0.0, 1./12., 1./6., 3./16., 1./4.]

def Reset_SPs(SPs, GridFactor, ParCount ):             
    GridFactor=2.0
    ParCount = 3

    SPs['Do_from_GUI'] = True
    
    SPs['fname0'] = 'FreqD-LBM-Output'
    SPs['ProblemType']  = 'StiffParticles' #'SoftParticles', 'StiffParticles'        
    SPs['VEPars_Choice']      = 'etaabs_tandel'  #'etaabs_tandel','from_J', 'Maxwell'

    SPs['Par1str'] = 'etaabscenSphmPas'; 
    SPs['Par2str'] = 'tandelcenSph'; 
    SPs['Par3str'] = 'CovTarget'

    if SPs['ProblemType'] in ['SoftParticles','StiffParticles','SFA'] : SPs['dimensions'] = 3
    if SPs['ProblemType'] in ['Roughness']      : SPs['dimensions'] = 2
    if SPs['ProblemType'] in ['FilmResonance']  : SPs['dimensions'] = 1
    
    if SPs['ProblemType'] in ['StiffParticles','Roughness''SFA'] : 
        SPs['Do_OscBnd'] = True
        if SPs['ProblemType'] == 'StiffParticles': 
            SPs['OscBndLocked']   = False
            SPs['OscBndLockedTo'] = 'Zero' 
        if SPs['ProblemType'] == 'Roughness': 
            SPs['OscBndLocked']   = True
            SPs['OscBndLockedTo'] = 'Substrate' 
        if SPs['ProblemType'] == ['SFA']: 
            SPs['OscBndLocked']   = True
            SPs['OscBndLockedTo'] = 'Zero' 
    else : 
        SPs['Do_OscBnd']      = False
        SPs['OscBndLocked']   = False
        SPs['OscBndLockedTo'] = 'Zero' 

    SPs['etaabsBulk']  = 1
    SPs['tandelBulk']  = 1e33
    SPs['f0_SI']      = 5e6
    SPs['Zq_SI']      = 8.8e6
    SPs['delta0_nm']  = 252 * SPs['etaabsBulk']**0.5
    SPs['Gap_P2P']    = 1
    SPs['rhoSph']     = 1.35
    SPs['Gap2TopbyR'] = 1.5

    # only for SPs['ProblemType'] = 'StiffParticles'
    SPs['UpdateMotionFac']  = 0.02; #a number between 0.0005 and 0.1

    SPs['Do_SavePlots']       = False
    SPs['Do_Plot_MotionPars'] = True
    SPs['Do_Plot_RingIns']    = True
    #The number of times ring-ing is plotted
    SPs['PrintIntervalFac']   = 2 # a number between 0.05 and 3; 1 : once per ring-in time 
    SPs['Do_Adapt_Dx2delta']  = False
    
    SPs['RSph_nm']      = 5.0 #RSph_nms[0]
    SPs['Dx_nm'] = SPs['RSph_nm']/GridFactor  # important
    

    # concerns when the ring-in in terminated and whether ring-in data are smoothed
    SPs['TargetSlopeFitResults'] = 0.1 #a number between 0.1 and 10
    SPs['MaxtbytRI']             = 100 #a number between 50 and 200 
    SPs['SigSmoothDfcbynsFac']   = 1e-2 # a number between 1e-3 and 1e-1   
    
    # concerns the collision step
    SPs['Lambda_TRT']            = 1./4. #0. ,     1./12,  1./6.,   3./16., or   1./4. # 0: BGK collision
    
    
    # Roughness
    # SPs['Roughn_VertScale_nm']   = 0
    # SPs['Roughn_Width_nm']       = 10
    # SPs['Roughn_nFourierComps']  = 1
    # SPs['Roughn_Width_nm']       = 10
    
    # Filmresonance
    # SPs['FilmThickness_nm']     = 0


    #parameters which often are looped, "s" at the end of a name indicates an array
    # RSph_nms          = np.array([10])
    # ySphbyRs          = np.array([0.9])
    # CovTargets        = np.linspace(0.30,0.1,8) 
    # ns                = np.array([3,5,7,9])
    # rhoSphs           = np.array([1])  
    # etaabscenSphmPass = np.array([10,30]) # 1e4 corresponds to 1.9 GPa
    # tandelcenSphs     = np.array([0.1,1.,10.])    
    
    SPs['nSph']       = 3  
    SPs['navg']       = 5
      
    SPs['ySphbyR']      = 0.0 #ySphbyRs[0] 
    SPs['CovTarget']    = 0.3 #CovTargets[0] 
    SPs['etaabscenSphmPas'] = 1e4 #etaabscenSphmPass[0]    
    SPs['tandelcenSph'] = 0.1 #tandelcenSphs[0]
    
    # only for SPs['VEPars_Choice' = 'from_J'
    SPs['Jp_FacSph']  = 0.1
    SPs['Jpp_FacSph'] = 0.1

    SPs['betap_Sph']    = 0 # betap_Sphs[0]
    SPs['betapp_Sph']   = 0# betapp_Sphs[0]

    # only for SPs['VEPars_Choice' = 'Maxwell'
    SPs['MaxwellRelaxRate_MHz'] = 1
    SPs['tau'] = 6e-7

   
    
    SPs['ns']           = np.array([7])
    SPs['n']            = SPs['ns'][0]
    SPs['Corrug_nm']    = 0 #Corrug_nms[0] 
        
    LoopStartMap['etaabs_tandel'][(SPsLoopMap['etaabs_tandel'].index('RSph_nm'))] = SPs['RSph_nm']
    LoopStartMap['from_J'][(SPsLoopMap['etaabs_tandel'].index('RSph_nm'))]= SPs['RSph_nm']
    LoopStartMap['Maxwell'][(SPsLoopMap['etaabs_tandel'].index('RSph_nm'))]= SPs['RSph_nm']
    return SPs, GridFactor, ParCount 

def Write_Config_Interface(): 
    config['Main']={}
    config['Main']['WindowHeight']      = str(WindowHeight)
    config['Main']['WindowWidth']       = str(WindowWidth)
    config['Main']['WindowTop']         = str(WindowTop)
    config['Main']['WindowLeft']        = str(WindowLeft)
    with open('FBLMDefaults.ini','w+') as configfile: config.write(configfile)

def Read_Config_Interface():
    global config, WindowHeight, WindowWidth, WindowTop, WindowLeft
    config = configparser.ConfigParser()
    config.read('FBLMDefaults.ini')
    WindowHeight      = int(config.get('Main','WindowHeight'              ,fallback=780))
    WindowWidth       = int(config.get('Main','WindowWidth'               ,fallback=1290))
    WindowTop         = int(config.get('Main','WindowTop'                 ,fallback=15))
    WindowLeft        = int(config.get('Main','WindowLeft'                ,fallback=15))

def ParallelSim(SPs):
    if not os.path.exists(SPs['folder']): os.mkdir(SPs['folder'])
    now = datetime.datetime.now()
    formatted = now.strftime("%Y-%m-%d-%H_%M_%S")    
    fname= formatted + '.flbmsym'
    np.save(SPs['folder'] +'\\'+ fname, SPs) 
    return SPs['folder'] +'\\' + fname + '.npy'

def Make_MotionParsDict_AusParsDict(MotionPars,MotionParTitles,AuxPars,AuxParTitles):
    MotionParsDict = {}
    for i in range(len(MotionPars)):
        MotionParsDict.update({MotionParTitles[i] + '_real': MotionPars[i].real}) 
        MotionParsDict.update({MotionParTitles[i] + '_imag': MotionPars[i].imag}) 
    AuxDict = {}
    for i in range(len(AuxPars)): AuxDict.update({AuxParTitles[i]: AuxPars[i]}) 
    return MotionParsDict,AuxDict

def Set_fname(SPs):
    if not os.path.exists(SPs['folder']): os.mkdir(SPs['folder'])
    now = datetime.datetime.now()
    formatted = now.strftime("%Y-%m-%d-%H_%M_%S")        
    SPs['fname'] = SPs['folder']+'/'+SPs['fname0'] +'_' + formatted +'.txt'

def Create_and_Write_Header(SPs,MotionParsDict,AuxDict):
    header = 'iavg'+'\t'+\
        'i'+SPs['Par1str']+'\t'+\
        'i'+SPs['Par2str']+'\t'+\
        'i'+SPs['Par3str']+'\t'+\
        'i'+'n'           +'\t'+\
            SPs['Par1str']+'\t'+\
            SPs['Par2str']+'\t'+\
            SPs['Par3str']+'\t'+\
            'n'           +'\t'+\
            'Dfbyn'       +'\t'+\
            'DGbyn'       +'\t'+\
            'nNodes'      +'\t'+\
            'steps'       +'\t'+\
            'tbytRI'      +'\t'+\
            'CompTimeMins'+'\t'+\
            'ProblemFlag' +'\t'+\
            'CoverageTrue'+'\t'+\
            'Dfratio_Ref' +'\t'+\
            'Dx_nm'       +'\t'
    for key in list(MotionParsDict.keys()): header += key + '\t'
    f = open(SPs['fname'],'w+'); f.write(header + '\n'); f.close()
    headerAux = 'iavg'+'\t'+\
        'i'+SPs['Par1str']+'\t'+\
        'i'+SPs['Par2str']+'\t'+\
        'i'+SPs['Par3str']+'\t'+\
        'i'+'n'           +'\t'+\
            SPs['Par1str']+'\t'+\
            SPs['Par2str']+'\t'+\
            SPs['Par3str']+'\t'+\
            'n'           +'\t'
    for key in list(AuxDict.keys()): headerAux += key + '\t'
    f = open(SPs['fname'][:-4]+'_Aux.txt','w+'); f.write(headerAux + '\n'); f.close()

def Save(SPs,MotionParsDict,AuxDict):
    if not os.path.isfile(SPs['fname']): 
        print(os.path.basename(SPs['fname']))
        Create_and_Write_Header(SPs,MotionParsDict,AuxDict)
    line = \
        str(SPs['iavg'])+'\t'+str(SPs['iPar1'])+'\t'+str(SPs['iPar2'])+'\t'+\
        str(SPs['iPar3'])+'\t'+str(SPs['iovt' ])+'\t'+\
        str(SPs[SPs['Par1str']])+'\t'+str(SPs[SPs['Par2str']])+'\t'+\
        str(SPs[SPs['Par3str']])+'\t'+str(SPs['n'])+'\t'+\
        str(SPs['Dfcbyn_Extrapol'].real)+'\t'+str(SPs['Dfcbyn_Extrapol'].imag)+'\t' +\
        str(SPs['nNodes'])+'\t'+str(SPs['steps'])+'\t'+\
        str(np.round(SPs['tbytRI'],1))+'\t'+str(SPs['CompTimeMins'])+'\t'+\
        str(SPs['ProblemFlag'])+'\t'+ str(SPs['CoverageTrue'])+'\t'+\
        str(SPs['Dfratio_Ref'])+'\t' + str(SPs['Dx_nm'])+'\t' 
    for key in list(MotionParsDict.keys()): line += str(MotionParsDict[key]) + '\t'
    f = open(SPs['fname'],'a'); f.write(line +'\n'); f.close()   
    
    lineAux=\
        str(SPs['iavg'])+'\t'+str(SPs['iPar1'])+'\t'+str(SPs['iPar2'])+'\t'+\
                        str(SPs['iPar3'])+'\t'+str(SPs['iovt'])+'\t'+\
                        str(SPs[SPs['Par1str']])+'\t'+str(SPs[SPs['Par2str']])+'\t'+\
                        str(SPs[SPs['Par3str']])+'\t'+str(SPs['n'])+'\t'
    for key in list(AuxDict.keys()): lineAux += str(AuxDict[key]) + '\t'
    f = open(SPs['fname'][:-4]+'_Aux.txt','a'); f.write(lineAux +'\n'); f.close()  
    # print('AuxDict',AuxDict)

def Write_Config(SPs): 
    config = configparser.RawConfigParser()
    config.optionxform = str
    config['Main'] = {}
    for key in SPs.keys(): config['Main'][key] = str(SPs[key])
    with open(SPs['fname'][:-4]+'.cfg','w') as configfile: config.write(configfile)