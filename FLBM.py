# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:59:49 2024

@author: IlyaReviakine
"""
import os
import sys
import time
import queue
import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
import webbrowser
import multiprocessing as mpr
from subprocess import Popen, CREATE_NEW_CONSOLE
from tkinter import ttk
from tkinter.filedialog import askdirectory
from tkinter.filedialog import asksaveasfilename
from tkinter.filedialog import askopenfilename
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tooltip import ToolTip  

import warnings
warnings.filterwarnings("ignore")

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_PATH);

from Libs import Lib_IO           as IO
from Libs import Lib_General       as Gen
from Libs import Lib_Soft as Geo_SoftPt_3D
from Libs import Lib_OscBnd as Geo_OscBnd_3D
from Libs import Lib_SingleSim as Single_Sim_3D
from Libs import Lib_Plots_from_GUI as Plots_3D
from Libs  import Lib_DResults

global Version
global bgrc, fgrc, wwidth, GridFactor, ParCount
global RSButton, SimType, VEParsFormatValues, ProblemTypeValues
global fig1, fig2, fig3, fig4, fig5
global canvas1, canvas2, canvas3, canvas4, canvas5
global menubar, APShow, LimShow
global DisplayPar, llines, DisplayPar


FBLMVersion= 'Version 1.0 2024.11.14'
SimType='Regular'
llines = (mpr.cpu_count())
VEParsFormatValues=['|\u03B7|, tan(\u03B4)','J\', J\"','Maxwell'] #['|\u03B7|, tan(\u03B4)','J\', J\"','Maxwell']
ProblemTypeValues=['Stiff Particles', 'Soft Particles']
bgrc = 'white'
fgrc = 'black'
APShow = False
LimShow = False

def CheckProblemType(SPs):
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
    return SPs        

def Initialize(SPs, flag):
    global wwidth, GridFactor, ParCount
    wwidth = 10
    SPs = {}
    GridFactor=0.0
    ParCount = 0.0
    SPs['PlotSimProgress'] = 0
    if flag=='reset':
        SPs['folder'] = os.path.abspath(__file__)[:-3]
        fname = os.path.dirname(os.path.abspath(__file__)) + '\\Defaults.npy'
        if os.path.isfile(fname): os.remove(fname)
        SPs, GridFactor, ParCount = IO.Reset_SPs(SPs, GridFactor, ParCount)
    elif flag=='default':
        SPs['folder'] = os.path.abspath(__file__)[:-3]
        fname = os.path.dirname(os.path.abspath(__file__)) + '\\Defaults.npy'
        if os.path.isfile(fname): 
            tempdata={}
            try:
                tempdata = np.load(fname, allow_pickle='TRUE').item()
                SPs = tempdata['SPs']
                if not os.path.exists(SPs['folder']): SPs['folder'] = os.path.abspath(__file__)[:-3]
                IO.Limits[SPs['VEPars_Choice']]=tempdata['Limits']
                IO.SPsLoopMap[SPs['VEPars_Choice']] = tempdata['SPsLoopMap'] 
                IO.LoopStartMap[SPs['VEPars_Choice']] = tempdata['LoopStartMap'] 
                IO.LoopStepMap[SPs['VEPars_Choice']] = tempdata['LoopStepMap']
                GridFactor = tempdata['GridFactor']  
            except:    
                print('Something didn\'t work')
                SPs, GridFactor, ParCount = IO.Reset_SPs(SPs, GridFactor, ParCount)
        else: 
            print('No Defaults file, reseting.')
            SPs, GridFactor, ParCount = IO.Reset_SPs(SPs, GridFactor, ParCount)
    SPs = CheckProblemType(SPs)
    return SPs
    
def UpdateSpherePlot(cntr, SPs):
    global fig1, canvas1   
    pfolder = SPs['folder'] + '\\tmpplot'
    if not os.path.exists(pfolder): os.mkdir(pfolder) 

    plt.rcParams["figure.figsize"] = (2.5,2.0); 
    plt.rcParams['lines.markersize'] = 2.0
    plt.rcParams["font.size"] = 7;        

    fig1 = plt.figure()                   
    canvas1 = FigureCanvasTkAgg(fig1,cntr); 
    canvas1.get_tk_widget().grid(row=1,column=20, rowspan = 16, columnspan = 2, sticky='N')
    canvas1.draw()
    
    SPs['Dx_nm']     =  SPs['RSph_nm']/GridFactor    
    if SPs['Do_Adapt_Dx2delta'] : SPs['Dx_nm'] = SPs['Dx_nm'] / SPs['n']**0.5
    else                        : SPs['Dx_nm'] = SPs['Dx_nm']
        
    nu_for_om = 1./6.
    SPs['delta'] = SPs['delta0_nm'] / SPs['Dx_nm'] / SPs['n']**0.5  
    SPs['om']    = 2*nu_for_om/SPs['delta']**2

    Gen.Calc_tauInvBulk_ZBulk(SPs)
    Gen.Calc_etaabstandel(SPs)
    if SPs['ProblemType'] in ['SoftParticles' ,'StiffParticles']: #Other flags are now ignored.
        Single_Sim_3D.Handle_Geometry_Spheres(SPs)
        FracVolSph = Geo_SoftPt_3D.Calc_FracVolSph_3D(SPs)
        FracVolSph,tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym,rhos = Geo_SoftPt_3D.Set_RelaxPars(SPs)   
        OscBndPars = Geo_OscBnd_3D.Setup_Boundaries_3D(SPs)                 
        tempData={}
        tempData['xLs'] = OscBndPars['xLs']
        tempData['yLs'] = OscBndPars['yLs'] 
        tempData['zLs'] = OscBndPars['zLs']
        tempData['color1'] = OscBndPars['yLs']*SPs['Dx_nm']
        tempData['title1'] = ' '
        fname = 'temp_Plot_LinkProps_3D'
        np.save(fname, tempData) 
        fname = os.path.realpath(fname)+'.npy'
        Plots_3D.Plot_LinkProps_3D(fname, fig1)
        canvas1.draw()
        plt.close("all")
        if os.path.isfile(fname): os.remove(fname)

#=================================================================================================
# Root Window houskeeping functions
def onClose():
    try: Lib_DResults.DRWinroot.destroy() 
    except: pass

    try: Root.destroy()
    except: pass

def on_mousewheel(event):
    MainCanvas.yview_scroll(int(-1*(event.delta/120)), "units")

def ScrollFunction(event):
    w = Root.winfo_width()
    h = Root.winfo_height()
    MainCanvas.configure(scrollregion=MainCanvas.bbox("all"),width=w-25,height=h-25)

    IO.WindowWidth  = w
    IO.WindowHeight = h
    IO.WindowTop    = Root.winfo_y()
    IO.WindowLeft   = Root.winfo_x()
    IO.Write_Config_Interface()
    
def onConfigure(event):
    def update_size():
        w = Root.winfo_width()
        h = Root.winfo_height()
        MainCanvas.configure(scrollregion=MainCanvas.bbox("all"),width=w-25,height=h-25)

    w = Root.winfo_width()
    h = Root.winfo_height()
    if w > 1500: 
        w=1500
        Root.geometry(str(int(w)) +'x' + str(int(h)) + '+15+15')
    if event.widget == Root:     
            if getattr(Root, "_after_id", None): 
                Root.after_cancel(Root._after_id)
            Root._after_id=Root.after(100, update_size)
    IO.WindowWidth  = w
    IO.WindowHeight = h+50
    IO.WindowTop    = Root.winfo_y()
    IO.WindowLeft   = Root.winfo_x()
    IO.Write_Config_Interface()

#==============================================================================

def MakePlotSim1Frame(cntr, SPs):
    global fig2, fig3, fig4, fig5
    global canvas2, canvas3, canvas4, canvas5

    for widget in cntr.winfo_children(): widget.destroy()
    tk.Label(cntr,text=' ',anchor='n').grid(row=0 ,column=0,sticky='W')  
    
    if SPs['ProblemType'] =='StiffParticles':
        plt.rcParams['figure.figsize'] = (4,4)
        plt.rcParams['font.size'] = 8
        fig2 = plt.figure()
        canvas2 = FigureCanvasTkAgg(fig2,cntr); 
        canvas2.get_tk_widget().grid(row=0,column=0, columnspan = 2, sticky='N')            
        Plots_3D.Plot_RI(' ', fig2)
        canvas2.draw() 
        plt.close('all')    
    
        plt.rcParams['figure.figsize'] = (5,4); plt.rcParams['font.size'] = 8
        plt.rcParams['lines.markersize'] = 3
        fig3 = plt.figure()
        canvas3 = FigureCanvasTkAgg(fig3,cntr); 
        canvas3.get_tk_widget().grid(row=0,column=2, columnspan = 2, sticky='N')            
        Plots_3D.Plot_MotionPars_RI(' ', fig3)
        canvas3.draw() 
        plt.close('all')    
        
        plt.rcParams["figure.figsize"] = (9,1.8)
        plt.rcParams["font.size"] = 8;
        fig4 = plt.figure()
        canvas4 = FigureCanvasTkAgg(fig4,cntr); 
        canvas4.get_tk_widget().grid(row=1,column=0, columnspan = 4, sticky='N')            
        Plots_3D.Plot_Fields_Horizontal(' ', fig4)
        canvas4.draw()
        plt.close('all')    

        plt.rcParams["figure.figsize"] = (9,1.8)
        plt.rcParams["font.size"] = 8;
        fig5 = plt.figure()
        Plots_3D.Plot_Fields_Vertical(' ', fig5)
        canvas5 = FigureCanvasTkAgg(fig5,cntr); 
        canvas5.get_tk_widget().grid(row=2,column=0, columnspan = 4, sticky='N')            
        canvas5.draw()
        plt.close('all')    
        
    elif SPs['ProblemType'] =='SoftParticles':
        plt.rcParams['figure.figsize'] = (6,3)
        plt.rcParams['font.size'] = 8
        fig2 = plt.figure()
        canvas2 = FigureCanvasTkAgg(fig2,cntr); 
        canvas2.get_tk_widget().grid(row=0,column=0, columnspan = 2, sticky='N')            
        Plots_3D.Plot_RI(' ', fig2)
        canvas2.draw()
        plt.close('all')    
        
        plt.rcParams["figure.figsize"] = (6,1.8)
        plt.rcParams["font.size"] = 8;
        fig4 = plt.figure()
        canvas4 = FigureCanvasTkAgg(fig4,cntr); 
        canvas4.get_tk_widget().grid(row=1,column=0, columnspan = 4, sticky='N')            
        Plots_3D.Plot_Fields_Horizontal(' ', fig4)
        canvas4.draw() 
        plt.close('all')    
    
        plt.rcParams["figure.figsize"] = (6,1.8)
        plt.rcParams["font.size"] = 8;
        fig5 = plt.figure()
        Plots_3D.Plot_Fields_Vertical(' ', fig5)
        canvas5 = FigureCanvasTkAgg(fig5,cntr); 
        canvas5.get_tk_widget().grid(row=2,column=0, columnspan = 4, sticky='N')            
        canvas5.draw()
        plt.close('all')    
    
def MakeModelParamFrame(cntr,SPs):

    global RSButton, SimType, bgrc, fgrc, APShow, LimShow
    for widget in cntr.winfo_children(): widget.destroy()    
    
    def ProblemSelect(SPs):
        global bgrc, fgrc
        SPs['ProblemType']=ProblemTypeCbox.get().replace(" ", "")
        if  SPs['ProblemType']=='SoftParticles': 
            bgrc = 'white' 
            fgrc = 'blue'
            SimTrackFrame.config(width=1160)
        elif SPs['ProblemType']=='StiffParticles': 
            bgrc = 'white' 
            fgrc = 'black'
            SimTrackFrame.config(width=1450)
        CheckProblemType(SPs)
        MakeSphereParamFrame(SphereParamFrame,SPs)
        MakeModelParamFrame(cntr, SPs)
        MakePlotSim1Frame(PlotSim1Frame, SPs)
        UpdateSpherePlot(SphereParamFrame, SPs)
    
    def SimTypeSelect():
        global SimType
        SimType=SimTypeOption.get()
        MakeModelParamFrame(cntr, SPs)

    def FormatUpdate(SPs):
        if FormatCbox.get() == FormatCbox['values'][0]: 
            Root.nametowidget('.!menu').entryconfig(2, state="normal")
            SPs['VEPars_Choice'] = 'etaabs_tandel' 
            if SPs['Par1str'] == 'Jp_FacSph': SPs['Par1str'] = 'etaabscenSphmPas'
            if SPs['Par2str'] == 'Jp_FacSph': SPs['Par2str'] = 'etaabscenSphmPas'
            if SPs['Par3str'] == 'Jp_FacSph': SPs['Par3str'] = 'etaabscenSphmPas'
            if SPs['Par1str'] == 'Jpp_FacSph': SPs['Par1str'] = 'tandelcenSph'
            if SPs['Par2str'] == 'Jpp_FacSph': SPs['Par2str'] = 'tandelcenSph'
            if SPs['Par3str'] == 'Jpp_FacSph': SPs['Par3str'] = 'tandelcenSph'                       
            if SPs['Par1str'] == 'tau': SPs['Par1str'] = 'tandelcenSph'
            if SPs['Par2str'] == 'tau': SPs['Par2str'] = 'tandelcenSph'
            if SPs['Par3str'] == 'tau': SPs['Par3str'] = 'tandelcenSph'                                   
        if FormatCbox.get()== FormatCbox['values'][1]:
            Root.nametowidget('.!menu').entryconfig(2, state="normal")
            SPs['VEPars_Choice'] = 'from_J' 
            if SPs['Par1str'] == 'etaabscenSphmPas': SPs['Par1str'] = 'Jp_FacSph'
            if SPs['Par2str'] == 'etaabscenSphmPas': SPs['Par2str'] = 'Jp_FacSph'
            if SPs['Par3str'] == 'etaabscenSphmPas': SPs['Par3str'] = 'Jp_FacSph'
            if SPs['Par1str'] == 'tandelcenSph': SPs['Par1str'] = 'Jpp_FacSph'
            if SPs['Par2str'] == 'tandelcenSph' : SPs['Par2str'] = 'Jpp_FacSph'
            if SPs['Par3str'] == 'tandelcenSph': SPs['Par3str'] ='Jpp_FacSph'
            if SPs['Par1str'] == 'tau': SPs['Par1str'] = 'Jpp_FacSph'
            if SPs['Par2str'] == 'tau': SPs['Par2str'] = 'Jpp_FacSph'
            if SPs['Par3str'] == 'tau': SPs['Par3str'] = 'Jpp_FacSph'                       
        if FormatCbox.get()== FormatCbox['values'][2]:   
            Root.nametowidget('.!menu').entryconfig(2, state="disabled")
            SPs['VEPars_Choice'] = 'Maxwell'
            if SPs['Par1str'] == 'tandelcenSph': SPs['Par1str'] = 'tau'
            if SPs['Par2str'] == 'tandelcenSph' : SPs['Par2str'] = 'tau'
            if SPs['Par3str'] == 'tandelcenSph': SPs['Par3str'] ='tau'
            if SPs['Par1str'] == 'Jpp_FacSph': SPs['Par1str'] = 'tau'
            if SPs['Par2str'] == 'Jpp_FacSph': SPs['Par2str'] = 'tau'
            if SPs['Par3str'] == 'Jpp_FacSph': SPs['Par3str'] = 'tau'
            
        InterfaceUpdate(SPs)
        
        
    def Enter_MPars(event,Selector):
        global GridFactor
        if Selector  == 'GridR': 
            string = GridR_Enter.get(); string.strip(); 
            try:
                GridFactor = np.float64(string)
            except:
                GridFactor = 5.0
            if GridFactor < 1.5:  GridFactor = 5.0
            SPs['Dx_nm'] =  SPs['RSph_nm']/GridFactor
            MakeModelParamFrame(cntr, SPs)
            UpdateSpherePlot(SphereParamFrame, SPs)
            
    
    def OvertoneSelect(SPs, CheckVar):
        global VEParsFormatValues, ProblemTypetValues
        SPs['ns']=[]
        for i in range(4): 
            if CheckVar[i].get()==1: SPs['ns']=np.append(SPs['ns'],(2*(i+1)+1))
        #check if at least one overtone is selected
        if len(SPs['ns']) == 0: 
            SPs['ns']=np.append(SPs['ns'],7)
            MakeModelParamFrame(cntr, SPs)  

    
    def Enter_Av(SPs):
        string = Averages_Enter.get(); string.strip();
        try:
            SPs['navg']=np.int32(string)
            if SPs['navg'] < 1: SPs['navg']=1
        except:
            SPs['navg']=1
        if  SPs['navg'] < 0: SPs['navg'] = 1 
        MakeModelParamFrame(cntr, SPs)

    def AdvParUpdate(SPs):
        def string_convert(string):
            s = -200
            try:
                s = np.float64(string)
            except:
                s=-200.00
            return s
        
        SPs['UpdateMotionFac'] = string_convert((UpdateMotionFac_enter.get().strip()))
        SPs['MaxtbytRI'] = string_convert((MaxtbytRI_enter.get().strip()))       
        SPs['PrintIntervalFac']   = string_convert((PrintIntervalFac_enter.get().strip()))       
        SPs['TargetSlopeFitResults'] = string_convert((TargetSlopeFitResults_enter.get().strip()))            
        SPs['SigSmoothDfcbynsFac'] = string_convert((SigSmoothDfcbynsFacenter.get().strip()))
        SPs['Gap_P2P'] = int(string_convert((UpdateGapP2P_enter.get().strip())))
        
        IO.Limits[SPs['VEPars_Choice']][0] = string_convert((betap_Sph_min_enter.get().strip()))
        IO.Limits[SPs['VEPars_Choice']][1] = string_convert((betap_Sph_max_enter.get().strip()))
        IO.Limits[SPs['VEPars_Choice']][2] = string_convert((betapp_Sph_min_enter.get().strip()))
        IO.Limits[SPs['VEPars_Choice']][3] = string_convert((betapp_Sph_max_enter.get().strip()))

        if SPs['UpdateMotionFac'] < 0.0005 :SPs['UpdateMotionFac'] = 0.02 #a number between 0.0005 and 0.1
        if SPs['MaxtbytRI']<=0 : SPs['MaxtbytRI'] = 100 #a number between 50 and 200
        if SPs['PrintIntervalFac']<= 0.05:SPs['PrintIntervalFac'] = 2 #a number between 0.05 and 3
        if SPs['TargetSlopeFitResults'] <= 0 : SPs['TargetSlopeFitResults'] = 0.1 #a number between 0.1 and 10
        if SPs['SigSmoothDfcbynsFac']<= 0:SPs['SigSmoothDfcbynsFac'] = 1e-2 # a number between 1e-3 and 1e-1       
        if SPs['Gap_P2P'] <= 0: SPs['Gap_P2P']=1
        
        if IO.Limits[SPs['VEPars_Choice']][0] == -200:
            if SPs['VEPars_Choice'] == 'etaabs_tandel':
                IO.Limits['etaabs_tandel'][0]=-2.0
            elif SPs['VEPars_Choice'] == 'from_J':
                IO.Limits['from_J'][0]=-2.0
            elif SPs['VEPars_Choice'] == 'Maxwell':
                IO.Limits['Maxwell'][0] = 0.0                    
        if IO.Limits[SPs['VEPars_Choice']][1] == -200:            
            if SPs['VEPars_Choice'] == 'etaabs_tandel':
                IO.Limits['etaabs_tandel'][0]=-0.0
            elif SPs['VEPars_Choice'] == 'from_J':
                IO.Limits['from_J'][1]=0.0
            elif SPs['VEPars_Choice'] == 'Maxwell':            
                IO.Limits['Maxwell'][1] = 0.0                    
        if IO.Limits[SPs['VEPars_Choice']][2] == -200:                        
            if SPs['VEPars_Choice'] == 'etaabs_tandel':
                IO.Limits['etaabs_tandel'][0]=-2.0
            elif SPs['VEPars_Choice'] == 'from_J':
                IO.Limits['from_J'][2]=-1.0
            elif SPs['VEPars_Choice'] == 'Maxwell':
                IO.Limits['Maxwell'][2] = 0.0                                    
        if IO.Limits[SPs['VEPars_Choice']][3] == -200:                        
            if SPs['VEPars_Choice'] == 'etaabs_tandel':
                IO.Limits['etaabs_tandel'][0]=2.0
            elif SPs['VEPars_Choice'] == 'from_J':
                IO.Limits['from_J'][3]=1.0        
            elif SPs['VEPars_Choice'] == 'Maxwell':            
                IO.Limits['Maxwell'][3] = 0.0                    
        MakeModelParamFrame(ModelParamFrame,SPs)
        
    def Lambda_TRT_Selection(SPs):
       s = Lambda_TRT_cbox.get()
       SPs['Lambda_TRT']=IO.Lambda_TRTvals[IO.Lambda_TRT.index(s)]
       
    def PlotProgress(SPs, s):
        if s.get() == 1: SPs['PlotSimProgress'] = True
        else: SPs['PlotSimProgress'] = False
    
    
    #--------------------------------------------------------------------------------------------------------------
    tk.Label(cntr,text='Problem',anchor='n').grid(row=0 ,column=0,columnspan=2,sticky='N')
    tk.Label(cntr,text='Grid factor, Rsph/X',anchor='n').grid(row=0 ,column=2,columnspan=2,sticky='N')
    tk.Label(cntr,text='\u0394x, nm',anchor='n').grid(row=0 ,column=4,columnspan=2,sticky='N')
    tk.Label(cntr,text='        Averages      '             ,anchor='n').grid(row=0 ,column=6,columnspan=4,sticky='N')    
    tk.Label(cntr,text='Overtones',anchor='n').grid(row=0 ,column=10,columnspan=4,sticky='N')
    
    tk.Label(cntr,text='Format',anchor='n').grid(row=2,column=0,columnspan=2,sticky='N')            
    tk.Label(cntr,text='Sim Type',anchor='n').grid(row=2 ,column=3, columnspan=3,sticky='N')    
    SimProgPlot=tk.Label(cntr,text='Progress Plots',anchor='n')
    SimProgPlot.grid_remove()
    
    ProgressPlotOptionVar = tk.IntVar(value = 0)
    ProgressPlotOption = tk.Checkbutton(cntr, text='',variable = ProgressPlotOptionVar, onvalue = 1, offvalue = 0, command=lambda : PlotProgress(SPs, ProgressPlotOptionVar))
    ProgressPlotOption.grid_remove()
    
    if SimType=='Parallel':
        SimProgPlot.grid(row=2 ,column=6, columnspan=1,sticky='N')    
        ProgressPlotOption.grid(row=3 ,column=6, columnspan=1,sticky='N') 
        ffgr='red'
        btntxt='Run Parallel \n Simulation'
    else:
        SimProgPlot.grid_remove()
        ProgressPlotOption.grid_remove()
        ffgr='black'
        btntxt='Run Simulation'
    
    RSButton=tk.Button(cntr,command=lambda : Run_Sim(SPs),foreground=ffgr,text=btntxt,justify='center', padx=25, pady=10)
    RSButton.grid(row=2,column=10,columnspan=4,rowspan=2,sticky='N')
    
    temp = tk.StringVar(cntr, SPs['ProblemType'])
    ProblemType = temp.get()
    ProblemTypeCbox=ttk.Combobox(cntr, width=12, state='readonly', justify="center", values = ProblemTypeValues, textvariable=ProblemType)
    ProblemTypeCbox.grid(row=1,column=0,columnspan=2,sticky='S')
    ProblemTypeCbox.bind('<<ComboboxSelected>>', lambda _ : ProblemSelect(SPs))
    if SPs['ProblemType'] == 'StiffParticles': i = 0 
    elif SPs['ProblemType'] == 'SoftParticles': i = 1   
    else: i=0
    ProblemTypeCbox.current(i)
    
    temp = tk.StringVar(cntr, VEParsFormatValues[0])
    VEParsFormat = temp.get()
    FormatCbox=ttk.Combobox(cntr, width=12, state='readonly', justify="center", values = VEParsFormatValues, textvariable=VEParsFormat)
    if SPs['VEPars_Choice']=='etaabs_tandel':i=0
    elif SPs['VEPars_Choice']=='from_J':i=1
    elif SPs['VEPars_Choice']=='Maxwell':i=2
    else: i = 0    
    FormatCbox.current(i)
    FormatCbox.grid(row=3,column=0,columnspan=2,sticky='N')
    FormatCbox.bind('<<ComboboxSelected>>', lambda _ : FormatUpdate(SPs))

    SPs['Dx_nm']=SPs['RSph_nm']/GridFactor
    GridR=tk.StringVar(cntr, GridFactor)
       
    GridR_Enter = tk.Entry(cntr,width=6, textvariable=GridR)
    GridR_Enter.grid(row=1,column=2,columnspan=2,sticky='S')
    GridR_Enter.bind('<Return>',lambda event: Enter_MPars(event,'GridR'))
    GridR_Enter.configure({"justify" : "center"})
    
    if SPs['Dx_nm'] >= SPs['RSph_nm']:  GridR_Enter.configure({"foreground": "red"}) 
    else: GridR_Enter.configure({"foreground": "black"})
    
    tk.Label(cntr,text=str("{:.2f}".format(SPs['Dx_nm'])),anchor='n').grid(row=1 ,column=4,columnspan=2,sticky='S')

    AveragesValue=tk.IntVar(cntr, SPs['navg'])
    Averages_Enter=tk.Entry(cntr,width=wwidth-4, background=bgrc, foreground=fgrc, state='normal', textvariable=AveragesValue)
    Averages_Enter.grid(row=1,column=6,columnspan=4,sticky='S')
    Averages_Enter.configure({"justify" : "center"})
    Averages_Enter.bind('<Return>',lambda event: Enter_Av(SPs))    
    
    CheckVar=[]
    for i in range(4):
        CheckVar.append(tk.IntVar(value=0))
    for i in range(4):                        
        if i < len(SPs['ns']): 
            CheckVar[ int((SPs['ns'][i]-1.0)/2.0-1.0) ].set(1)
        
    tk.Checkbutton(cntr, text='3',variable = CheckVar[0],onvalue = 1, offvalue = 0, command=lambda : OvertoneSelect(SPs,CheckVar)).grid(row=1,column=10,columnspan=1,sticky='N')
    tk.Checkbutton(cntr, text='5',variable = CheckVar[1],onvalue = 1, offvalue = 0, command=lambda : OvertoneSelect(SPs,CheckVar)).grid(row=1,column=11,columnspan=1,sticky='N')
    tk.Checkbutton(cntr, text='7',variable = CheckVar[2],onvalue = 1, offvalue = 0, command=lambda : OvertoneSelect(SPs,CheckVar)).grid(row=1,column=12,columnspan=1,sticky='N')
    tk.Checkbutton(cntr, text='9',variable = CheckVar[3],onvalue = 1, offvalue = 0, command=lambda : OvertoneSelect(SPs,CheckVar)).grid(row=1,column=13,columnspan=1,sticky='N')
        
    SimTypeOption = tk.StringVar(cntr, SimType)
    tk.Radiobutton(cntr, text='Regular', value='Regular', variable=SimTypeOption, command=lambda:SimTypeSelect()).grid(row=3,column=3,columnspan=1,sticky='N')
    tk.Radiobutton(cntr, text='Parallel', value='Parallel', variable=SimTypeOption, foreground='red',command=lambda:SimTypeSelect()).grid(row=3,column=5,columnspan=1,sticky='N')

    AdvancedSimPars=tk.Label(ModelParamFrame,text='Advanced Simulation Parameters', fg = 'maroon', anchor='n')    
    UpdateMotionFacLabel=tk.Label(ModelParamFrame,text='UpdateMotionFac',anchor='n')  
    ToolTip(widget = UpdateMotionFacLabel, text = "A number between 0.0005 and 0.1 ")
    TargetSlopeFitResultsLabel = tk.Label(ModelParamFrame,text='TargetSlopeFitResults',anchor='n')        
    ToolTip(widget = TargetSlopeFitResultsLabel, text = "A number between 0.1 and 10")
    PrintIntervalFacLabel = tk.Label(ModelParamFrame,text='PrintIntervalFac',anchor='n')  
    ToolTip(widget = PrintIntervalFacLabel, text = "A number between 0.05 and 3")
    MaxtbytRILabel=tk.Label(ModelParamFrame,text='max(t/t_RI)',anchor='n')
    ToolTip(widget = MaxtbytRILabel, text = "A number between 50 and 200")
    Lambda_TRTLabel=tk.Label(ModelParamFrame,text='\u039B_TRT',anchor='n')
    ToolTip(widget = Lambda_TRTLabel, text = "0., 1./12, 1./6.,  3./16., or  1./4.")

    temp = tk.IntVar()
    VEParsFormat = temp.get()
    Lambda_TRT_cbox=ttk.Combobox(cntr, width=8, state='readonly', justify="center", values = IO.Lambda_TRT, textvariable=VEParsFormat)
    Lambda_TRT_cbox.current(4)
    Lambda_TRT_cbox.option_add('*TCombobox*Listbox.Justify', 'center') 
    Lambda_TRT_cbox.bind('<<ComboboxSelected>>', lambda _ : Lambda_TRT_Selection(SPs))

    SigSmoothDfcbynsFac_label=tk.Label(ModelParamFrame,text='SigSmoothDfcbynsFac',anchor='n')     
    ToolTip(widget = SigSmoothDfcbynsFac_label, text = "A number between 1e-3 and 1e-1")
    
    GapP2P_label=tk.Label(ModelParamFrame,text='Particle gap',anchor='n')     
    ToolTip(widget = GapP2P_label, text = "Gap between particles to prevent contact. Typically, 1. ")
    
    
    UpdateMotionFac = tk.StringVar(ModelParamFrame, SPs['UpdateMotionFac'])
    UpdateMotionFac_enter = tk.Entry(ModelParamFrame,width=6, textvariable=UpdateMotionFac)        
    UpdateMotionFac_enter.configure({"justify" : "center"})   
    UpdateMotionFac_enter.bind('<Return>',lambda event: AdvParUpdate(SPs))        
    ToolTip(widget = UpdateMotionFac_enter, text = "A number between 0.0005 and 0.1 ")
    
    TargetSlopeFitResults=tk.StringVar(ModelParamFrame, SPs['TargetSlopeFitResults'])
    TargetSlopeFitResults_enter=tk.Entry(ModelParamFrame,width=6, textvariable=TargetSlopeFitResults)        
    TargetSlopeFitResults_enter.configure({"justify" : "center"})   
    TargetSlopeFitResults_enter.bind('<Return>',lambda event: AdvParUpdate(SPs))        
    ToolTip(widget = TargetSlopeFitResults_enter, text = "A number between 0.1 and 10")
    
    PrintIntervalFac=tk.StringVar(ModelParamFrame, SPs['PrintIntervalFac'])
    PrintIntervalFac_enter=tk.Entry(ModelParamFrame,width=6, textvariable=PrintIntervalFac)        
    PrintIntervalFac_enter.configure({"justify" : "center"})   
    PrintIntervalFac_enter.bind('<Return>',lambda event: AdvParUpdate(SPs))        
    ToolTip(widget = PrintIntervalFac_enter, text = "A number between 0.05 and 3")

    MaxtbytRI=tk.StringVar(ModelParamFrame, SPs['MaxtbytRI'])
    MaxtbytRI_enter=tk.Entry(ModelParamFrame,width=6, textvariable=MaxtbytRI)        
    MaxtbytRI_enter.configure({"justify" : "center"})   
    MaxtbytRI_enter.bind('<Return>',lambda event: AdvParUpdate(SPs))        
    ToolTip(widget = MaxtbytRI_enter, text = "A number between 50 and 200")
    
    SigSmoothDfcbynsFac=tk.StringVar(ModelParamFrame, SPs['SigSmoothDfcbynsFac'])
    SigSmoothDfcbynsFacenter=tk.Entry(ModelParamFrame,width=6, textvariable=SigSmoothDfcbynsFac)        
    SigSmoothDfcbynsFacenter.configure({"justify" : "center"})   
    SigSmoothDfcbynsFacenter.bind('<Return>',lambda event: AdvParUpdate(SPs))        
    ToolTip(widget = SigSmoothDfcbynsFacenter, text = "A number between 1e-3 and 1e-1")
    
    
    UpdateGapP2P = tk.StringVar(ModelParamFrame, SPs['Gap_P2P'])
    UpdateGapP2P_enter = tk.Entry(ModelParamFrame,width=6, textvariable=UpdateGapP2P )        
    UpdateGapP2P_enter.configure({"justify" : "center"})   
    UpdateGapP2P_enter.bind('<Return>',lambda event: AdvParUpdate(SPs))        
    ToolTip(widget = UpdateGapP2P_enter, text = "Gap between particles to prevent contact. Typically, 1. ")
    
    if APShow:
        AdvancedSimPars.grid(row=6 ,column=0,columnspan=4,sticky='W')
        UpdateMotionFacLabel.grid(row=7 ,column=0,columnspan=3,sticky='N')
        TargetSlopeFitResultsLabel.grid(row=7 ,column=2,columnspan=4,sticky='N')
        PrintIntervalFacLabel.grid(row=7 ,column=5,columnspan=5,sticky='N')
        MaxtbytRILabel.grid(row=7 ,column=9,columnspan=2,sticky='N')
        Lambda_TRTLabel.grid(row=7 ,column=11,columnspan=3,sticky='N')
        GapP2P_label.grid(row=9 ,column=0,columnspan=7,sticky='N')
        SigSmoothDfcbynsFac_label.grid(row=9 ,column=7,columnspan=7,sticky='N')
        
        UpdateMotionFac_enter.grid(row=8,column=0,columnspan=3,sticky='N') 
        UpdateMotionFac_enter.grid_configure(pady = 5)        
        TargetSlopeFitResults_enter.grid(row=8,column=2,columnspan=4,sticky='N')        
        TargetSlopeFitResults_enter.grid_configure(pady = 5)        
        PrintIntervalFac_enter.grid(row=8,column=5,columnspan=5,sticky='N')        
        PrintIntervalFac_enter.grid_configure(pady = 5)        
        MaxtbytRI_enter.grid(row=8,column=9,columnspan=2,sticky='N')        
        MaxtbytRI_enter.grid_configure(pady = 5)        
        Lambda_TRT_cbox.grid(row=8,column=11,columnspan=3,sticky='N')        
        Lambda_TRT_cbox.grid_configure(pady = 5)        
        UpdateGapP2P_enter.grid(row=10,column=0,columnspan=7,sticky='N')        
        UpdateGapP2P_enter.grid_configure(pady = 5) 
        SigSmoothDfcbynsFacenter.grid(row=10,column=7,columnspan=14,sticky='N')        
        SigSmoothDfcbynsFacenter.grid_configure(pady = 5) 

    elif not APShow:
        AdvancedSimPars.grid_remove()
        UpdateMotionFacLabel.grid_remove()
        TargetSlopeFitResultsLabel.grid_remove()
        PrintIntervalFacLabel.grid_remove()
        MaxtbytRILabel.grid_remove()
        Lambda_TRTLabel.grid_remove()
        GapP2P_label.grid_remove()
        SigSmoothDfcbynsFac_label.grid_remove()

        UpdateMotionFac_enter.grid_remove()
        TargetSlopeFitResults_enter.grid_remove()
        PrintIntervalFac_enter.grid_remove()
        MaxtbytRI_enter.grid_remove()
        Lambda_TRT_cbox.grid_remove()
        UpdateGapP2P_enter.grid_remove()
        SigSmoothDfcbynsFacenter.grid_remove()


    LimitingValuesLabel=tk.Label(ModelParamFrame,text='Limits', fg = 'maroon', anchor='n')    
    betap_Sph_min_Label=tk.Label(ModelParamFrame,text='min', anchor='n')    
    betap_Sph_max_Label=tk.Label(ModelParamFrame,text='max', anchor='n')    
    betap_Sph_Label=tk.Label(ModelParamFrame,text='\u03B2\'', anchor='n')    
    betapp_Sph_min_Label=tk.Label(ModelParamFrame,text='min', anchor='n')
    betapp_Sph_max_Label=tk.Label(ModelParamFrame,text='max', anchor='n')   
    betapp_Sph_Label=tk.Label(ModelParamFrame,text='\u03B2\'\'', anchor='n') 
    
    betap_Sph_min=tk.StringVar(ModelParamFrame, IO.Limits[SPs['VEPars_Choice']][0])
    betap_Sph_min_enter=tk.Entry(ModelParamFrame,width=6, textvariable=betap_Sph_min)        
    betap_Sph_min_enter.configure({"justify" : "center"})   
    betap_Sph_min_enter.bind('<Return>',lambda event: AdvParUpdate(SPs))        

  
    betap_Sph_max=tk.StringVar(ModelParamFrame, IO.Limits[SPs['VEPars_Choice']][1])
    betap_Sph_max_enter=tk.Entry(ModelParamFrame,width=6, textvariable=betap_Sph_max)        
    betap_Sph_max_enter.configure({"justify" : "center"})   
    betap_Sph_max_enter.bind('<Return>',lambda event: AdvParUpdate(SPs))        
    
    betapp_Sph_min=tk.StringVar(ModelParamFrame, IO.Limits[SPs['VEPars_Choice']][2])
    betapp_Sph_min_enter=tk.Entry(ModelParamFrame,width=6, textvariable=betapp_Sph_min)        
    betapp_Sph_min_enter.configure({"justify" : "center"})   
    betapp_Sph_min_enter.bind('<Return>',lambda event: AdvParUpdate(SPs))            
    
    
    betapp_Sph_max=tk.StringVar(ModelParamFrame, IO.Limits[SPs['VEPars_Choice']][3])
    betapp_Sph_max_enter=tk.Entry(ModelParamFrame,width=6, textvariable=betapp_Sph_max)        
    betapp_Sph_max_enter.configure({"justify" : "center"})   
    betapp_Sph_max_enter.bind('<Return>',lambda event: AdvParUpdate(SPs))                

     
    
    if LimShow: 
        LimitingValuesLabel.grid(row=6 ,column=0,columnspan=12,sticky='W')
        betap_Sph_min_Label.grid(row=7 ,column=1,columnspan=2,sticky='W')
        betap_Sph_max_Label.grid(row=7 ,column=3,columnspan=2,sticky='W')
        betapp_Sph_min_Label.grid(row=7 ,column=9,columnspan=2,sticky='W')
        betapp_Sph_max_Label.grid(row=7 ,column=11,columnspan=2,sticky='W')
        
        betap_Sph_Label.grid(row=8 ,column=0,columnspan=1,sticky='W')
        betapp_Sph_Label.grid(row=8 ,column=6,columnspan=1,sticky='W')
        
        betap_Sph_min_enter.grid(row=8 ,column=1,columnspan=2,sticky='W')
        betap_Sph_max_enter.grid(row=8 ,column=3,columnspan=2,sticky='W')
        betapp_Sph_min_enter.grid(row=8 ,column=9,columnspan=2,sticky='W')
        betapp_Sph_max_enter.grid(row=8 ,column=11,columnspan=2,sticky='W')
        
    elif not LimShow:
        LimitingValuesLabel.grid_remove()
        betap_Sph_Label.grid_remove()
        betapp_Sph_Label.grid_remove()
        betap_Sph_min_Label.grid_remove()
        betap_Sph_max_Label.grid_remove()
        betapp_Sph_max_Label.grid_remove()
        betapp_Sph_max_Label.grid_remove()
        betap_Sph_min_enter.grid_remove()
        betap_Sph_max_enter.grid_remove()
        betapp_Sph_min_enter.grid_remove()
        betapp_Sph_max_enter.grid_remove()
        
def MakeBulkParamFrame(cntr, SPs):
    for widget in cntr.winfo_children(): widget.destroy()
    
    f_fund=tk.StringVar(cntr, str("{:.0f}".format(SPs['f0_SI']/1e6)))
    Zq=tk.StringVar(cntr, str("{:.2e}".format(SPs['Zq_SI'])))
    rholiq=tk.StringVar(cntr, '1')
    etaliq=tk.StringVar(cntr, str(SPs['etaabsBulk']))
    tandelBulk=tk.StringVar(cntr, str(SPs['tandelBulk']))
    
    def Enter_BPars(event, choice):
        if choice=='f_fund':
            string = f_Enter.get(); string.strip(); 
            try:
                SPs['f0_SI'] = np.float64(string)*1e6
            except:
                SPs['f0_SI'] = 5e6 
        if choice=='Zq':
            string = Zq_Enter.get(); string.strip(); 
            try:
                SPs['Zq_SI'] = np.float64(string)
            except:
                SPs['f0_SI'] = 8.8e6
        if choice=='etaliq':
            string = etaliq_Enter.get(); string.strip(); 
            try:
                SPs['etaabsBulk'] = np.float64(string)
            except:
                SPs['etaabsBulk'] = 1
        if choice=='tandelBulk':
            string = tandelBulk_Enter.get(); string.strip(); 
            try:
                SPs['tandelBulk'] = np.float64(string)
            except:
                SPs['tandelBulk'] = 1e33
        SPs['delta0_nm']  = 252 * SPs['etaabsBulk']**0.5
        MakeBulkParamFrame(BulkParamFrame, SPs)
    
    tk.Label(cntr,text='f0, MHz',width=wwidth,anchor='n').grid(row=0 ,column=0,columnspan=2,sticky='N')
    tk.Label(cntr,text='Zq, kg/(m\u00B2s)', width=wwidth+2, anchor='n').grid(row=0 ,column=2,columnspan=2,sticky='N')        
    tk.Label(cntr,text='\u03C1 liq, g/cm\u00B3', width=wwidth+2, anchor='n').grid(row=0 ,column=4,columnspan=2,sticky='N')        
    tk.Label(cntr,text='|\u03B7| liq', width=wwidth+2, anchor='n').grid(row=0 ,column=6,columnspan=2,sticky='N')        
    tk.Label(cntr,text='tan(\u03B4) liq', width=wwidth+2, anchor='n').grid(row=0 ,column=8,columnspan=2,sticky='N')      

    f_Enter = tk.Entry(cntr,width=wwidth-5, background='white', foreground='black', textvariable=f_fund)
    f_Enter.grid(row=1,column=0,columnspan=2,sticky='N')
    f_Enter.grid_configure(pady=5)
    f_Enter.bind('<Return>',lambda event: Enter_BPars(event,'f_fund'))
    f_Enter.configure({"justify" : "center"})

    ToolTip(widget = f_Enter, text = "Fundamental frequency")
    
    Zq_Enter = tk.Entry(cntr,width=wwidth, background='black', foreground='white', textvariable=Zq)
    Zq_Enter.grid(row=1,column=2,columnspan=2,sticky='N')
    Zq_Enter.grid_configure(pady=5)
    Zq_Enter.bind('<Return>',lambda event: Enter_BPars(event,'Zq'))
    Zq_Enter.configure({"justify" : "center"})
    ToolTip(widget = Zq_Enter, text = "Acoustic impedance of quartz")
         
    rholiq_Enter = tk.Entry(cntr,width=wwidth - 5 , background='white', foreground='black', state='disabled', textvariable=rholiq)
    rholiq_Enter.grid(row=1,column=4,columnspan=2,sticky='N')
    rholiq_Enter.grid_configure(pady=5)
    rholiq_Enter.bind('<Return>',lambda event: Enter_BPars(event,'rholiq'))
    rholiq_Enter.configure({"justify" : "center"})

    etaliq_Enter = tk.Entry(cntr,width=wwidth-5, background='white', foreground='black', textvariable=etaliq)
    etaliq_Enter.grid(row=1,column=6,columnspan=2,sticky='N')
    etaliq_Enter.grid_configure(pady=5)
    etaliq_Enter.bind('<Return>',lambda event: Enter_BPars(event,'etaliq'))
    etaliq_Enter.configure({"justify" : "center"}) 
    ToolTip(widget = etaliq_Enter, text = "Bulk liquid viscosity")      

    tandelBulk_Enter = tk.Entry(cntr,width=wwidth, background='white', foreground='black', textvariable=tandelBulk)
    tandelBulk_Enter.grid(row=1,column=8,columnspan=2,sticky='N')
    tandelBulk_Enter.grid_configure(pady=5)
    tandelBulk_Enter.bind('<Return>',lambda event: Enter_BPars(event,'tandelBulk'))
    tandelBulk_Enter.configure({"justify" : "center"})        
    ToolTip(widget = tandelBulk_Enter, text = "Loss tangent of the bulk liquid")
    
def MakeSphereParamFrame(cntr,SPs):
    global bgrc, fgrc
        
    def LoopParSelect(SPs,LoopPar):
        global ParCount
        for i in range(len(IO.Parameters[SPs['VEPars_Choice']])):
            if (LoopPar[i].get()==1): 
                if ParCount < 3:
                    if (SPs['Par1str'] == ' ') and (IO.SPsLoopMap[SPs['VEPars_Choice']][i]!=SPs['Par2str']) and (IO.SPsLoopMap[SPs['VEPars_Choice']][i]!=SPs['Par3str']) :
                        SPs['Par1str'] = IO.SPsLoopMap[SPs['VEPars_Choice']][i]
                        ParCount=ParCount+1 
                    elif SPs['Par2str'] == ' ' and (IO.SPsLoopMap[SPs['VEPars_Choice']][i]!=SPs['Par1str']) and (IO.SPsLoopMap[SPs['VEPars_Choice']][i]!=SPs['Par3str']) :
                        SPs['Par2str'] = IO.SPsLoopMap[SPs['VEPars_Choice']][i]
                        ParCount=ParCount+1
                    elif SPs['Par3str'] == ' ' and (IO.SPsLoopMap[SPs['VEPars_Choice']][i]!=SPs['Par1str']) and (IO.SPsLoopMap[SPs['VEPars_Choice']][i]!=SPs['Par2str']) :
                        SPs['Par3str'] = IO.SPsLoopMap[SPs['VEPars_Choice']][i] 
                        ParCount=ParCount+1
            elif (LoopPar[i].get()==0) :                
                if IO.SPsLoopMap[SPs['VEPars_Choice']][i]==SPs['Par1str']:
                    ParCount=ParCount-1
                    SPs['Par1str'] = ' '
                    IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]
                    IO.LoopStepMap[SPs['VEPars_Choice']][i]=1
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i]==SPs['Par2str']: 
                    ParCount=ParCount-1
                    SPs['Par2str'] = ' '
                    IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]
                    IO.LoopStepMap[SPs['VEPars_Choice']][i]=1
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i]==SPs['Par3str']: 
                    ParCount=ParCount-1
                    SPs['Par3str'] = ' '
                    IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]
                    IO.LoopStepMap[SPs['VEPars_Choice']][i]=1    
            
        MakeSphereParamFrame(SphereParamFrame,SPs)
        UpdateSpherePlot(SphereParamFrame, SPs)    
    
    def Enter_SPars(event):
        for i in range(len(IO.Parameters[SPs['VEPars_Choice']])):
            string = EntryArray[i].get(); string.strip();            
            if IO.SPsLoopMap[SPs['VEPars_Choice']][i] =='nSph':
                try:
                    SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]=np.int32(string)
                    if SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]<1: SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]=1
                except:
                    SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]=1
            else:
                try:
                    SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]=np.float64(string)
                except:
                    SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]=-200
            
            if i<11:
                if np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i]) == np.float64(SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]): 
                    IO.LoopStepMap[SPs['VEPars_Choice']][i] = 1
                else:
                    if (LoopPar[i].get()==1) and (np.int32(IO.LoopStepMap[SPs['VEPars_Choice']][i]) == 1): IO.LoopStepMap[SPs['VEPars_Choice']][i] = 2
                    elif(LoopPar[i].get()==0) : IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]
                    
        if SPs['rhoSph']<0.: SPs['rhoSph']=1.0
        if SPs['RSph_nm']<0.: SPs['RSph_nm']=10.0
        if (SPs['ySphbyR']<0.) or (SPs['ySphbyR']>1.0): SPs['ySphbyR']=0.0
        if SPs['nSph'] <1: SPs['nSph'] = 1
        if SPs['etaabscenSphmPas']<0.0: SPs['etaabscenSphmPas']=1.00e4
        if SPs['tandelcenSph']<0.0: SPs['tandelcenSph']=0.1
        if SPs['Jp_FacSph']<0.0 : SPs['Jp_FacSph'] = 0.01 
        if SPs['Jpp_FacSph']<0.0 : SPs['Jpp_FacSph'] = 0.01 
        if SPs['MaxwellRelaxRate_MHz']<0.0 : SPs['MaxwellRelaxRate_MHz'] = 1
        if SPs['tau'] <0.0  : SPs['tau'] = 6e-7
        
        
        if (np.float64(SPs['betap_Sph']) < np.float64(IO.Limits[SPs['VEPars_Choice']][0])) or (np.float64(SPs['betap_Sph']) > np.float64(IO.Limits[SPs['VEPars_Choice']][1])): SPs['betap_Sph'] = 0.0 
        if (np.float64(SPs['betapp_Sph']) < np.float64(IO.Limits[SPs['VEPars_Choice']][2])) or (np.float64(SPs['betapp_Sph']) > np.float64(IO.Limits[SPs['VEPars_Choice']][3])): SPs['betapp_Sph'] = 0.0 

        if (SPs['CovTarget']<=0.0) or (SPs['CovTarget']>1.0): SPs['CovTarget']=0.3 

        InterfaceUpdate(SPs)
        
    def LoopStartSelect(event):
        for i in range(len(IO.Parameters[SPs['VEPars_Choice']])):
            string = LoopStartPar[i].get(); string.strip();
            if IO.SPsLoopMap[SPs['VEPars_Choice']][i] =='nSph' :
                try:
                    IO.LoopStartMap[SPs['VEPars_Choice']][i]=np.int32(string)
                    if IO.LoopStartMap[SPs['VEPars_Choice']][i]<1: IO.LoopStartMap[SPs['VEPars_Choice']][i]=1
                except:
                    IO.LoopStartMap[SPs['VEPars_Choice']][i]=1
            else:
                try:
                    IO.LoopStartMap[SPs['VEPars_Choice']][i]=np.float64(string)
                except:
                    IO.LoopStartMap[SPs['VEPars_Choice']][i]=SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]                    
                    
            if ((IO.SPsLoopMap[SPs['VEPars_Choice']][i] =='rhoSph') and (IO.LoopStartMap[SPs['VEPars_Choice']][i] < 0.)): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['rhoSph']
            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i]  =='rhoSph') and  (IO.LoopStartMap[SPs['VEPars_Choice']][i] < 0)): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['RSph_nm']
            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i] =='ySphbyR') and (IO.LoopStartMap[SPs['VEPars_Choice']][i] <0.0 or IO.LoopStartMap[SPs['VEPars_Choice']][i] >1.0)): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['ySphbyR']
            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i]  == 'nSph') and (IO.LoopStartMap[SPs['VEPars_Choice']][i] <1)): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['nSph'] 
            
            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'etaabscenSphmPas') and (IO.LoopStartMap[SPs['VEPars_Choice']][i] <0.0)): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['etaabscenSphmPas']
            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i]  ==  'tandelcenSph') and (IO.LoopStartMap[SPs['VEPars_Choice']][i] <0.0)): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['tandelcenSph']
            
            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'Jp_FacSph') and (IO.LoopStartMap[SPs['VEPars_Choice']][i] <0.0)): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['Jp_FacSph']
            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'Jpp_FacSph') and (IO.LoopStartMap[SPs['VEPars_Choice']][i] <0.0)) : IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['Jpp_FacSph']

            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'MaxwellRelaxRate_MHz') and (IO.LoopStartMap[SPs['VEPars_Choice']][i] <0.0)): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['MaxwellRelaxRate_MHz']
            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i]  ==  'tau') and (IO.LoopStartMap[SPs['VEPars_Choice']][i] <0.0)): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['tau']
            
            elif IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'betap_Sph':
                if np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i])<np.float64(IO.Limits[SPs['VEPars_Choice']][0]) or np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i])>np.float64(IO.Limits[SPs['VEPars_Choice']][1]): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['betap_Sph']
            elif IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'betapp_Sph':
                if np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i])<np.float64(IO.Limits[SPs['VEPars_Choice']][2]) or np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i])>np.float64(IO.Limits[SPs['VEPars_Choice']][3]): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['betapp_Sph']                
        
            elif ((IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'CovTarget') and ((IO.LoopStartMap[SPs['VEPars_Choice']][i] <=0.0) or (IO.LoopStartMap[SPs['VEPars_Choice']][i] >1.0))): IO.LoopStartMap[SPs['VEPars_Choice']][i] = SPs['CovTarget']
            if np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i]) == np.float64(SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]): 
                IO.LoopStepMap[SPs['VEPars_Choice']][i] = 1
            else:
                if (LoopPar[i].get()==1) and (np.int32(IO.LoopStepMap[SPs['VEPars_Choice']][i]) == 1): IO.LoopStepMap[SPs['VEPars_Choice']][i] = 2

        MakeSphereParamFrame(cntr, SPs)
        UpdateSpherePlot(cntr, SPs)
        
    def LoopStepsSelect(event):
        for i in range(len(IO.Parameters[SPs['VEPars_Choice']])):
            temp=LoopStepPar[i].get()
            try:
                IO.LoopStepMap[SPs['VEPars_Choice']][i]=np.int32(temp)
                if IO.LoopStepMap[SPs['VEPars_Choice']][i]<1: IO.LoopStepMap[SPs['VEPars_Choice']][i]=1
            except:
                IO.LoopStepMap[SPs['VEPars_Choice']][i]=1
            if np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i]) == np.float64(SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]) : 
                IO.LoopStepMap[SPs['VEPars_Choice']][i] = 1
            else:
                if (LoopPar[i].get()==1) and (np.int32(IO.LoopStepMap[SPs['VEPars_Choice']][i]) == 1): IO.LoopStepMap[SPs['VEPars_Choice']][i] = 2
        MakeSphereParamFrame(cntr, SPs)
        UpdateSpherePlot(cntr, SPs)
                
                
    #==============================================================================            
    for widget in cntr.winfo_children(): widget.destroy()
    
    LoopPar=[]              #checkbuton 
    EntryArray=[]           #SPs
    LoopStartArray=[]       #Loop starting values
    LoopStepArray=[]        #LoopSteps
    
    LoopTargetPar=[]
    LoopStartPar=[]
    LoopStepPar=[]
    LStateMap=[]
      
    tk.Label(cntr,text='Loop'             ,anchor='n').grid(row=0 ,column=2,columnspan=2,sticky='N')    
    tk.Label(cntr,text='Start'             ,anchor='n').grid(row=0 ,column=4,columnspan=4,sticky='N')        
    tk.Label(cntr,text='Stop'             ,anchor='n').grid(row=0 ,column=8,columnspan=4,sticky='N')
    tk.Label(cntr,text='Steps'             ,anchor='n').grid(row=0 ,column=12,columnspan=4,sticky='N')
    tk.Label(cntr,text='Spheres Plot'             ,anchor='n').grid(row=0 ,column=20,columnspan=4,sticky='N')
    
    llabel=[]
    for i in range(len(IO.Parameters[SPs['VEPars_Choice']])):
        llabel.append(tk.Label(cntr,text=IO.Parameters[SPs['VEPars_Choice']][i],anchor='n'))
        llabel[i].grid(row=i+1 ,column=0,columnspan=2,sticky='N')
        
        if IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'etaabscenSphmPas':
            LoopTargetPar.append(tk.StringVar(cntr,"{:.2e}".format(SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]])))
            LoopStartPar.append(tk.StringVar(cntr, "{:.2e}".format(np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i]))))
            ToolTip(widget = llabel[i], text = "1e4 corresponds to 1.9 GPa")      
        elif IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'betap_Sph' or IO.SPsLoopMap[SPs['VEPars_Choice']][i] ==  'betapp_Sph':
            LoopTargetPar.append(tk.StringVar(cntr,"{:.2f}".format(SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]])))
            LoopStartPar.append(tk.StringVar(cntr,"{:.2f}".format(np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i]))))        
            ToolTip(widget = llabel[i], text = "Values constraiend by limits")      
        elif IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'tau':
            LoopTargetPar.append(tk.StringVar(cntr,"{:.2e}".format(np.float64(SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]))))
            LoopStartPar.append(tk.StringVar(cntr,"{:.2e}".format(np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i]))))
            ToolTip(widget = llabel[i], text = "~ e-8, the material is mostly liquid; ~ e-6, mostly solid")      
        elif IO.SPsLoopMap[SPs['VEPars_Choice']][i] =='nSph':
            LoopTargetPar.append(tk.StringVar(cntr,"{:.0f}".format(np.int32(SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]))))
            LoopStartPar.append(tk.StringVar(cntr,"{:.0f}".format(np.int32(IO.LoopStartMap[SPs['VEPars_Choice']][i]))))
        else:
            LoopTargetPar.append(tk.StringVar(cntr,"{:.2f}".format(np.float64(SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i]]))))
            LoopStartPar.append(tk.StringVar(cntr,"{:.2f}".format(np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i]))))
        if IO.SPsLoopMap[SPs['VEPars_Choice']][i]==SPs['Par1str'] or IO.SPsLoopMap[SPs['VEPars_Choice']][i]==SPs['Par2str'] or IO.SPsLoopMap[SPs['VEPars_Choice']][i]==SPs['Par3str']:      
            if IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'etaabscenSphmPas':LoopStartPar[i].set("{:.2e}".format(np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i])))  
            elif IO.SPsLoopMap[SPs['VEPars_Choice']][i] == 'tau':LoopStartPar[i].set("{:.2e}".format(np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i])))  
            elif IO.SPsLoopMap[SPs['VEPars_Choice']][i] =='nSph':LoopStartPar[i].set("{:.0f}".format(int(IO.LoopStartMap[SPs['VEPars_Choice']][i])))  
            else: LoopStartPar[i].set("{:.2f}".format(np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][i]))) #Here, set the starting loop prameter to the default. 
            LoopPar.append(tk.IntVar(value=1))            
        else:
            LoopPar.append(tk.IntVar(value=0))
        if (LoopPar[i].get()==1) : LState='normal' 
        elif (LoopPar[i].get()==0): LState='disabled'
        LStateMap.append('normal')
        
        LoopStepPar.append(tk.IntVar(cntr,np.int32(IO.LoopStepMap[SPs['VEPars_Choice']][i])))
        tk.Checkbutton(cntr, variable = LoopPar[i], onvalue = 1, offvalue = 0,  state=LStateMap[i], command=lambda : LoopParSelect(SPs,LoopPar)).grid(row=i+1,column=2,columnspan=2,sticky='S')
        
        EntryArray.append(tk.Entry(cntr, width=wwidth, background=bgrc, foreground=fgrc, state=LStateMap[i], textvariable=LoopTargetPar[i]))
        EntryArray[i].grid(row=i+1,column=4,columnspan=1,sticky='S')                
        EntryArray[i].configure({"justify" : "center"})
        EntryArray[i].bind('<Return>',lambda event: Enter_SPars(event))
        
        LoopStartArray.append(tk.Entry(cntr, width=wwidth, background=bgrc, foreground=fgrc, state=LState, textvariable=LoopStartPar[i]))
        LoopStartArray[i].grid(row=i+1,column=8,columnspan=1,sticky='S')
        LoopStartArray[i].configure({"justify" : "center"})
        LoopStartArray[i].bind('<Return>',lambda event: LoopStartSelect(event))
        
        LoopStepArray.append(tk.Entry(cntr,width=wwidth-2, background=bgrc, foreground=fgrc, state=LState, textvariable=LoopStepPar[i]))
        LoopStepArray[i].grid(row=i+1,column=12,columnspan=1,sticky='S')
        LoopStepArray[i].configure({"justify" : "center"})
        LoopStepArray[i].bind('<Return>',lambda event: LoopStepsSelect(event))

#------------------------------------------------------------------------------------------
    
def MakeSimTrackFrame(cntr, SPs, inPars, flag):
    global wwidth, c, a, b
    ffg = 'black'
    cconst=60
    cplus = 200
    rheader = 20
    rcurrent = 50
    rprevious = 80
    if flag == 'Clear': 
        for widget in cntr.winfo_children(): widget.destroy()
        a =  np.zeros(len(IO.Parameters[SPs['VEPars_Choice']])+3, dtype=int) 
        b =  np.zeros(len(IO.Parameters[SPs['VEPars_Choice']])+3, dtype=int) 
        cntr.update_idletasks()
        c=tk.Canvas(cntr,width=cntr.winfo_width()-20, height=cntr.winfo_height()-30) #, bg='#ffffff'
        c.grid(column=0, row=0)
        c.create_text(100, rheader, text="Parameter:     ", fill='black')
        c.create_text(100, rcurrent, text="Current step:  ", fill='black')
        c.create_text(100, rprevious, text="Previous step: ", fill='black')
        for i in range(0, len(IO.Parameters[SPs['VEPars_Choice']])+3):
            if i == 0 : 
                if SPs['navg'] == 1: ffg = 'black'
                else: ffg = 'red'
                c.create_text(cconst*i+cplus, rheader,  text="    Average    ", fill=ffg)                
            elif i == 1 : 
                if SPs['novt'] == 1: ffg = 'black'
                else: ffg = 'red'
                c.create_text(cconst*i+cplus, rheader,  text="  n  ", fill=ffg)
            elif i == 2 :     
                ffg = 'black'
                c.create_text(cconst*i+cplus, rheader,  text="  Dx, nm  ", fill=ffg)
            elif i > 2:
                if IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par1str']:
                    ffg = 'red'
                    c.create_text(cconst*i+cplus, rheader,  text=IO.Parameters[SPs['VEPars_Choice']][i-3], fill=ffg)
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par2str']:
                    ffg = 'red'
                    c.create_text(cconst*i+cplus, rheader,  text=IO.Parameters[SPs['VEPars_Choice']][i-3], fill=ffg) 
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par3str'] :
                    ffg = 'red'
                    c.create_text(cconst*i+cplus, rheader,  text=IO.Parameters[SPs['VEPars_Choice']][i-3], fill=ffg)
                else:
                    ffg = 'black'
                    c.create_text(cconst*i+cplus, rheader,  text=IO.Parameters[SPs['VEPars_Choice']][i-3], fill=ffg) 
    elif flag == 'Current':         
        for i in range(0, len(IO.Parameters[SPs['VEPars_Choice']])+3):
            if i == 0 : 
                if SPs['navg'] == 1: ffg = 'black'
                else: ffg = 'red'             
                if a[i] != 0: c.delete(a[i])
                a[i] = c.create_text(cconst*i+cplus, rcurrent,  text=inPars[i], fill=ffg)                                                
            elif i == 1 : 
                if SPs['novt'] == 1: ffg = 'black'
                else: ffg = 'red'
                if a[i] != 0: c.delete(a[i])
                a[i] = c.create_text(cconst*i+cplus, rcurrent,  text=inPars[i], fill=ffg)                
            elif i == 2 :     
                ffg = 'black'
                if a[i] != 0: c.delete(a[i]) 
                a[i] = c.create_text(cconst*i+cplus, rcurrent,  text="{:.2f}".format(SPs['Dx_nm']), fill=ffg)                
            elif i > 2:
                if IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par1str']:
                    ffg = 'red'
                    if a[i] != 0: c.delete(a[i]) 
                    a[i] = c.create_text(cconst*i+cplus, rcurrent,  text=inPars[4], fill=ffg)                    
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par2str']:
                    ffg = 'red'
                    if a[i] != 0: c.delete(a[i]) 
                    a[i] = c.create_text(cconst*i+cplus, rcurrent,  text=inPars[3], fill=ffg)                     
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par3str'] :
                    ffg = 'red'
                    if a[i] != 0: c.delete(a[i]) 
                    a[i] = c.create_text(cconst*i+cplus, rcurrent,  text=inPars[2], fill=ffg)                    
                else:
                    ffg = 'black'
                    if a[i] != 0: c.delete(a[i]) 
                    a[i] = c.create_text(cconst*i+cplus, rcurrent,  text=SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]], fill=ffg)
    elif flag == 'Previous':         
        for i in range(0, len(IO.Parameters[SPs['VEPars_Choice']])+3):
            if i == 0 : 
                if SPs['navg'] == 1: ffg = 'black'
                else: ffg = 'red'
                if b[i] != 0: c.delete(b[i]) 
                b[i] = c.create_text(cconst*i+cplus, rprevious,  text=inPars[i], fill=ffg)                                                
            elif i == 1 : 
                if SPs['novt'] == 1: ffg = 'black'
                else: ffg = 'red'
                if b[i] != 0: c.delete(b[i]) 
                b[i] = c.create_text(cconst*i+cplus, rprevious,  text=inPars[i], fill=ffg)                
            elif i == 2 :     
                ffg = 'black'
                if b[i] != 0: c.delete(b[i]) 
                b[i] = c.create_text(cconst*i+cplus, rprevious,  text="{:.2f}".format(SPs['Dx_nm']), fill=ffg)                
            elif i > 2:
                if IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par1str']:
                    ffg = 'red'
                    if b[i] != 0: c.delete(b[i]) 
                    b[i] = c.create_text(cconst*i+cplus, rprevious,  text=inPars[4], fill=ffg)                    
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par2str']:
                    ffg = 'red'
                    if b[i] != 0: c.delete(b[i]) 
                    b[i] = c.create_text(cconst*i+cplus, rprevious,  text=inPars[3], fill=ffg)                     
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par3str'] :
                    ffg = 'red'
                    if b[i] != 0: c.delete(b[i]) 
                    b[i] = c.create_text(cconst*i+cplus, rprevious,  text=inPars[2], fill=ffg)                    
                else:
                    ffg = 'black'
                    if b[i] != 0: c.delete(b[i]) 
                    b[i] = c.create_text(cconst*i+cplus, rprevious,  text=SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]], fill=ffg)                                                                    


# =============================================================================

def MakeSimTrackFramePR(cntr, SPs, inParsC, inParsP, flag):
    global wwidth, w, c, a, b, DisplayPar, llines
    
    def MakeDisplayCheckFrame(cntr):
        global DisplayPar, llines
    
        cntr.configure(width = 50, height=58*llines + 25)
        cntr.update_idletasks()
        DisplayPar=tk.IntVar()
        tk.Label(cntr, text='Plot').place(x = 20, y = 25)
        for i in range(llines):
            tk.Radiobutton(cntr, text='', value=i, variable=DisplayPar).place(x = 25, y = 55+55*i)
        
   
    ffg = 'black'
    cconst=60
    cplus = 200
    rheader = 20
    rcurrent = 50
    rprevious = 72
    rconst = 55

    if flag == 'Clear': 
        for widget in cntr.winfo_children(): widget.destroy()
        a =  np.zeros((len(IO.Parameters[SPs['VEPars_Choice']])+3, llines), dtype=int) 
        b =  np.zeros((len(IO.Parameters[SPs['VEPars_Choice']])+3, llines), dtype=int) 
        cntr.configure(height=58*llines + 25)
        cntr.update_idletasks()
        DisplayCheckFrame = tk.LabelFrame(cntr, padx=5, borderwidth  = 0, highlightthickness=0)
        DisplayCheckFrame.grid(row=0,column=0, sticky='N')
        MakeDisplayCheckFrame(DisplayCheckFrame)
        c=tk.Canvas(cntr,width=cntr.winfo_width()-20, height=cntr.winfo_height()-30) #, bg='#ffffff'
        c.grid(column=1, row=0)
        
        c.create_text(100, rheader, text="Parameter:     ", fill='black')
        for i in range(llines):           
            c.create_text(100, rconst*i + rcurrent, text="Current step:  ", fill='black')
            c.create_text(100, rconst*i + rprevious, text="Previous step: ", fill='black')                       

        for i in range(0, len(IO.Parameters[SPs['VEPars_Choice']])+3):
            if i == 0 : 
                if SPs['navg'] == 1: ffg = 'black'
                else: ffg = 'red'
                c.create_text(cconst*i+cplus, rheader,  text="    Average    ", fill=ffg)                
            elif i == 1 : 
                if SPs['novt'] == 1: ffg = 'black'
                else: ffg = 'red'
                c.create_text(cconst*i+cplus, rheader,  text="  n  ", fill=ffg)
            elif i == 2 :     
                ffg = 'black'
                c.create_text(cconst*i+cplus, rheader,  text="  Dx, nm  ", fill=ffg)
            elif i > 2:
                if IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par1str']:
                    ffg = 'red'
                    c.create_text(cconst*i+cplus, rheader,  text=IO.Parameters[SPs['VEPars_Choice']][i-3], fill=ffg)
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par2str']:
                    ffg = 'red'
                    c.create_text(cconst*i+cplus, rheader,  text=IO.Parameters[SPs['VEPars_Choice']][i-3], fill=ffg) 
                elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par3str'] :
                    ffg = 'red'
                    c.create_text(cconst*i+cplus, rheader,  text=IO.Parameters[SPs['VEPars_Choice']][i-3], fill=ffg)
                else:
                    ffg = 'black'
                    c.create_text(cconst*i+cplus, rheader,  text=IO.Parameters[SPs['VEPars_Choice']][i-3], fill=ffg)          
    elif flag == 'Current':   
        for j in range(len(inParsC)):
            for i in range(0, len(IO.Parameters[SPs['VEPars_Choice']])+3):
                if i == 0 : 
                    if SPs['navg'] == 1: ffg = 'black'
                    else: ffg = 'red'             
                    if a[i,j] != 0: c.delete(a[i,j])
                    if b[i,j] != 0: c.delete(b[i,j])
                    if inParsC[j] != []: a[i,j] = c.create_text(cconst*i+cplus, rconst*j + rcurrent,  text=inParsC[j][i], fill=ffg)                                                                   
                    if inParsP[j] != []: b[i,j] = c.create_text(cconst*i+cplus, rconst*j + rprevious,  text=inParsP[j][i], fill=ffg)
                elif i == 1 : 
                   if SPs['novt'] == 1: ffg = 'black'
                   else: ffg = 'red'
                   if a[i,j] != 0: c.delete(a[i,j])
                   if b[i,j] != 0: c.delete(b[i,j])
                   if inParsC[j] != []: a[i,j]= c.create_text(cconst*i+cplus, rconst*j + rcurrent,  text=inParsC[j][i], fill=ffg)                
                   if inParsP[j] != []: b[i,j] = c.create_text(cconst*i+cplus, rconst*j + rprevious,  text=inParsP[j][i], fill=ffg)
                elif i == 2 :     
                    ffg = 'black'
                    if a[i,j] != 0: c.delete(a[i,j])
                    if b[i,j] != 0: c.delete(b[i,j])
                    if inParsC[j] != []: a[i,j] = c.create_text(cconst*i+cplus, rconst*j + rcurrent,  text="{:.2f}".format(SPs['Dx_nm']), fill=ffg)
                    if inParsP[j] != []: b[i,j] = c.create_text(cconst*i+cplus, rconst*j + rprevious,  text="{:.2f}".format(SPs['Dx_nm']), fill=ffg)
                elif i > 2:
                    if IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par1str']:
                        ffg = 'red'
                        if a[i,j] != 0: c.delete(a[i,j])
                        if b[i,j] != 0: c.delete(b[i,j])
                        if inParsC[j] != []: a[i,j]=c.create_text(cconst*i+cplus, rconst*j + rcurrent,  text=inParsC[j][4], fill=ffg)
                        if inParsP[j] != []: b[i,j] = c.create_text(cconst*i+cplus, rconst*j + rprevious,  text=inParsP[j][4], fill=ffg)
                    elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par2str']:
                        ffg = 'red'
                        if a[i,j] != 0: c.delete(a[i,j])
                        if b[i,j] != 0: c.delete(b[i,j])
                        if inParsC[j] != []: a[i,j]=c.create_text(cconst*i+cplus, rconst*j + rcurrent,  text=inParsC[j][3], fill=ffg)
                        if inParsP[j] != []: b[i,j] = c.create_text(cconst*i+cplus, rconst*j + rprevious,  text=inParsP[j][3], fill=ffg)
                    elif IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]==SPs['Par3str']:
                        ffg = 'red'
                        if a[i,j] != 0: c.delete(a[i,j])
                        if b[i,j] != 0: c.delete(b[i,j])
                        if inParsC[j] != []: a[i,j]=c.create_text(cconst*i+cplus, rconst*j + rcurrent,  text=inParsC[j][2], fill=ffg)
                        if inParsP[j] != []: b[i,j] = c.create_text(cconst*i+cplus, rconst*j + rprevious,  text=inParsP[j][2], fill=ffg)
                    else:
                        ffg = 'black'
                        if a[i,j] != 0: c.delete(a[i,j])
                        if b[i,j] != 0: c.delete(b[i,j])
                        if inParsC[j] != []: a[i,j] = c.create_text(cconst*i+cplus, rconst*j + rcurrent,  text=SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]], fill=ffg)
                        if inParsP[j] != []: b[i,j] = c.create_text(cconst*i+cplus, rconst*j + rprevious,  text=SPs[IO.SPsLoopMap[SPs['VEPars_Choice']][i-3]], fill=ffg)

#======================================================================================================    
def InterfaceUpdate(SPs):
    MakeModelParamFrame(ModelParamFrame,SPs)
    MakeBulkParamFrame(BulkParamFrame, SPs)
    MakeSphereParamFrame(SphereParamFrame,SPs)
    UpdateSpherePlot(SphereParamFrame, SPs)
    MakePlotSim1Frame(PlotSim1Frame, SPs)

def Make_FLBM_Menu(SPs):
    global menubar
    #File (Open, Save, Set Folder)
    #Limits
    #Advanced Parameters
    #Display Results
    #Help (About, Manual)
    
    def SaveSimPars(SPs):
        global GridFactor 
        tempf=asksaveasfilename(filetypes=[("All files", ".*")], initialdir=SPs['folder'], confirmoverwrite=False)
        tempdata={}
        tempdata['SPs'] = SPs
        tempdata['Limits']=IO.Limits[SPs['VEPars_Choice']] 
        tempdata['SPsLoopMap'] = IO.SPsLoopMap[SPs['VEPars_Choice']]
        tempdata['LoopStartMap'] = IO.LoopStartMap[SPs['VEPars_Choice']]
        tempdata['LoopStepMap'] = IO.LoopStepMap[SPs['VEPars_Choice']]
        tempdata['GridFactor'] = GridFactor
        
        np.save(tempf, tempdata, allow_pickle=True)
    
    def OpenSimPars(SPs):
        global GridFactor, menubar
        tempf=askopenfilename(filetypes=[("Simulation Parameter Files", ".npy")], initialdir=SPs['folder'])
        try:
            tempdata = np.load(tempf, allow_pickle='TRUE').item()
            SPs = tempdata['SPs'] 
            IO.Limits[SPs['VEPars_Choice']] = tempdata['Limits']
            IO.SPsLoopMap[SPs['VEPars_Choice']] = tempdata['SPsLoopMap'] 
            IO.LoopStartMap[SPs['VEPars_Choice']] = tempdata['LoopStartMap'] 
            IO.LoopStepMap[SPs['VEPars_Choice']] = tempdata['LoopStepMap']
            GridFactor = tempdata['GridFactor'] 
            InterfaceUpdate(SPs)
        except:
            print('Invalide file type or content')
    
    def SetSimFolder(SPs):
        USerPath = askdirectory(title='Select Folder') 
        SPs['folder']  = USerPath  

    def ResetDefault(SPs):
        SPs = {} 
        SPs = Initialize(SPs, 'reset')
        InterfaceUpdate(SPs)
    
    def MakeDefault(SPs):
        """
        Creates a defalut.npy file with the current simulation parameters to be read by the initialization routine
        """
        global GridFactor 
        tempfolder = os.path.dirname(os.path.abspath(__file__))

        tempf=tempfolder + '\\Defaults'
        tempdata={}
        tempdata['SPs'] = SPs
        tempdata['Limits']=IO.Limits[SPs['VEPars_Choice']]
        tempdata['SPsLoopMap'] = IO.SPsLoopMap[SPs['VEPars_Choice']]
        tempdata['LoopStartMap'] = IO.LoopStartMap[SPs['VEPars_Choice']]
        tempdata['LoopStepMap'] = IO.LoopStepMap[SPs['VEPars_Choice']]
        tempdata['GridFactor'] = GridFactor
        
        np.save(tempf, tempdata, allow_pickle=True)    
        
    def LimitsOpen():
        global APShow, LimShow 
        if LimShow == False:
            LimShow   = True
            APShow = False
        elif LimShow == True:
            LimShow   = False
        
        MakeModelParamFrame(ModelParamFrame, SPs)
    
    def AdvParOpen(SPs):
        global APShow, LimShow 
        
        if APShow == False:
            APShow  = True
            LimShow = False
        elif APShow == True:
            APShow  = False
            
        MakeModelParamFrame(ModelParamFrame, SPs)
                        
    def SimResDisplay(SPs):
        fname  = askopenfilename(filetypes=[("Simulation Results", ".txt")], initialdir=SPs['folder'])
        if fname != '':   
            menubar.entryconfig(4, state="disabled")
            Lib_DResults.DRWin(Root, fname, os.path.dirname(os.path.abspath(__file__)))
        
        
    
    def Help_About(): 
        global FBLMVersion
        def callback(url): webbrowser.open_new(url)
    
        AboutWin = tk.Toplevel()
        AboutWin.title("About fLBSim")
        AboutWin.iconbitmap("flbm.ico")
        
        Root.title('fLBM simulation interface '+ FBLMVersion)
        
        AboutText = tk.Label(AboutWin, text='Frequeny-Domain Lattice Boltzman Simulations \n Core by Diethelm Johannsmann, TU Clausthal \n GUI and parallelization by Ilya Reviakine, AWSensors and University of Washington, Seattle \n Acceleration by Viktor Vanoppen and Paul Hausner, Uppsala \n \n https://www.pc.tu-clausthal.de')                               
        AboutText.pack(ipadx=50, ipady=10, fill='both', expand=True)
        AboutText.bind("<Button-1>", lambda e: callback("https://www.pc.tu-clausthal.de/en/research"))                
        AboutWinbutton = tk.Button(AboutWin, text="OK", command=AboutWin.destroy)
        AboutWinbutton.pack(pady=10, padx=10, ipadx=20, side='bottom')
        
    
    def ManualOpen(SPs):
        fname = os.path.dirname(os.path.abspath(__file__)) + '\\FLBMHelp.pdf'
        webbrowser.open_new(fname)
    
    menubar = tk.Menu(Root)
    Root.config(menu=menubar) 
    file_menu = tk.Menu(menubar)
    file_menu.add_command(label='Save SimPars',command= lambda : SaveSimPars(SPs))
    file_menu.add_command(label='Open SimPars',command=lambda : OpenSimPars(SPs))
    file_menu.add_command(label='Set Folder',command=lambda : SetSimFolder(SPs))
    file_menu.add_command(label='Make Default',command=lambda : MakeDefault(SPs))
    file_menu.add_command(label='Reset SimPars',command=lambda : ResetDefault(SPs))
    menubar.add_cascade(label='File',menu=file_menu)
    
    menubar.add_command(label='Limits'       ,command=LimitsOpen)
    menubar.add_command(label='Advanced SimPars',command=lambda : AdvParOpen(SPs))
    menubar.add_command(label='DisplaySimResults',command=lambda : SimResDisplay(SPs))
    
    help_menu = tk.Menu(menubar)
    help_menu.add_command(label='About',command=Help_About)
    help_menu.add_command(label='Manual',command=lambda : ManualOpen(SPs))
    menubar.add_cascade(label='Help',menu=help_menu)
    
    

#======================================================================================================        
def Run_Sim(SPs):
    global RSButton, SimType
    global fig1, fig2, fig3, fig4, fig5
    global canvas1, canvas2, canvas3, canvas4, canvas5
    global PlotsMap, llines, DisplayPar
    
    #global toUpdate, CurrentParsPars, PreviousParsPars, FileListPar 
    
    PlotsMap=['Plot_LinkProps_3D','Plot_RI','Plot_MotionPars_RI','Plot_Fields_Vertical','Plot_Fields_Horizontal']
    
    def fnMask(s, s1):
        #compares params from file names excluding the overtone par
        answ=False
        if s[0]==s1[0] and s[2] == s1[2] and s[3] == s1[3] and s[4] == s1[4]: answ=True
        return answ

    def fnParse(fn):
        s = os.path.basename(fn)[:-4]
        s=s.split('_') #Need the last 5 elements of this list
        s = s[-5:]
        return s
    
    def tAppend(f, inPars, CurrentParsPars, PreviousParsPars, FileListPar, plotkey):
        global PlotsMap, llines
        
        toUpdate = False
        fIndex  = -1
        if plotkey=='Plot_LinkProps_3D':
            FileListPar['Plot_LinkProps_3D'].append(f)
            toUpdate = False
        elif plotkey == 'Plot_RI': 
            if inPars not in CurrentParsPars: 
                try:
                    fIndex = CurrentParsPars.index([])
                except: 
                    CurrentParsPars.append(inPars)
                    FileListPar['Plot_RI'].append(f)
                    for c in PlotsMap:
                        if c != 'Plot_RI' and c != 'Plot_LinkProps_3D':
                            FileListPar[c].append('')
                    toUpdate = True
                else:
                    fIndex = CurrentParsPars.index([])
                    CurrentParsPars[fIndex] = inPars
                    FileListPar['Plot_RI'][fIndex] = f                                    
                    toUpdate = True                    
            else:   
                fIndex = CurrentParsPars.index(inPars)
                FileListPar['Plot_RI'][fIndex] = f                                    
                toUpdate = True
        elif plotkey == 'Plot_MotionPars_RI':
            fIndex = CurrentParsPars.index(inPars)
            FileListPar['Plot_MotionPars_RI'][fIndex] = f
        elif plotkey == 'Plot_Fields_Horizontal':
            try:
                fIndex = CurrentParsPars.index(inPars)
                PreviousParsPars[fIndex] = inPars
                CurrentParsPars[fIndex] = []            
                FileListPar['Plot_RI'][fIndex]=[]
                FileListPar['Plot_MotionPars_RI'][fIndex]=[]
                FileListPar['Plot_Fields_Horizontal'][fIndex]=f
                toUpdate = True
            except:
                pass                
        elif plotkey == 'Plot_Fields_Vertical':        
            try:
                fIndex = PreviousParsPars.index(inPars)
                FileListPar['Plot_Fields_Vertical'][fIndex]=f    
            except:
                pass
        elif plotkey == 'Error':        
            fIndex = CurrentParsPars.index(inPars)  
            CurrentParsPars[fIndex] = []            
            FileListPar['Plot_RI'][fIndex]=[]
            FileListPar['Plot_MotionPars_RI'][fIndex]=[]
            toUpdate = True            
        return CurrentParsPars, PreviousParsPars, FileListPar, toUpdate
                     
    def pplot(observer, q, CurrentParsPars, PreviousParsPars, FileListPar):
        global SimType, DisplayPar
        global PlotsMap, llines
        
        toUpdate = False    
        if not observer.is_alive():
            print('Observer dead')
            return
        try:
            file_name = q.get_nowait()
        except queue.Empty:
            pass                        
        else:
            if SimType == 'Regular':
                if 'STOP' in file_name:
                    print('Observer stopping')
                    observer.stop()
                    observer.join()
                    print('Deleting ', file_name)
                    while True:
                        try:
                            os.remove(file_name)
                            print(file_name, ' Removed')
                            break
                        except:
                            time.sleep(0.2)
                             
                    for i in range(menubar.index('end')):
                        menubar.entryconfig(i+1, state="normal")
                    MakeModelParamFrame(ModelParamFrame,SPs)
                    MakeBulkParamFrame(BulkParamFrame, SPs)
                    MakeSphereParamFrame(SphereParamFrame,SPs)
                    UpdateSpherePlot(SphereParamFrame, SPs)
                    return
                elif 'Plot_Fields_Horizontal' in file_name:
                    previousPars =[]
                    previousPars = fnParse(file_name)
                    MakeSimTrackFrame(SimTrackFrame, SPs, previousPars, 'Previous')
                    Plots_3D.Plot_Fields_Horizontal(file_name, fig4)
                    canvas4.draw()
                    plt.close("all")
                elif 'Plot_Fields_Vertical' in file_name:                  
                    Plots_3D.Plot_Fields_Vertical(file_name, fig5)
                    canvas5.draw()
                    plt.close("all")
                elif 'Plot_RI' in file_name:
                    currentPars =[]
                    currentPars = fnParse(file_name)
                    MakeSimTrackFrame(SimTrackFrame, SPs, currentPars, 'Current')
                    Plots_3D.Plot_RI(file_name, fig2)
                    canvas2.draw()
                    plt.close("all")
                elif 'Plot_MotionPars_RI' in file_name:
                    Plots_3D.Plot_MotionPars_RI(file_name, fig3)
                    canvas3.draw()
                    plt.close("all")
                elif 'Plot_LinkProps_3D' in file_name:                
                    Plots_3D.Plot_LinkProps_3D(file_name, fig1)               
                    canvas1.draw()
                    plt.close("all")                
            elif SimType == 'Parallel':
                toplt = DisplayPar.get()                                        
                if 'STOP' in file_name:
                    print('Observer stopping')
                    observer.stop()
                    observer.join()
                    print('Deleting ', file_name)
                    while True:
                        try:
                            os.remove(file_name)
                            print(file_name, ' Removed')
                            break
                        except:
                            time.sleep(0.2)
                             
                    for i in range(menubar.index('end')):
                        menubar.entryconfig(i+1, state="normal")
                    MakeModelParamFrame(ModelParamFrame,SPs)
                    MakeBulkParamFrame(BulkParamFrame, SPs)
                    MakeSphereParamFrame(SphereParamFrame,SPs)
                    UpdateSpherePlot(SphereParamFrame, SPs)
                    return                
                elif 'Plot_LinkProps_3D' in file_name:
                    currentPars =[]
                    currentPars = fnParse(file_name)
                    tAppend(file_name, currentPars, CurrentParsPars, PreviousParsPars, FileListPar, 'Plot_LinkProps_3D')
                    if SPs['PlotSimProgress']: 
                        Plots_3D.Plot_LinkProps_3D(file_name, fig1)
                        canvas1.draw()
                        plt.close("all")
                elif 'Plot_RI' in file_name:
                    currentPars =[]
                    currentPars = fnParse(file_name)
                    CurrentParsPars, PreviousParsPars, FileListPar, toUpdate = tAppend(file_name, currentPars, CurrentParsPars, PreviousParsPars, FileListPar, 'Plot_RI')
                    if SPs['PlotSimProgress']: 
                        if (len(FileListPar['Plot_RI']) > 0):
                            try:
                                if (FileListPar['Plot_RI'][toplt] != []):
                                    Plots_3D.Plot_RI(FileListPar['Plot_RI'][toplt], fig2)
                                    canvas2.draw()
                                    plt.close("all")
                                    temp1 = fnParse(FileListPar['Plot_RI'][toplt])
                                    temp1.pop(1)
                                    temp = np.empty(len(FileListPar['Plot_LinkProps_3D']), dtype = float)
                                    temp = list(map(fnParse, FileListPar['Plot_LinkProps_3D']))
                                    for count, i  in enumerate(temp):
                                        if len(i)>0:
                                            del i[1]
                                    findex = temp.index(temp1)
                                    if SPs['PlotSimProgress']: 
                                        Plots_3D.Plot_LinkProps_3D(FileListPar['Plot_LinkProps_3D'][findex], fig1)
                                        canvas1.draw()
                                        plt.close("all")                            
                            except: pass
                elif 'Plot_MotionPars_RI' in file_name:
                    currentPars =[]
                    currentPars = fnParse(file_name)
                    CurrentParsPars, PreviousParsPars, FileListPar, toUpdate = tAppend(file_name, currentPars, CurrentParsPars, PreviousParsPars, FileListPar, 'Plot_MotionPars_RI')
                    if SPs['PlotSimProgress']: 
                        if (len(FileListPar['Plot_MotionPars_RI']) > 0):
                            try:
                                if (FileListPar['Plot_MotionPars_RI'][toplt] != []):
                                    Plots_3D.Plot_MotionPars_RI(FileListPar['Plot_MotionPars_RI'][toplt], fig3)
                                    canvas3.draw()
                                    plt.close("all")
                            except: pass
                elif 'Plot_Fields_Horizontal' in file_name:                        
                    currentPars =[]
                    currentPars = fnParse(file_name)
                    CurrentParsPars, PreviousParsPars, FileListPar, toUpdate = tAppend(file_name, currentPars, CurrentParsPars, PreviousParsPars, FileListPar, 'Plot_Fields_Horizontal')
                    if SPs['PlotSimProgress']: 
                        if (len(FileListPar['Plot_Fields_Horizontal']) > 0):
                            try:
                                if (FileListPar['Plot_Fields_Horizontal'][toplt] != [] and FileListPar['Plot_Fields_Horizontal'][toplt] !=''):
                                    Plots_3D.Plot_Fields_Horizontal(FileListPar['Plot_Fields_Horizontal'][toplt], fig4)
                                    canvas4.draw()
                                    plt.close("all")
                            except: pass
                elif 'Plot_Fields_Vertical' in file_name:                  
                    currentPars =[]
                    currentPars = fnParse(file_name)
                    CurrentParsPars, PreviousParsPars, FileListPar, toUpdate = tAppend(file_name, currentPars, CurrentParsPars, PreviousParsPars, FileListPar, 'Plot_Fields_Vertical')                        
                    if SPs['PlotSimProgress']: 
                        if (len(FileListPar['Plot_Fields_Vertical']) > 0):
                            try:
                                if (FileListPar['Plot_Fields_Vertical'][toplt] != [] and FileListPar['Plot_Fields_Vertical'][toplt] !=''):
                                    Plots_3D.Plot_Fields_Vertical(FileListPar['Plot_Fields_Vertical'][toplt], fig5)
                                    canvas5.draw()
                                    plt.close("all")
                            except: pass
                elif 'SimError' in file_name:
                    currentPars =[]
                    currentPars = fnParse(file_name)
                    CurrentParsPars, PreviousParsPars, FileListPar, toUpdate = tAppend(file_name, currentPars, CurrentParsPars, PreviousParsPars, FileListPar, 'Error')                        
                if toUpdate: 
                    MakeSimTrackFramePR(SimTrackFrame, SPs, CurrentParsPars, PreviousParsPars, 'Current')
                    toUpdate=False

        Root.after(2, pplot, observer, q, CurrentParsPars, PreviousParsPars, FileListPar)                   
        
    class MyHandler(FileSystemEventHandler):        
        def __init__(self, q):
            self._q = q
            super().__init__()
            
            
        def on_created(self,event):
            if event:
              file_name = event.src_path
              #print('file created: ', file_name)              
              self._q.put(file_name)
    
    if SPs['Par1str'] == ' ' or SPs['Par2str'] == ' ' or SPs['Par3str'] == ' ':    
        print('Need three looping parameters. Nothing to simulate.')

    elif SPs['Par1str'] != ' ' and SPs['Par2str'] != ' ' and SPs['Par3str'] != ' ':
        ns=SPs['ns']
        SPs['novt']=len(ns);
        IO.Set_fname(SPs)
        MakePlotSim1Frame(PlotSim1Frame, SPs)
        
        for child in ModelParamFrame.winfo_children():
            child.configure(state = 'disable')
        for child in BulkParamFrame.winfo_children():
            child.configure(state = 'disable') 
        for child in SphereParamFrame .winfo_children():                            
            if not isinstance(child, tk.Canvas): 
                child.configure(state = 'disable')
        for i in range(menubar.index('end')):
            menubar.entryconfig(i+1, state="disabled")
        

        
        sstart = np.float64(SPs[SPs['Par1str']])
        sstop = np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][IO.SPsLoopMap[SPs['VEPars_Choice']].index((SPs['Par1str']))])
        sstep = np.int32(IO.LoopStepMap[SPs['VEPars_Choice']][IO.SPsLoopMap[SPs['VEPars_Choice']].index((SPs['Par1str']))])
        Pars1=np.linspace(sstart,sstop, sstep)

        sstart = np.float64(SPs[SPs['Par2str']])
        sstop = np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][IO.SPsLoopMap[SPs['VEPars_Choice']].index((SPs['Par2str']))])
        sstep = np.int32(IO.LoopStepMap[SPs['VEPars_Choice']][IO.SPsLoopMap[SPs['VEPars_Choice']].index((SPs['Par2str']))])
        Pars2=np.linspace(sstart,sstop, sstep)

        sstart = np.float64(SPs[SPs['Par3str']])
        sstop = np.float64(IO.LoopStartMap[SPs['VEPars_Choice']][IO.SPsLoopMap[SPs['VEPars_Choice']].index((SPs['Par3str']))])
        sstep = np.int32(IO.LoopStepMap[SPs['VEPars_Choice']][IO.SPsLoopMap[SPs['VEPars_Choice']].index((SPs['Par3str']))])
        Pars3=np.linspace(sstart,sstop, sstep)
        
        SPs['nPar1']=len(Pars1);
        SPs['nPar2']=len(Pars2);
        SPs['nPar3']=len(Pars3);        
        ns=SPs['ns']
        SPs[SPs['Par1str'] + 's'] = Pars1
        SPs[SPs['Par2str'] + 's'] = Pars2
        SPs[SPs['Par3str'] + 's'] = Pars3               
        SPs['novt']=len(ns);
              
        IO.Set_fname(SPs)
        fname = '"'+IO.ParallelSim(SPs)+'"'
        print('Simulation parameters file: ',fname)
        SimTrackFrame.grid(row=4,column=0,columnspan=6,sticky='N')
        if SimType == 'Regular': 
            scripttoexecute="fblmn.py"
            MakeSimTrackFrame(SimTrackFrame, SPs, [],'Clear') 
        elif SimType == 'Parallel': 
            scripttoexecute="fblmpar.py" 
            MakeSimTrackFramePR(SimTrackFrame, SPs, [], [],'Clear')    
        
        command = f"python {scripttoexecute} {fname}"
        pfolder = SPs['folder'] + '\\tmpplot'
        if not os.path.exists(pfolder): os.mkdir(pfolder)   
        f=open(pfolder+'\\errors.txt', "w")        
        
        Popen(command, creationflags=CREATE_NEW_CONSOLE, stderr=f)
        
        if __name__ == "__main__":
            FileListPar = {}
            FileListPar['Plot_LinkProps_3D'] = []
            FileListPar['Plot_RI'] = []
            FileListPar['Plot_MotionPars_RI'] = []
            FileListPar['Plot_Fields_Vertical'] = []
            FileListPar['Plot_Fields_Horizontal'] = []
            
            CurrentParsPars = []
            PreviousParsPars = [[] for _ in range(llines)]

            observer=Observer()
            q = queue.Queue()                              
            observer.schedule(MyHandler(q), pfolder, recursive=False)   
            observer.start()
            print('Observer starting')
            Root.after(1, pplot, observer, q, CurrentParsPars, PreviousParsPars, FileListPar)
                           
                            
#======================================================================================================
IO.Read_Config_Interface()
Root = tk.Tk()
string_for_size =  str(IO.WindowWidth)+'x'+\
                   str(IO.WindowHeight)+'+'+\
                   str(IO.WindowLeft)+'+'+\
                   str(IO.WindowTop)
Root.geometry(string_for_size)
Root.protocol("WM_DELETE_WINDOW",onClose)
Root.title('fLBM simulation interface '+ FBLMVersion)
Root.iconbitmap("flbm.ico"); 

MainFr=tk.Frame(Root)
MainFr.grid(row=0, column=0,sticky='nw')
MainCanvas=tk.Canvas(MainFr)
MainCanvas.grid(row=0,column=0,sticky='nw')
scrlbrV=tk.Scrollbar(MainFr, orient="vertical",command=MainCanvas.yview)
scrlbrH=tk.Scrollbar(MainFr, orient="horizontal",command=MainCanvas.xview)
MainCanvas.configure(yscrollcommand=scrlbrV.set)
MainCanvas.configure(xscrollcommand=scrlbrH.set)
scrlbrV.grid(row=0,column=2,sticky='ns')
scrlbrH.grid(row=1,column=0,sticky='we')
innerFrame=tk.Frame(MainCanvas)
innerFrame.grid(row=0, column=0, sticky='nw')
innerFrame.rowconfigure(1, weight=1)

MainCanvas.create_window((0,0),window=innerFrame,anchor='nw')
MainFr.bind('<Configure>',ScrollFunction)
MainCanvas.bind_all('<MouseWheel>', on_mousewheel)
Root.bind('<Configure>',onConfigure) 

SPs = {} 
SPs = Initialize(SPs, 'default') 

ModelParamFrame   = tk.LabelFrame(innerFrame, padx=15, text='Simulation Parameters'         ,foreground= 'maroon'); ModelParamFrame.grid(row=0,column=0,sticky='N')
BulkParamFrame = tk.LabelFrame(innerFrame, padx=15, text='Bulk and Sensor Parameters'         ,foreground= 'maroon'); BulkParamFrame.grid(row=1,column=0,sticky='N')
SphereParamFrame   = tk.LabelFrame(innerFrame, padx=15, text='Sphere Parameters'         ,foreground= 'maroon'); SphereParamFrame.grid(row=2,column=0,sticky='N')
PlotSim1Frame = tk.LabelFrame(innerFrame, padx=5, text='Simulation Plots'         ,foreground= 'maroon'); PlotSim1Frame.grid(row=0,column=1, rowspan=3,sticky='N')
SimTrackFrame = tk.LabelFrame(innerFrame, width = 1100, height = 120, padx=5, text='Simulation Progress'         ,foreground= 'maroon'); 
SimTrackFrame.grid(row=4,column=0,columnspan=6,sticky='N'); SimTrackFrame.grid_propagate(False)
SimTrackFrame.grid_remove()

InterfaceUpdate(SPs)
Make_FLBM_Menu(SPs)

Root.mainloop()