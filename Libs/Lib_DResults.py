import os
import pandas as pd
import tkinter as tk
import matplotlib.pyplot as plt
import configparser
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Libs import Lib_IO as IO

global DRWMainCanvas, PSymbols, PColors, OvtLabels, Overtones, OvertoneCheckVar, xaCheckVar, ParCheckVar, ParCheckVar1, ParCheckVar2
global dataForPlotting, canvas, fig 

PColors = ['r', 'grey', 'b', 'brown', 'navy', 'orange','olive','fuchsia', 'cyan','darksalmon']
PSymbols = ['o', 's', 'v', 'd', 'x', '+', '*', '.', 'p', '^', 'h']
OSymbols = ['x', 's', 'o', '^', ]


def onClose(cntr):
    cntr.nametowidget('.!menu').entryconfig(4, state="normal")
    try: DRWinroot.destroy() 
    except: pass

def on_mousewheel(event):
    global DRWMainCanvas
    DRWMainCanvas.yview_scroll(int(-1*(event.delta/120)), "units")

def ScrollFunction(event):
    global DRWMainCanvas
    w = DRWinroot.winfo_width()
    h = DRWinroot.winfo_height()
    DRWMainCanvas.configure(scrollregion=DRWMainCanvas.bbox("all"),width=w-25,height=h-25)

    # IO.WindowWidth  = w
    # IO.WindowHeight = h
    # IO.WindowTop    = Root.winfo_y()
    # IO.WindowLeft   = Root.winfo_x()
    # IO.Write_Config_Interface()
    
def onConfigure(event):
    global DRWMainCanvas
    def update_size():
        global DRWMainCanvas
        w = DRWinroot.winfo_width()
        h = DRWinroot.winfo_height()
        DRWMainCanvas.configure(scrollregion=DRWMainCanvas.bbox("all"),width=w-25,height=h-25)

    w = DRWinroot.winfo_width()
    h = DRWinroot.winfo_height()
    if w > 1300: 
        w=1300
        DRWinroot.geometry(str(int(w)) +'x' + str(int(h)) + '+15+15')
    if event.widget == DRWinroot:     
            if getattr(DRWinroot, "_after_id", None): 
                DRWinroot.after_cancel(DRWinroot._after_id)
            DRWinroot._after_id=DRWinroot.after(100, update_size)
    # IO.WindowWidth  = w
    # IO.WindowHeight = h+50
    # IO.WindowTop    = Root.winfo_y()
    # IO.WindowLeft   = Root.winfo_x()
    # IO.Write_Config_Interface()


def MakeFParamFrame(cntr, f):
    for widget in cntr.winfo_children(): widget.destroy()    
    w = DRWinroot.winfo_width()
    p, f1 = os.path.split(f)
    tk.Label(cntr,text=f1, justify="center", width = int(w/8.5), anchor='n').grid(row=0 ,column=0, sticky='N')
    
def MakeSParamFrame (cntr, SPs):
    for widget in cntr.winfo_children(): widget.destroy()    

    tk.Label(cntr,text='Problem:',anchor='n').grid(row=0 ,column=0,columnspan=2,sticky='N')
    tk.Label(cntr,text= SPs['problemtype'],fg = 'darkblue', anchor='n').grid(row=1 ,column=0,columnspan=2,sticky='N')
    
    tk.Label(cntr,text='Averages:   ',anchor='n').grid(row=0 ,column=2,columnspan=2,sticky='N')
    tk.Label(cntr,text= SPs['navg'],fg = 'darkblue', anchor='n').grid(row=1 ,column=2,columnspan=2,sticky='N')
    
    tk.Label(cntr,text='\u0394x, nm',anchor='n').grid(row=0 ,column=4,columnspan=2,sticky='N')
    tk.Label(cntr,text="{:.2f}".format(float(SPs['dx_nm'])),fg = 'darkblue', anchor='n').grid(row=1 ,column=4,columnspan=2,sticky='N')
    
    tk.Label(cntr,text='f0, MHz',anchor='n').grid(row=0 ,column=6,columnspan=2,sticky='N')
    tk.Label(cntr,text=str("{:.0f}".format(float(SPs['f0_SI'.lower()])/1e6)),fg = 'darkblue', anchor='n').grid(row=1 ,column=6,columnspan=2,sticky='N')
       
    tk.Label(cntr,text='Zq, kg/(m\u00B2s)',anchor='n').grid(row=0 ,column=8,columnspan=2,sticky='N')
    tk.Label(cntr,text=str("{:.2e}".format(float(SPs['Zq_SI'.lower()]))), fg = 'darkblue', anchor='n').grid(row=1 ,column=8,columnspan=2,sticky='N')

    tk.Label(cntr,text='\u03C1 liq, g/cm\u00B3',anchor='n').grid(row=0 ,column=10,columnspan=2,sticky='N')
    tk.Label(cntr,text='1', fg = 'blue', anchor = 'n').grid(row=1 ,column=10,columnspan=2,sticky='N')

    tk.Label(cntr,text='|\u03B7| liq',anchor='n').grid(row=0 ,column=12,columnspan=2,sticky='N')
    tk.Label(cntr,text=str(SPs['etaabsBulk'.lower()]), fg = 'darkblue', anchor='n').grid(row=1 ,column=12,columnspan=2,sticky='N')

    tk.Label(cntr,text='tan(\u03B4) liq',anchor='n').grid(row=0 ,column=14,columnspan=2,sticky='N')
    tk.Label(cntr,text=str(SPs['tandelBulk'.lower()]), fg = 'darkblue', anchor='n').grid(row=1 ,column=14,columnspan=2,sticky='N')
       
    tt = 0
    for i in range(len(IO.Parameters[SPs['vepars_choice']])):
        if IO.SPsLoopMap[SPs['vepars_choice']][i] != SPs['par1str'] and IO.SPsLoopMap[SPs['vepars_choice']][i] != SPs['par2str'] and IO.SPsLoopMap[SPs['vepars_choice']][i] != SPs['par3str']:
            tk.Label(cntr, text=IO.Parameters[SPs['vepars_choice']][i], width = 10, anchor='n').grid(row=2 ,column=2*i - tt,columnspan=2,sticky='N')
            l = SPs[IO.SPsLoopMap[SPs['vepars_choice']][i].lower()]
            tk.Label(cntr, text=l, width = 10, fg = 'darkblue', anchor='n').grid(row=3 ,column=2*i - tt,columnspan=2,sticky='N')
        else:
            tt = tt + 2

def MakePParamFrame(cntr, SPs):
    global Colors, OvtLabels, OvertoneCheckVar, xaCheckVar
    for widget in cntr.winfo_children(): widget.destroy()    

    def OvertoneSelect():
        RD(SPs)
            
    def xaSelect(SPs):
        MakePPFrame1(DRWinroot.nametowidget('.!frame').nametowidget('.!frame.!canvas').nametowidget('.!frame.!canvas.!frame').nametowidget('.!frame.!canvas.!frame.!labelframe4'), SPs)
        RD(SPs)
        

    #----------------------------------------------------------------------------------------
    w = 10            
    novt  = int(SPs['novt'])            
    tk.Label(cntr,text='Overtones: ', width = w, anchor='n').grid(row=0 ,column=0,columnspan=2,sticky='N')
    tk.Label(cntr,text='X-axis: ', width = w, anchor='n').grid(row=1 ,column=0,columnspan=2,sticky='N')
                    
    OvertoneCheckButton=[]
    OvertoneCheckVar=[] 
    for i in range(novt):
        temp = tk.IntVar(cntr, 1)
        OvertoneCheckVar.append(temp)
        OvertoneCheckButton.append(tk.Checkbutton(cntr, width = w-2, text=SPs['ns'][i],variable = OvertoneCheckVar[i],onvalue = 1, offvalue = 0, command=OvertoneSelect))
        OvertoneCheckButton[i].grid(row=0,column=i+2,columnspan=1,sticky='N')

    temp = tk.IntVar(cntr, 2)
    xaCheckVar=temp 
    xaheckButton=[]
    for i in range(3):
        s = 'par'+str(i+1)+'str'
        label = (IO.Parameters[SPs['vepars_choice']][IO.SPsLoopMap[SPs['vepars_choice']].index(SPs[s])])
        xaheckButton.append(tk.Radiobutton(cntr,width = w-2, text=label, variable = xaCheckVar, value=i, command = lambda : xaSelect(SPs)))
        xaheckButton[i].grid(row=1 ,column=i + 2,columnspan=1,sticky='N')


def MakePPFrame1(cntr, SPs):
    global ParCheckVar, ParCheckVar1, ParCheckVar2
    for widget in cntr.winfo_children(): widget.destroy()    
    w = 10

    def ParSelect(SPs):
        RD(SPs)

    def Par1Select(SPs):
        RD(SPs)

    def Par2Select(SPs):
        RD(SPs)

    
    s = xaCheckVar.get()+1
    if s == 3:
        s1 = 1
        s2 = 2
    elif s == 2:
        s1 = 1
        s2 = 3
    elif s == 1:
        s1 = 2
        s2 = 3             
    
    s1 = 'par'+str(s1)+'str'
    s2 = 'par'+str(s2)+'str'
    s = 'par'+str(s)+'str'
    
    label = (IO.Parameters[SPs['vepars_choice']][IO.SPsLoopMap[SPs['vepars_choice']].index(SPs[s1])]) + ': '
    tk.Label(cntr,text=label,anchor='n').grid(row=0, column=0, columnspan=2,sticky='W')
    
    ParCheckButton=[]
    ParCheckVar=[]
    for i in range(len(SPs[SPs[s1].lower()+'s'])):
        temp = tk.IntVar(cntr, 1)
        ParCheckVar.append(temp)
        label = SPs[SPs[s1].lower()+'s'][i]
        ParCheckButton.append(tk.Checkbutton(cntr, width = w-2, text=label, fg = 'darkblue', anchor='n', variable = ParCheckVar[i],onvalue = 1, offvalue = 0, command=lambda : ParSelect(SPs)))
        ParCheckButton[i].grid(row=0, column=2 + i, columnspan=1,sticky='W')

    label = (IO.Parameters[SPs['vepars_choice']][IO.SPsLoopMap[SPs['vepars_choice']].index(SPs[s2])]) + ': '
    tk.Label(cntr,text=label,anchor='n').grid(row=0, column=6 + (len(SPs[SPs[s1].lower()+'s'])), columnspan=2,sticky='W')

    ParCheckButton1=[]
    ParCheckVar1=[]
    for i in range(len(SPs[SPs[s2].lower()+'s'])):
        temp = tk.IntVar(cntr, 1)
        ParCheckVar1.append(temp)
        label = SPs[SPs[s2].lower()+'s'][i]
        ParCheckButton1.append(tk.Checkbutton(cntr, width = w-2, text=label, fg = 'darkblue', anchor='n', variable = ParCheckVar1[i],onvalue = 1, offvalue = 0, command= lambda : Par1Select(SPs)))
        tk.Label(cntr,text=label,anchor='n')
        ParCheckButton1[i].grid(row=0, column=8 + (len(SPs[SPs[s1].lower()+'s'])) + i, columnspan=1,sticky='W')

    label = (IO.Parameters[SPs['vepars_choice']][IO.SPsLoopMap[SPs['vepars_choice']].index(SPs[s])]) + ': '
    tk.Label(cntr,text=label,anchor='n').grid(row=0, column=10 + (len(SPs[SPs[s1].lower()+'s'])) + (len(SPs[SPs[s2].lower()+'s'])), columnspan=2,sticky='W')

    ParCheckButton2=[]
    ParCheckVar2=[]
    for i in range(len(SPs[SPs[s].lower()+'s'])):
        temp = tk.IntVar(cntr, 1)
        ParCheckVar2.append(temp)
        label = SPs[SPs[s].lower()+'s'][i]
        ParCheckButton2.append(tk.Checkbutton(cntr, width = w-2, text=label, fg = 'darkblue',  anchor='n', variable = ParCheckVar2[i],onvalue = 1, offvalue = 0, command= lambda : Par2Select(SPs)))
        tk.Label(cntr,text=label,anchor='n')
        ParCheckButton2[i].grid(row=0, column=12 + (len(SPs[SPs[s1].lower()+'s']))+(len(SPs[SPs[s2].lower()+'s'])) + i, columnspan=1,sticky='W')
    

    
def dLoad(fname, SPs):

    data = pd.read_csv(fname, sep='\t')
    data = data.sort_values(['n', SPs['par1str'], SPs['par2str'], SPs['par3str'], 'iavg'])
    data['CovTarget'] = data['CoverageTrue'].round(3)
    data['CovTarget'] = data['CovTarget'].astype(float)
 
    data.drop(labels = 'iavg', axis = 1, inplace=True) 
    
    SPs['navg'] = int(SPs['navg'])  
    SPs['ns'] = list(map(float, SPs['ns'][1:-1].split(' ')))    
    
    if SPs['par1str'].lower()+'s' == 'covtargets': 
        SPs[SPs['par1str'].lower()+'s'] = data['CovTarget'].unique().round(3)
    else:
        SPs[SPs['par1str'].lower()+'s'] = list(map(float, SPs[SPs['par1str'].lower()+'s'][1:-1].strip().split()))

    if SPs['par2str'].lower()+'s' == 'covtargets': 
        SPs[SPs['par2str'].lower()+'s'] = data['CovTarget'].unique().round(3)
    else:
        SPs[SPs['par2str'].lower()+'s']= list(map(float, SPs[SPs['par2str'].lower()+'s'][1:-1].strip().split()))
        
    if SPs['par3str'].lower()+'s' == 'covtargets': 
        SPs[SPs['par3str'].lower()+'s'] = data['CovTarget'].unique().round(3)
    else:
        SPs[SPs['par3str'].lower()+'s'] = list(map(float, SPs[SPs['par3str'].lower()+'s'][1:-1].strip().split()))
    
    data1=[]
    
    for ovt in SPs['ns']:
        temp = data.loc[(data['n'] == ovt)] 
        for SPs['ipar1'],SPs[SPs['par1str'].lower()] in enumerate(SPs[SPs['par1str'].lower()+'s']):
            temp1 = temp.loc[temp[SPs['par1str']] == SPs[SPs['par1str'].lower()+'s'][SPs['ipar1']]]
            for SPs['ipar2'],SPs[SPs['par2str'].lower()] in enumerate(SPs[SPs['par2str'].lower()+'s']):
                temp2 = temp1.loc[temp1[SPs['par2str']] == SPs[SPs['par2str'].lower()+'s'][SPs['ipar2']]]
                temp2.insert(9,'DfbynErr', 0)
                temp2.insert(11,'DGbynErr', 0)                    
                temp2.insert(13,'Dfratio', 0)
                temp2.insert(14,'DfrErr', 0)                
                temp2.insert(15,'ffratio', 0)
                temp2.insert(16,'ffrErr', 0)
                for SPs['ipar3'],SPs[SPs['par3str'].lower()] in enumerate(SPs[SPs['par3str'].lower()+'s']):
                    temp3 = temp2[temp2[SPs['par3str']] == SPs[SPs['par3str'].lower()+'s'][SPs['ipar3']]]
                    if not temp3.empty: 
                        temp3 = temp2.loc[temp2[SPs['par3str']] == SPs[SPs['par3str'].lower()+'s'][SPs['ipar3']]].mean().transpose()                
                        tErr = temp2.loc[temp2[SPs['par3str']] == SPs[SPs['par3str'].lower()+'s'][SPs['ipar3']]].std()            
                        temp3['DfbynErr'] = tErr['Dfbyn'] 
                        temp3['DGbynErr'] = tErr['DGbyn'] 
                        temp3['Dfratio'] = -1 * temp3['DGbyn'] / temp3['Dfbyn']
                        temp3['DfrErr'] = abs(temp3['Dfratio'])*(abs(temp3['DGbynErr']/temp3['DGbyn']) + abs(temp3['DfbynErr']/temp3['Dfbyn']))
                        data1.append(temp3)
    
    data1=pd.DataFrame(data1)                
    c_drop = ['ietaabscenSphmPas', 'itandelcenSph', 'iCovTarget', 'in', 'nNodes', 'steps', 'tbytRI', 'CompTimeMins']
    data1.drop(labels = c_drop, axis = 1, inplace=True)
    data1.sort_values(['n', SPs['par1str'], SPs['par2str'], SPs['par3str']], inplace=True)

    for SPs['ipar1'],SPs[SPs['par1str'].lower()] in enumerate(SPs[SPs['par1str'].lower()+'s']):
        for SPs['ipar2'],SPs[SPs['par2str'].lower()] in enumerate(SPs[SPs['par2str'].lower()+'s']):        
            for SPs['ipar3'],SPs[SPs['par3str'].lower()] in enumerate(SPs[SPs['par3str'].lower()+'s']):
                for ovt in SPs['ns']:
                    t = data1['Dfbyn'][(data1['n']==ovt) & (data1[SPs['par1str']] == SPs[SPs['par1str'].lower()+'s'][SPs['ipar1']]) & (data1[SPs['par2str']] == SPs[SPs['par2str'].lower()+'s'][SPs['ipar2']]) & (data1[SPs['par3str']] == SPs[SPs['par3str'].lower()+'s'][SPs['ipar3']])].to_numpy()
                    t1 = data1['Dfbyn'][(data1['n']==SPs['ns'][0]) & (data1[SPs['par1str']] == SPs[SPs['par1str'].lower()+'s'][SPs['ipar1']]) & (data1[SPs['par2str']] == SPs[SPs['par2str'].lower()+'s'][SPs['ipar2']]) & (data1[SPs['par3str']] == SPs[SPs['par3str'].lower()+'s'][SPs['ipar3']]) ].to_numpy()
                    et = data1['DfbynErr'][(data1['n']==ovt) & (data1[SPs['par1str']] == SPs[SPs['par1str'].lower()+'s'][SPs['ipar1']]) & (data1[SPs['par2str']] == SPs[SPs['par2str'].lower()+'s'][SPs['ipar2']]) & (data1[SPs['par3str']] == SPs[SPs['par3str'].lower()+'s'][SPs['ipar3']])].to_numpy()
                    et1 = data1['DfbynErr'][(data1['n']==SPs['ns'][0]) & (data1[SPs['par1str']] == SPs[SPs['par1str'].lower()+'s'][SPs['ipar1']]) & (data1[SPs['par2str']] == SPs[SPs['par2str'].lower()+'s'][SPs['ipar2']]) & (data1[SPs['par3str']] == SPs[SPs['par3str'].lower()+'s'][SPs['ipar3']]) ].to_numpy()
                    rratio = 1 - t/t1
                    erratio = abs(rratio)*(abs(et/t)+abs(et1/t1))                    
                    if len(rratio) != 0 :
                        data1['ffratio'][(data1['n']==ovt) & (data1[SPs['par1str']] == SPs[SPs['par1str'].lower()+'s'][SPs['ipar1']]) & (data1[SPs['par2str']] == SPs[SPs['par2str'].lower()+'s'][SPs['ipar2']]) & (data1[SPs['par3str']] == SPs[SPs['par3str'].lower()+'s'][SPs['ipar3']])] = rratio[0]
                        data1['ffrErr'][(data1['n']==ovt) & (data1[SPs['par1str']] == SPs[SPs['par1str'].lower()+'s'][SPs['ipar1']]) & (data1[SPs['par2str']] == SPs[SPs['par2str'].lower()+'s'][SPs['ipar2']]) & (data1[SPs['par3str']] == SPs[SPs['par3str'].lower()+'s'][SPs['ipar3']])] = erratio[0]
    data1.to_clipboard()    
    return data1, SPs

def RD(SPs):
    global PSymbols, PColors, OvtLabels, Overtones, OvertoneCheckVar, xaCheckVar, ParCheckVar, ParCheckVar1, ParCheckVar2
    global dataForPlotting, canvas, fig 
    
    fig.clear()
    
    data = dataForPlotting 
#    print(SPs['par1str'],SPs['par2str'],SPs['par3str'])
    novt  = int(SPs['novt']) 

    if novt>1:        
        ax4 = fig.add_subplot(3,3,1)
        ax5 = fig.add_subplot(3,3,4)
        ax6 = fig.add_subplot(3,3,7) 
    
        ax1 = fig.add_subplot(3,3,2)
        ax2 = fig.add_subplot(3,3,5)
        ax3 = fig.add_subplot(3,3,8) 
        
        ax7 = fig.add_subplot(3,3,3)
        ax8 = fig.add_subplot(3,3,6)
        ax9 = fig.add_subplot(3,3,9)

    else:
        ax1 = fig.add_subplot(2,3,1)
        ax2 = fig.add_subplot(2,3,2)
        ax3 = fig.add_subplot(2,3,3) 
        ax4 = fig.add_subplot(2,3,4) 
                
    
    
    s = xaCheckVar.get()+1
    if s == 3:
        s1 = 1
        s2 = 2
    elif s == 2:
        s1 = 1
        s2 = 3
    elif s == 1:
        s1 = 2
        s2 = 3             

    is1 = 'ipar' + str(s1)
    s1 = 'par'+str(s1)+'str'
    
    l1 = SPs[s1]
    m1 = SPs[s1]

    is2 = 'ipar' + str(s2)
    s2 = 'par'+str(s2)+'str'    
    l2 = SPs[s2]
    m2 = SPs[s2]
    
    is3 = 'ipar' + str(s)        
    s = 'par' + str(s).strip()+'str'
    l = SPs[s]
    m = SPs[s]

    if l == 'CovTarget': 
        l = 'Coverage'
        m = 'CoverageTrue'
    else: 
        l = IO.Parameters[SPs['vepars_choice']][IO.SPsLoopMap[SPs['vepars_choice']].index(l)]

    if l1 == 'CovTarget': 
        l1 = 'Coverage'
        m1 = 'CoverageTrue'
    else: 
        l1 = IO.Parameters[SPs['vepars_choice']][IO.SPsLoopMap[SPs['vepars_choice']].index(l1)]
        
    if l2 == 'CovTarget': 
        l2 = 'Coverage'
        m2 = 'CoverageTrue'
    else: 
        l2 = IO.Parameters[SPs['vepars_choice']][IO.SPsLoopMap[SPs['vepars_choice']].index(l2)]
    
    psindex = 0
    for SPs[is1],SPs[SPs[s1].lower()] in enumerate(SPs[SPs[s1].lower()+'s']):
        if psindex > len(PColors)-1: psindex = 0
        temp = data.loc[data[SPs[s1]] == SPs[SPs[s1].lower()+'s'][SPs[is1]]]    
        l11 = ' '  + str(SPs[SPs[s1].lower()+'s'][SPs[is1]])
        for SPs[is2],SPs[SPs[s2].lower()] in enumerate(SPs[SPs[s2].lower()+'s']):
            if psindex > len(PColors)-1: psindex = 0
            temp1 = temp.loc[temp[SPs[s2]] == SPs[SPs[s2].lower()+'s'][SPs[is2]]]    
            l12 = ' '  + str(SPs[SPs[s2].lower()+'s'][SPs[is2]])
            for iovt in range(novt):                
                c = SPs['ns'][iovt]                
                if OvertoneCheckVar[iovt].get() == 1 and ParCheckVar[SPs[is1]].get() == 1 and ParCheckVar1[SPs[is2]].get() == 1:
                    ax1.errorbar(temp1[temp1.n == c][m], temp1[temp1.n == c]['Dfbyn'], temp1[temp1.n == c]['DfbynErr'], fmt=OSymbols[iovt], mfc='w', ms = 6, label = str("{:.0f}".format(c*float(SPs['f0_SI'.lower()])/1e6)) + ' MHz' + l11 + l12,color = PColors[psindex])
                    ax2.errorbar(temp1[temp1.n == c][m], temp1[temp1.n == c]['DGbyn'], temp1[temp1.n == c]['DGbynErr'], fmt=OSymbols[iovt], mfc='w', ms = 6, label = str("{:.0f}".format(c*float(SPs['f0_SI'.lower()])/1e6)) + ' MHz' + l11 + l12,color = PColors[psindex])
                    ax3.errorbar((-1*temp1[temp1.n == c]['Dfbyn']), temp1[temp1.n == c]['Dfratio'], temp1[temp1.n == c]['DfrErr'], fmt=OSymbols[iovt], mfc='w', ms = 6, label = str("{:.0f}".format(c*float(SPs['f0_SI'.lower()])/1e6)) + ' MHz' + l11 + l12,color = PColors[psindex])
                    if temp1[temp1.n == c]['ffratio'].all() != 0:
                        ax7.errorbar(-1*temp1['Dfbyn'], temp1['ffratio'], temp1['ffrErr'], fmt=OSymbols[iovt], mfc='w', ms = 6, label = l11 + l12,color = 'blue')                        
            psindex += 1                       
    if novt>1:        
        ax4.set_ylabel('$\Delta f/n$ [Hz]'); 
        ax4.set_xticks
        ax4.set_xlabel(('Frequency [MHz]'))
        
        ax5.set_ylabel('$\Delta \Gamma /n$ [Hz]'); 
        ax5.set_xticks
        ax5.set_xlabel(('Frequency [MHz]'))
        
        ax6.set_ylabel('$\Delta \Gamma /(-\Delta f)$'); 
        ax6.set_xticks
        ax6.set_xlabel(('Frequency [MHz]'))

        ax7.set_ylabel('1-$\Delta fn /(\Delta fmin)$'); 
        ax7.set_xticks
        ax7.set_xlabel(('- $\Delta f/n$ [Hz]'))

        
        cindex = 0
        sindex = len(PSymbols)-1
        
        for SPs[is1],SPs[SPs[s1].lower()] in enumerate(SPs[SPs[s1].lower()+'s']):
            if cindex > len(PColors)-1: cindex = 0
            if sindex ==0: sindex = len(PSymbols)-1
            temp = data.loc[data[SPs[s1]] == SPs[SPs[s1].lower()+'s'][SPs[is1]]]    
            l11 = ' '  + str(SPs[SPs[s1].lower()+'s'][SPs[is1]])
            for SPs[is2],SPs[SPs[s2].lower()] in enumerate(SPs[SPs[s2].lower()+'s']):
                if cindex > len(PColors)-1: cindex = 0
                if sindex ==0: sindex = len(PSymbols)-1                
                temp1 = temp.loc[temp[SPs[s2]] == SPs[SPs[s2].lower()+'s'][SPs[is2]]]    
                l12 = ' '  + str(SPs[SPs[s2].lower()+'s'][SPs[is2]])
                for SPs[is3],SPs[SPs[s].lower()] in enumerate(SPs[SPs[s].lower()+'s']):
                    if cindex > len(PColors)-1: cindex = 0
                    if sindex ==0: sindex = len(PSymbols)-1
                    temp2 = temp1.loc[temp1[SPs[s]] == SPs[SPs[s].lower()+'s'][SPs[is3]]]    
                    l13 = ' '  + str(SPs[SPs[s].lower()+'s'][SPs[is3]])
                    if ParCheckVar[SPs[is1]].get()==1 and ParCheckVar1[SPs[is2]].get()==1 and ParCheckVar2[SPs[is3]].get()==1:
                        if psindex > len(PColors)-1: psindex = 0
                        ax4.errorbar(temp2['n']*float(SPs['f0_SI'.lower()])/1e6, temp2['Dfbyn'], temp2['DfbynErr'], fmt=PSymbols[sindex], mfc='w', ms = 6, label = l11 + l12 + l13,color = PColors[cindex])
                        ax5.errorbar(temp2['n']*float(SPs['f0_SI'.lower()])/1e6, temp2['DGbyn'], temp2['DGbynErr'], fmt=PSymbols[sindex], mfc='w', ms = 6, label = l11 + l12 + l13,color = PColors[cindex])
                        negfreq = -1*temp2['Dfbyn']
                        dfratio = -1*temp2['DGbyn']/temp2['Dfbyn']
                        dferror = abs(dfratio)*(abs(temp2['DfbynErr']/negfreq) + abs(temp2['DGbynErr'] / temp2['DGbyn']))
                        ax6.errorbar(temp2['n']*float(SPs['f0_SI'.lower()])/1e6, dfratio, dferror, fmt=PSymbols[sindex], mfc='w', ms = 6, label = l11 + l12 + l13,color = PColors[cindex])                                                                                        
                        cindex+=1
                        sindex-=1                    

        label_params = ax1.get_legend_handles_labels() 
        label_params1 = ax4.get_legend_handles_labels() 
        ax8.legend(*label_params1, loc="center", ncols = 3, prop={"size":7})
        ax8.axis('off')
        ax9.legend(*label_params, loc="center", ncols = 3,  prop={"size":7})
        ax9.axis('off')



    else:
        label_params = ax1.get_legend_handles_labels() 
        ax4.legend(*label_params, loc="center", ncols = 3, prop={"size":8})
        ax4.axis('off')
        
    ax1.set_ylabel('$\Delta f/n$ [Hz]'); 
    ax1.set_xticks
    ax1.set_xlabel(l)
    
    ax2.set_ylabel('$\Delta \Gamma /n$ [Hz]'); 
    ax2.set_xticks
    ax2.set_xlabel(l)
    
    ax3.set_ylabel('$\Delta \Gamma /(-\Delta f)$'); 
    ax3.set_xticks
    ax3.set_xlabel(('- $\Delta f/n$ [Hz]'))        
    
    fig.tight_layout(); 
    plt.close('all')
    canvas.draw()
    
def DRWin(cntr, fname, path):
    global DRWMainCanvas, DRWinroot, dataForPlotting, canvas, fig 
    SPs = {}
    
    DRWinroot = tk.Tk();
    DRWinroot.geometry('1300x800')
    DRWinroot.title('Display Simulation Results')
    DRWinroot.iconbitmap(path + r'\flbm.ico')
    DRWinroot.protocol("WM_DELETE_WINDOW", lambda : onClose(cntr))
    
    DRWMainFr=tk.Frame(DRWinroot)
    DRWMainFr.grid(row=0, column=0,sticky='nw')
    DRWMainCanvas=tk.Canvas(DRWMainFr)
    DRWMainCanvas.grid(row=0,column=0,sticky='nw')
    scrlbrV=tk.Scrollbar(DRWMainFr, orient="vertical",command=DRWMainCanvas.yview)
    scrlbrH=tk.Scrollbar(DRWMainFr, orient="horizontal",command=DRWMainCanvas.xview)
    DRWMainCanvas.configure(yscrollcommand=scrlbrV.set)
    DRWMainCanvas.configure(xscrollcommand=scrlbrH.set)
    scrlbrV.grid(row=0,column=2,sticky='ns')
    scrlbrH.grid(row=1,column=0,sticky='we')
    DRWinnerFrame=tk.Frame(DRWMainCanvas)
    DRWinnerFrame.grid(row=0, column=0, sticky='nw')
    DRWinnerFrame.rowconfigure(1, weight=1)

    DRWMainCanvas.create_window((0,0),window=DRWinnerFrame,anchor='nw')
    DRWMainFr.bind('<Configure>', ScrollFunction)
    DRWMainCanvas.bind_all('<MouseWheel>', on_mousewheel)
    DRWinroot.bind('<Configure>', onConfigure) 
    
    config = configparser.ConfigParser()
    config.read(fname[:-4]+'.cfg')
    config.optionxform = str    
    for section in config.sections():
        for option in config.options(section): SPs[option]=config.get(section,option)                

    dataForPlotting, SPs = dLoad(fname, SPs)
    
    SParamFrame = tk.LabelFrame(DRWinnerFrame, padx=15, text='Simulation Parameters',foreground= 'maroon')
    SParamFrame .grid(row=1,column=0,sticky ='N')
    PParamFrame = tk.LabelFrame(DRWinnerFrame, padx=15, text='Plot Parameters',foreground= 'maroon')
    PParamFrame.grid(row=1,column=1, sticky ='N')
    
    MakePParamFrame(PParamFrame, SPs)
    MakeSParamFrame(SParamFrame, SPs)
    
    FParamFrame = tk.LabelFrame(DRWinnerFrame, padx=15, text='Simulation File Name',foreground= 'maroon')
    FParamFrame.grid(row=0,column=0, columnspan = 2, sticky ='N')
    DRWinroot.update_idletasks()
    MakeFParamFrame(FParamFrame, fname)

    PPFrame1 = tk.LabelFrame(DRWinnerFrame, padx=15)
    PPFrame1.grid(row=3,column=0, columnspan = 2, sticky ='N')
    MakePPFrame1(PPFrame1, SPs)

    plt.rcParams["figure.figsize"] = (12,8);
    plt.rcParams["font.size"] = 11;
    fig = plt.figure()                   
    
    
    DRWPlotFrame = tk.LabelFrame(DRWinnerFrame) 
    DRWPlotFrame.grid(row=4,column=0,columnspan = 2, sticky ='N')
    
    canvas = FigureCanvasTkAgg(fig,DRWPlotFrame); 
    canvas.get_tk_widget().grid(row=0, column=0, columnspan = 4, sticky='N')
    RD(SPs)
    canvas.draw()
            
    DRWinroot.mainloop()

