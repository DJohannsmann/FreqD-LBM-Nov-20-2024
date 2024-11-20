import numpy as np
from datetime import datetime

def SimError(SPs):
    
    pfolder = SPs['folder'] + '\\tmpplot'
    fname = 'SimError'
    now = datetime.now()
    formatted = now.strftime("%Y-%m-%d-%H_%M_%S")
    fname=formatted+'_'+fname    
    s = '_'
    try:
        s = s + str(SPs['iavg']+1) + '_' + str(SPs['n']) + '_' + str("{:.2f}".format(SPs[SPs['Par3str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par2str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par1str']]))
    except:
        pass
    fname = fname + s
    print('Simulation error, creating: ', fname)  
    
    tempData={}
    try:
        tempData['SPs'] = SPs
    except:
        print('Error in append')
    
    np.save(pfolder +'\\'+ fname, tempData)  


def Plot_LinkProps_3D(xLs,yLs,zLs,color1,title1,SPs):
    """
    This function saves the plot data for the 3D spheres plot in an .npy file in the tmpplot folder.
    The folder is assumed to exist.
    """
    
    pfolder = SPs['folder'] + '\\tmpplot'
    fname = 'Plot_LinkProps_3D'
    now = datetime.now()
    formatted = now.strftime("%Y-%m-%d-%H_%M_%S")
    fname=formatted+'_'+fname    
    s = '_'
    try:
        s = s + str(SPs['iavg']+1) + '_' + str(SPs['n']) + '_' + str("{:.2f}".format(SPs[SPs['Par3str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par2str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par1str']]))
    except:
        pass
    fname = fname + s
    
    tempData={}
    try:
        tempData['xLs'] = xLs
        tempData['yLs'] = yLs
        tempData['zLs'] = zLs
        tempData['color1'] = yLs*SPs['Dx_nm']
        tempData['title1'] = ' '
    except:
        print('Error in append')
    
    np.save(pfolder +'\\'+ fname, tempData)  
    
    
    
def Plot_Sph(nx,ny,nz,FracVolSph,SPs):
    """
    This function saves the plot data for the 2D spheres plot in an .npy file in the tmpplot folder.
    The folder is assumed to exist.
    UNUSED
    I saw that the plot is called at some point. Not sure where.
    """    
    pfolder = SPs['folder'] + '\\tmpplot'
    fname = 'Plot_Sphere'
    now = datetime.now()
    formatted = now.strftime("%Y-%m-%d-%H_%M_%S")
    fname=formatted+'_'+fname
    s = '_'
    try:
        s = s + str(SPs['iavg']+1) + '_' + str(SPs['n']) + '_' + str("{:.2f}".format(SPs[SPs['Par3str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par2str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par1str']]))
    except:
        pass
    fname = fname + s
    
    tempData={}
    try: 
        tempData['nx'] = nx
        tempData['ny'] = ny
        tempData['nz'] = nz
        tempData['FracVolSph'] = FracVolSph
        tempData['SPs'] = SPs

    except:
        print('Error in append')

    np.save(pfolder +'\\'+ fname, tempData) 


def Plot_RI(tbytRI_RIs4Fit,Dfcbyn_RIs4Fit,Dfcbyn_RI_Fit,tbytRI_Extrapols,Dfcbyn_Extrapols,countFits, DriftFitResults40perc,DriftFitResults20perc, SPs):
    """
    This function saves the plot data for the ring-in plots in an .npy file in the tmpplot folder.
    The folder is assumed to exist.
    Change from the original version: SPs is added as a paramter, needed for folder access.
    SPs are not saved as they are not needed for the plot itself.
    """
    pfolder = SPs['folder'] + '\\tmpplot'
    fname = 'Plot_RI'
    now = datetime.now()
    formatted = now.strftime("%Y-%m-%d-%H_%M_%S")
    fname=formatted+'_'+fname
    s = '_'
    try:
        s = s + str(SPs['iavg']+1) + '_' + str(SPs['n']) + '_' + str("{:.2f}".format(SPs[SPs['Par3str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par2str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par1str']]))
    except:
        pass
    fname = fname + s

    
    tempData={}
    try: 
        tempData['tbytRI_RIs4Fit'] = tbytRI_RIs4Fit
        tempData['Dfcbyn_RIs4Fit'] = Dfcbyn_RIs4Fit
        tempData['Dfcbyn_RI_Fit'] = Dfcbyn_RI_Fit
        tempData['tbytRI_Extrapols'] = tbytRI_Extrapols
        tempData['Dfcbyn_Extrapols'] = Dfcbyn_Extrapols
        tempData['countFits'] = countFits    
        tempData['DriftFitResults40perc'] = DriftFitResults40perc
        tempData['DriftFitResults20perc'] = DriftFitResults20perc
        tempData['SPs'] = SPs
    except:
        print('Error in append')

    np.save(pfolder +'\\'+ fname, tempData)


def Plot_MotionPars_RI(tbytHyd_RI,MotionPars_RI,MotionParsTitles,SPs):
    """
    This function saves the plot data for the motional paramters in plots in an .npy file in the tmpplot folder.
    The folder is assumed to exist.
    """

    pfolder = SPs['folder'] + '\\tmpplot'
    fname = 'Plot_MotionPars_RI'
    now = datetime.now()
    formatted = now.strftime("%Y-%m-%d-%H_%M_%S")
    fname=formatted+'_'+fname
    s = '_'
    try:
        s = s + str(SPs['iavg']+1) + '_' + str(SPs['n']) + '_' + str("{:.2f}".format(SPs[SPs['Par3str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par2str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par1str']]))
    except:
        pass
    fname = fname + s
    tempData={}
    try: 
        tempData['tbytHyd_RI']=tbytHyd_RI
        tempData['MotionPars_RI']= MotionPars_RI
        tempData['MotionParsTitles']=MotionParsTitles
        tempData['SPs']=SPs
    except:
        print('Error in append')
    np.save(pfolder +'\\'+ fname, tempData)    
    
    

    
def Plot_Fields_Horizontal(val1,val2,val3,val4,title1,title2,title3,title4,SPs,yplot):
#Plots.Plot_Fields_Horizontal(dr,ux,uy,uz,'Re($\Delta \\rho$)', 'Re(u$_{\mathrm{x}}$)','Re(u$_{\mathrm{y}}$)', 'Re(u$_{\mathrm{z}}$)',SPs,int(SPs['RSph']))
    """
    This function saves the plot data for the horizontal fields plots in plots in an .npy file in the tmpplot folder.
    The folder is assumed to exist.
    """
    pfolder = SPs['folder'] + '\\tmpplot'
    fname = 'Plot_Fields_Horizontal'
    now = datetime.now()
    formatted = now.strftime("%Y-%m-%d-%H_%M_%S")
    fname=formatted+'_'+fname
    s = '_'
    try:
        s = s + str(SPs['iavg']+1) + '_' + str(SPs['n']) + '_' + str("{:.2f}".format(SPs[SPs['Par3str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par2str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par1str']]))
    except:
        pass
    fname = fname + s
    
    tempData={}
    try: 
        tempData['dr'] = val1
        tempData['ux'] = val2
        tempData['uy'] = val3
        tempData['uz'] = val4
        tempData['title1'] = title1
        tempData['title2'] = title2
        tempData['title3'] = title3
        tempData['title4'] = title4
        tempData['SPs'] = SPs                                
        tempData['RSph'] = SPs['RSph']                                                
    except:
        print('Error in append')

    np.save(pfolder +'\\'+ fname, tempData) 

def Plot_Fields_Vertical(val1,val2,val3,val4,title1,title2,title3,title4,SPs):
    """
    This function saves the plot data for the vertical fields plots in plots in an .npy file in the tmpplot folder.
    The folder is assumed to exist.
    """
    pfolder = SPs['folder'] + '\\tmpplot'
    fname = 'Plot_Fields_Vertical'
    now = datetime.now()
    formatted = now.strftime("%Y-%m-%d-%H_%M_%S")
    fname=formatted+'_'+fname
    s = '_'
    try:
        s = s + str(SPs['iavg']+1) + '_' + str(SPs['n']) + '_' + str("{:.2f}".format(SPs[SPs['Par3str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par2str']])) + '_' + str("{:.2f}".format(SPs[SPs['Par1str']]))
    except:
        pass
    fname = fname + s    
                    
    tempData={}
    try: 
        tempData['dr'] = val1
        tempData['ux'] = val2
        tempData['uy'] = val3
        tempData['uz'] = val4
        tempData['title1'] = title1
        tempData['title2'] = title2
        tempData['title3'] = title3
        tempData['title4'] = title4
        tempData['SPs'] = SPs                                
    except:
        print('Error in append')
    np.save(pfolder +'\\'+ fname, tempData) 
    
    
"""
def Plot_Top(val1,val2,val3,val4,title1,title2,title3,title4,SPs):
    nplots = 0
    if title1 != '': nplots += 1
    if title2 != '': nplots += 1
    if title3 != '': nplots += 1
    if title4 != '': nplots += 1
    plt.rcParams["figure.figsize"] = (2.5*nplots,2.3)
    plt.rcParams["font.size"] = 10;
    X,Z = np.meshgrid(np.arange(SPs['nx']), np.arange(SPs['nz']))
    iplot = 1
    if title1 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Z, np.transpose(val1[:,-1,:].real)) 
        plt.title(title1); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar()
        iplot += 1
    if title2 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Z, np.transpose(val2[:,-1,:].real)) 
        plt.title(title2); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar()
        iplot += 1
    if title3 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Z, np.transpose(val3[:,-1,:].real))
        plt.title(title3); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar()
        iplot += 1
    if title4 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Z, np.transpose(val4[:,-1,:].real))
        plt.title(title4); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar()
    plt.suptitle('Top')    
    plt.tight_layout()
    if SPs['Do_SavePlots']:
        plt.savefig(str(SPs['iavg'])+'_'+\
                    str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
                    str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_Top.png')
    plt.show()


def Plot_Fx_on_Wall(Fx_on_Wall,SPs):
    nx = int(SPs['nx'])
    nz = int(SPs['nz'])
    plt.rcParams["figure.figsize"] = (9,2.3)
    plt.rcParams["font.size"] = 10;
    X,Z = np.meshgrid(np.arange(nx), np.arange(nz))
    plt.subplot(141); plt.pcolormesh(X,Z,np.transpose(Fx_on_Wall.real)) 
    plt.title('Re(F$_{x}^{S})$'); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar()
    plt.subplot(142); plt.pcolormesh(X,Z,np.transpose(Fx_on_Wall.imag)) 
    plt.title('Im(F$_{x}^{S})$'); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar()
    plt.suptitle('transverse stress at resonator surface')    
    plt.tight_layout()
    if SPs['Do_SavePlots']:
        plt.savefig(str(SPs['iavg'])+'_'+\
                    str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
                    str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_Fx_on_Wall.png')
    plt.show()

def Plot_hs(h,SPs):
    if SPs['dimensions'] == 2: nd = 9 
    if SPs['dimensions'] == 3: nd = 19 
    nx = int(SPs['nx'])
    nz = int(SPs['nz'])
    plt.rcParams["figure.figsize"] = (9,8)
    X,Z = np.meshgrid(np.arange(nx), np.arange(nz))
    for i in range(nd):
        plt.subplot(6,7,2*i  +1); plt.pcolormesh(X,Z, np.transpose(h[:,0,:,i].real)); plt.axis('off') 
        plt.text(nx/2,nz/2,str(i),color = 'w'); 
        plt.subplot(6,7,2*i+1+1); plt.pcolormesh(X,Z, np.transpose(h[:,0,:,i].imag)); plt.axis('off') 
    plt.tight_layout(); 
    if SPs['Do_SavePlots']:
        plt.savefig(str(SPs['iavg'])+'_'+\
                    str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
                    str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_hs.png')
    plt.show()

def Plot_DisplacementField_1D(h,SPs):
    nd,cxs,cys,czs,ibars,wi,i_ups,i_notups,i_downs,i_notdowns = General.ReadStencil(SPs['dimensions'])
    ny = SPs['ny']
    y_nm = np.linspace(0,ny,ny)*SPs['Dx_nm']
    ux = np.ones(ny,dtype=np.complex128)*np.nan
    for y in range(ny): ux[y] = np.sum(h[y]*cxs)

    plt.rcParams['figure.figsize'] = (3,1.6); plt.rcParams['font.size'] = 8
    plt.plot(y_nm,ux.real,label = 'real',color = 'r'); 
    plt.plot(y_nm,ux.imag,label = 'imag',color = 'b'); 
    plt.xlabel('y [nm]'); plt.ylabel('displacement')
    plt.legend(loc = 'upper right'); 
    if SPs['Do_SavePlots']:
        plt.savefig(str(SPs['iavg'])+'_'+\
                    str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
                    str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_DisplacementField_1D.png')
    plt.show()

"""