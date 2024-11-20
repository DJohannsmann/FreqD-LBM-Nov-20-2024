import numpy as np
import matplotlib.pyplot as plt
from Libs import Lib_General as General

def Plot_Sph(nx,ny,nz,FracVolSph,SPs):
    SphPoss = SPs['SphPoss']
    yplot = int(SphPoss[1,0])
    plt.rcParams["figure.figsize"] = (1.5,1.3)  
    plt.rcParams["font.size"] = 8;
    X,Z = np.meshgrid(np.arange(nx), np.arange(nz))
    plt.pcolormesh(X,Z,np.transpose(FracVolSph[:,yplot,:])); plt.axis('off') 
    plt.colorbar(); plt.title('$\phi_{\mathrm{Sph}}$'); plt.show()

def Plot_Fields_Horizontal(val1,val2,val3,val4,title1,title2,title3,title4,SPs,yplot):
    if yplot > SPs['ny']-1: yplot = SPs['ny']-1
    nplots = 0
    if title1 != '': nplots += 1
    if title2 != '': nplots += 1
    if title3 != '': nplots += 1
    if title4 != '': nplots += 1
    plt.rcParams["figure.figsize"] = (2.5*nplots,2.3)
    plt.rcParams["font.size"] = 10;
    X,Z = np.meshgrid(np.arange(int(SPs['nx'])), np.arange(int(SPs['nz'])))
    iplot = 1
    if title1 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Z, np.transpose(val1[:,yplot,:].real)) 
        plt.title(title1); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar()
        iplot += 1
    if title2 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Z, np.transpose(val2[:,yplot,:].real)) 
        plt.title(title2); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar()
        iplot += 1
    if title3 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Z, np.transpose(val3[:,yplot,:].real))
        plt.title(title3); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar()
        iplot += 1
    if title4 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Z, np.transpose(val4[:,yplot,:].real))
        plt.title(title4); plt.xlabel('x'); plt.ylabel('z'); plt.xticks([]); plt.yticks([]); plt.colorbar() 
    plt.suptitle('Horizontal cut, y = ' + str(yplot))    
    plt.tight_layout(); 
    if SPs['Do_SavePlots']:
        plt.savefig(str(SPs['iavg'])+'_'+\
                    str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
                    str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_Horizontal.png')
    plt.show()

def Plot_Fields_Vertical(val1,val2,val3,val4,title1,title2,title3,title4,SPs):
    nplots = 0
    if title1 != '': nplots += 1
    if title2 != '': nplots += 1
    if title3 != '': nplots += 1
    if title4 != '': nplots += 1
    plt.rcParams["figure.figsize"] = (2.5*nplots,2.3)
    plt.rcParams["font.size"] = 10;
    zmid = int(SPs['nz']/2)
    X,Y = np.meshgrid(np.arange(SPs['nx']), np.arange(SPs['ny']))
    iplot = 1
    if title1 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Y, np.transpose(val1[:,:,zmid].real)) 
        plt.title(title1); plt.xlabel('x'); plt.ylabel('y'); plt.xticks([]); plt.yticks([]); plt.colorbar()
        iplot += 1
    if title2 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Y, np.transpose(val2[:,:,zmid].real)) 
        plt.title(title2); plt.xlabel('x'); plt.ylabel('y'); plt.xticks([]); plt.yticks([]); plt.colorbar()
        iplot += 1
    if title3 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Y, np.transpose(val3[:,:,zmid].real))
        plt.title(title3); plt.xlabel('x'); plt.ylabel('y'); plt.xticks([]); plt.yticks([]); plt.colorbar()
        iplot += 1
    if title4 != '':
        plt.subplot(1,nplots,iplot); plt.pcolormesh(X,Y, np.transpose(val4[:,:,zmid].real))
        plt.title(title4); plt.xlabel('x'); plt.ylabel('y'); plt.xticks([]); plt.yticks([]); plt.colorbar()
    plt.suptitle('Vertical cut, z = ' + str(zmid))    
    plt.tight_layout()
    if SPs['Do_SavePlots']:
        plt.savefig(str(SPs['iavg'])+'_'+\
                    str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
                    str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_Vertical.png')
    plt.show()

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

def Plot_LinkProps_3D(xLs,yLs,zLs,color1,color2,color3,title1,title2,title3,SPs):
    nplots = 0
    if title1 != '': nplots += 1
    if title2 != '': nplots += 1
    if title3 != '': nplots += 1
    plt.rcParams["figure.figsize"] = (7.5*nplots,1.5)
    plt.rcParams['lines.markersize'] = 0.3
    plt.rcParams["font.size"] = 8;
    fig = plt.figure()
    iplot = 1
    if title1 != '':
        ax = fig.add_subplot(1,nplots,iplot,projection='3d')
        sc = ax.scatter(xLs,zLs,yLs,c=color1.real); plt.axis('off'); #plt.colorbar(sc); plt.title(title1)  
        iplot += 1
    if title2 != '':
        ax = fig.add_subplot(1,nplots,iplot,projection='3d')
        sc = ax.scatter(xLs,zLs,yLs,c=color2.real); plt.axis('off'); plt.colorbar(sc); plt.title(title2)  
        iplot += 1
    if title3 != '':
        ax = fig.add_subplot(1,nplots,iplot,projection='3d')
        sc = ax.scatter(xLs,zLs,yLs,c=color3.real); plt.axis('off'); plt.colorbar(sc); plt.title(title3)  
    plt.tight_layout(); 
    if SPs['Do_SavePlots']:
        plt.savefig(str(SPs['iavg'])+'_'+\
                    str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
                    str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_LinkProps.png')
    plt.show()

def Plot_LinkProps_Roughness_2D(xLs,yLs,color1,color2,color3,title1,title2,title3,SPs):
    nplots = 0
    if title1 != '': nplots += 1
    if title2 != '': nplots += 1
    if title3 != '': nplots += 1
    plt.rcParams["figure.figsize"] = (2.5*nplots,2.3)
    plt.rcParams['lines.markersize'] = 0.3
    plt.rcParams["font.size"] = 8;
    fig = plt.figure()
    iplot = 1
    if title1 != '':
        ax = fig.add_subplot(1,nplots,iplot)
        sc = ax.scatter(xLs,yLs,c=color1.real); plt.axis('off'); plt.colorbar(sc); plt.title(title1)  
        iplot += 1
    if title2 != '':
        ax = fig.add_subplot(1,nplots,iplot)
        sc = ax.scatter(xLs,yLs,c=color2.real); plt.axis('off'); plt.colorbar(sc); plt.title(title2)  
        iplot += 1
    if title3 != '':
        ax = fig.add_subplot(1,nplots,iplot)
        sc = ax.scatter(xLs,yLs,c=color3.real); plt.axis('off'); plt.colorbar(sc); plt.title(title3)  
    plt.tight_layout(); 
    if SPs['Do_SavePlots']:
        plt.savefig(str(SPs['iavg'])+'_'+\
                    str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
                    str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_LinkProps_Roughness_2D.png')
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

def Plot_MotionPars_RI(tbytHyd_RI,MotionPars_RI,MotionParsTitles,SPs):
    i_ini = int(len(MotionPars_RI)/3)
    plt.rcParams['figure.figsize'] = (7,2.5)
    plt.rcParams['font.size'] = 9
    for i_MotionPar in range(len(MotionPars_RI[0])):
        plt.subplot(2,6,i_MotionPar+1) 
        plt.plot(tbytHyd_RI[i_ini:],MotionPars_RI[i_ini:,i_MotionPar].real,label='Re')
        plt.plot(tbytHyd_RI[i_ini:],MotionPars_RI[i_ini:,i_MotionPar].imag,label='Im')
        plt.xlabel('$t$ / $t_{\mathrm{RI}}$',fontsize = 10); 
        plt.title(MotionParsTitles[i_MotionPar]) 
        if i_MotionPar == 0 : plt.legend(fontsize = 8)
    plt.tight_layout()
    if SPs['Do_SavePlots']:
        plt.savefig(str(SPs['iavg'])+'_'+\
                    str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
                    str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_MotionPars.png')
    plt.show()  
    
def Plot_RI(tbytRI_RIs4Fit,Dfcbyn_RIs4Fit,Dfcbyn_RI_Fit,\
            tbytRI_Extrapols,Dfcbyn_Extrapols,countFits):
    plt.rcParams['figure.figsize'] = (5,1.7); 
    plt.rcParams['font.size'] = 8

    plt.subplot(221); 
    plt.plot(tbytRI_RIs4Fit,Dfcbyn_RIs4Fit.real); 
    plt.plot(tbytRI_RIs4Fit,Dfcbyn_RI_Fit.real)
    plt.ylabel('$\Delta f/n$')
    plt.title('from individual steps')
    plt.xticks([])

    plt.subplot(223); 
    plt.plot(tbytRI_RIs4Fit,Dfcbyn_RIs4Fit.imag); 
    plt.plot(tbytRI_RIs4Fit,Dfcbyn_RI_Fit.imag)
    plt.xlabel('$t$ / $t_{\mathrm{RI}}$',fontsize = 10); 
    plt.ylabel('$\Delta \Gamma/n$')
 
    plt.subplot(222) 
    plt.plot(tbytRI_Extrapols[int(countFits*2./3.):countFits],\
             Dfcbyn_Extrapols[int(countFits*2./3.):countFits].real); 
    plt.title('extrapolated, fit results')
    plt.xticks([])

    plt.subplot(224) 
    plt.plot(tbytRI_Extrapols[int(countFits*2./3.):countFits],\
             Dfcbyn_Extrapols[int(countFits*2./3.):countFits].imag); 
    plt.xlabel('$t$ / $t_{\mathrm{RI}}$',fontsize = 10); 
    plt.tight_layout(); plt.show()
