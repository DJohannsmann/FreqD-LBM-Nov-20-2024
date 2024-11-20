import numpy as np

"""
When adding functions to this library, use the OO format of matplotlib (ax=fig. etc. instead of plt.)
These functions are called ONLY by the GUI. They read temporary npy files containing plot data and create plots in the GUI
The corresponding functions that create the files are in the Lip_Plots library. 
"""


def Plot_LinkProps_3D(fname, fig):
    #xLs,yLs,zLs,color1,title1
    """
    This function takes the values from an .npy file saved in the tmpplot folder
    and uses them to plot the spheres in the GUI frame in 3D.
    The files are created by the GUI interface during setup but also as the simulation is running
    whenever the sphere number, truncation, or arrangement is altered.
    This function is called only from the GUI interface.
    """
    #print('In Plot_LinkProps_3D ', fname)
    fig.clear()
    ax=fig.add_subplot(projection='3d')     
    ax.set_title('Spheres')  
    if fname!=' ':
        iParams={}
        try:        
            iParams = np.load(fname, allow_pickle='TRUE').item()
            xLs = iParams['xLs']
            yLs = iParams['yLs']
            zLs = iParams['zLs']
            color1 = iParams['color1']
            title1 = iParams['title1']
        except:
            print("Plot_LinkProps_3D No idea what happened here. File doesn't exist")
            return
        
        sc=ax.scatter(xLs,zLs,yLs,c=color1.real)
        ax.set_aspect('equal')        
        ax.figure.canvas.figure.colorbar(sc);
        ax.set_title(title1)  

    ax.axis("off")
    fig.tight_layout()
        #     if SPs['Do_SavePlots']:
        #         plt.savefig(str(SPs['iavg'])+'_'+\
        #                     str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
        #                     str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_LinkProps.png')
    #plt.close("all")



def Plot_Sph(fname, fig):
    #nx,ny,nz,FracVolSph,SPs
    """
    This function takes the values from an .npy file saved in the tmpplot folder
    and uses them to plot the spheres in the GUI frame in 2D.
    This function is called only from the GUI interface.
    UNUSED
    """
    #print('Plot_Sph')
    fig.clear()
    if fname!=' ':
        iParams={}
        try:        
            iParams = np.load(fname, allow_pickle='TRUE').item()
            nx=iParams['nx']
            ny=iParams['ny']
            nz=iParams['nz']
            FracVolSph=iParams['FracVolSph']
            SPs=iParams['SPs']
        except:
            print("Plot_Sph No idea what happened here. File doesn't exist")
            return
        
        SphPoss = SPs['SphPoss']
        yplot = int(SphPoss[1,0])
        X,Z = np.meshgrid(np.arange(nx), np.arange(nz))
        ax=fig.add_subplot()
        sc=ax.pcolormesh(X,Z,np.transpose(FracVolSph[:,yplot,:]))
        ax.figure.canvas.figure.colorbar(sc);
        ax.axis("off")
        ax.set_title('$\phi_{\mathrm{Sph}}$')  
        ax.axis("off")
        fig.tight_layout()

def Plot_RI(fname, fig):   
    #tbytRI_RIs4Fit,Dfcbyn_RIs4Fit,Dfcbyn_RI_Fit,tbytRI_Extrapols,Dfcbyn_Extrapols,countFits, SPs
    """
    This function takes the values from an .npy file saved in the tmpplot folder
    and uses them to plot the ring-in results in the GUI frame.
    The folder is assumed to exist.
    The files are created by the corresponding function in Lib_plots.py that is called as the simulation is running. 
    This function is called only from the GUI interface.
    """
    
    fig.clear()
    fig.suptitle('Ring In \n')    
    iParams={}
    
    if fname!=' ':
        try:        
            iParams = np.load(fname, allow_pickle='TRUE').item()

            tbytRI_RIs4Fit = iParams['tbytRI_RIs4Fit']
            Dfcbyn_RIs4Fit = iParams['Dfcbyn_RIs4Fit'] 
            Dfcbyn_RI_Fit = iParams['Dfcbyn_RI_Fit'] 
            tbytRI_Extrapols = iParams['tbytRI_Extrapols']
            Dfcbyn_Extrapols = iParams['Dfcbyn_Extrapols'] 
            countFits = iParams['countFits'] 
            DriftFitResults40perc= iParams['DriftFitResults40perc']
            DriftFitResults20perc = iParams['DriftFitResults20perc']
            SPs = iParams['SPs']
        except:
            print("Display_RingIn No idea what happened here. File doesn't exist")
            return
        
        ax1=fig.add_subplot(2,2,1)
        ax1.plot(tbytRI_RIs4Fit,Dfcbyn_RIs4Fit.real); 
        ax1.plot(tbytRI_RIs4Fit,Dfcbyn_RI_Fit.real)
        ax1.set_ylabel('$\Delta f/n$')
        ax1.set_title('From individual steps')
        ax1.set_xticks
        
        ax2=fig.add_subplot(2,2,3)
        ax2.plot(tbytRI_RIs4Fit,Dfcbyn_RIs4Fit.imag); 
        ax2.plot(tbytRI_RIs4Fit,Dfcbyn_RI_Fit.imag)
        ax2.set_xlabel('$t$ / $t_{RI}$'); 
        ax2.set_ylabel('$\Delta \Gamma/n$')
        
         
        ax3=fig.add_subplot(2,2,2)
        ax3.plot(tbytRI_Extrapols[int(countFits*2./3.):countFits],Dfcbyn_Extrapols[int(countFits*2./3.):countFits].real); 
        ax3.set_xticks([])
        
        if not np.isnan(DriftFitResults40perc) and not np.isnan(DriftFitResults20perc):
            s = 'Drift [Hz/$t_{\mathrm{RI}}$] 40%: ' + str("{:.2f}".format(float(np.abs(DriftFitResults40perc)))) +\
             ' 20%: '                              + str("{:.2f}".format(float(np.abs(DriftFitResults20perc)))) + \
             ' Target: ' + str(SPs['TargetSlopeFitResults'])
            s1 = str(SPs['iavg']+1) + ', ' + str(SPs['n']) + ', ' + str("{:.2f}".format(SPs[SPs['Par3str']])) + ', ' + str("{:.2f}".format(SPs[SPs['Par2str']])) + ', ' + str("{:.2f}".format(SPs[SPs['Par1str']]))
            fig.suptitle('Ring In:' + s1 + ' \n' + s)

        ax3.set_title('Extrapolated fit results')
        
        ax4=fig.add_subplot(2,2,4)
        ax4.plot(tbytRI_Extrapols[int(countFits*2./3.):countFits],Dfcbyn_Extrapols[int(countFits*2./3.):countFits].imag); 
        ax4.set_xlabel('$t$ / $t_{RI}$'); 
        
        fig.tight_layout()

    else:
        ax1=fig.add_subplot(2,2,1)
        ax1.set_ylabel('$\Delta f/n$')
        ax1.set_title('Individual steps')
        ax2=fig.add_subplot(2,2,3)
        ax2.set_xlabel('$t$ / $t_{RI}$'); 
        ax2.set_ylabel('$\Delta \Gamma/n$')
        ax3=fig.add_subplot(2,2,2)
        ax3.set_title('Extrapolated results')
        ax4=fig.add_subplot(2,2,4)
        ax4.set_xlabel('$t$ / $t_{RI}$'); 
        fig.tight_layout()


def Plot_MotionPars_RI(fname, fig):  
    #tbytHyd_RI,MotionPars_RI,MotionParsTitles,SPs
    """
    This function takes the values from an .npy file saved in the tmpplot folder
    and uses them to plot the motional parameters of the particles in the GUI frame.
    The folder is assumed to exist.
    The files are created by the corresponding function in Lib_plots.py that is called as the simulation is running. 
    This function is called only from the GUI interface.
    Only used for stiff particles?
    """
   
    fig.clear()
    fig.suptitle('Motion Pars')    
    iParams={}
    ax=[]    
    if fname!=' ':
        try:        
            iParams = np.load(fname, allow_pickle='TRUE').item()      

            tbytHyd_RI=iParams['tbytHyd_RI']
            MotionPars_RI=iParams['MotionPars_RI']
            MotionParsTitles=iParams['MotionParsTitles']
            SPs=iParams['SPs']
        except:
            print("Plot_MotionPars_RI No idea what happened here. File doesn't exist")
            return
    
        i_ini = int(len(MotionPars_RI)/3)
    
        for i_MotionPar in range(len(MotionPars_RI[0])):
            ax.append(fig.add_subplot(2,3,i_MotionPar+1))
            ax[i_MotionPar].plot(tbytHyd_RI[i_ini:],MotionPars_RI[i_ini:,i_MotionPar].real,label='Re')
            ax[i_MotionPar].plot(tbytHyd_RI[i_ini:],MotionPars_RI[i_ini:,i_MotionPar].imag,label='Im')
            ax[i_MotionPar].set_xlabel('$t$ / $t_{RI}$'); 
            ax[i_MotionPar].set_title(MotionParsTitles[i_MotionPar]) 
            if i_MotionPar == 0 : ax[i_MotionPar].legend(fontsize = 8)
        fig.tight_layout()
        s1 = str(SPs['iavg']+1) + ', ' + str(SPs['n']) + ', ' + str("{:.2f}".format(SPs[SPs['Par3str']])) + ', ' + str("{:.2f}".format(SPs[SPs['Par2str']])) + ', ' + str("{:.2f}".format(SPs[SPs['Par1str']]))
        fig.suptitle('Motion Pars: ' + s1)
    else:
        for i in range(6):
            ax.append(fig.add_subplot(2,3,i+1))
            ax[i].set_xlabel('$t$ / $t_{RI}$'); 
            ax[i].set_title(' ') 

        fig.tight_layout()
    
    # if SPs['Do_SavePlots']:
    #     plt.savefig(str(SPs['iavg'])+'_'+\
    #                 str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
    #                 str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'_MotionPars.png')
    

    
def Plot_Fields_Horizontal(fname, fig):
    # val1,val2,val3,val4,title1,title2,title3,title4,SPs,yplot   
    """
    This function takes the values from an .npy file saved in the tmpplot folder
    and uses them to plot horizontal fields of the particles in the GUI frame.
    The folder is assumed to exist.
    The files are created by the corresponding function in Lib_plots.py that is called as the simulation is running. 
    This function is called only from the GUI interface.
    """

    fig.clear()

    s_title='Horizontal cut, y ='  
    s1 = ''
    iParams={}
    
    if fname!=' ':
        try:        
            iParams = np.load(fname, allow_pickle='TRUE').item()      
            val1 = iParams['dr']
            val2 = iParams['ux']
            val3 = iParams['uy']
            val4 = iParams['uz']
            title1 = iParams['title1']
            title2 = iParams['title2'] 
            title3 = iParams['title3'] 
            title4=  iParams['title4'] 
            SPs = iParams['SPs'] 
            s1 = str(SPs['iavg']+1) + ', ' + str(SPs['n']) + ', ' + str("{:.2f}".format(SPs[SPs['Par3str']])) + ', ' + str("{:.2f}".format(SPs[SPs['Par2str']])) + ', ' + str("{:.2f}".format(SPs[SPs['Par1str']]))
            yplot = int(iParams['RSph'])
            
            if yplot > SPs['ny']-1: yplot = SPs['ny']-1 
            
            
            X,Z = np.meshgrid(np.arange(int(SPs['nx'])), np.arange(int(SPs['nz'])))
            if title1 != '':
                ax1 = fig.add_subplot(1,4,1); 
                sc = ax1.pcolormesh(X,Z, np.transpose(val1[:,yplot,:].real)) 
                ax1.set_title(title1); 
                ax1.set_xlabel('x'); 
                ax1.set_ylabel('z'); 
                ax1.set_xticks([]); 
                ax1.set_yticks([]); 
                fig.colorbar(sc)
            if title2 != '':
                ax2 = fig.add_subplot(1,4,2); 
                sc = ax2.pcolormesh(X,Z, np.transpose(val2[:,yplot,:].real)) 
                ax2.set_title(title2); 
                ax2.set_xlabel('x'); 
                ax2.set_ylabel('z'); 
                ax2.set_xticks([]); 
                ax2.set_yticks([]); 
                fig.colorbar(sc)
            if title3 != '':
                ax3 = fig.add_subplot(1,4,3); 
                sc = ax3.pcolormesh(X,Z, np.transpose(val3[:,yplot,:].real))
                ax3.set_title(title3); 
                ax3.set_xlabel('x'); 
                ax3.set_ylabel('z'); 
                ax3.set_xticks([]); 
                ax3.set_yticks([]); 
                fig.colorbar(sc)
            if title4 != '':
                ax4 = fig.add_subplot(1,4,4); 
                sc = ax4.pcolormesh(X,Z, np.transpose(val4[:,yplot,:].real))
                ax4.set_title(title4); 
                ax4.set_xlabel('x'); 
                ax4.set_ylabel('z'); 
                ax4.set_xticks([]); 
                ax4.set_yticks([]); 
                fig.colorbar(sc)
            s_title='Horizontal cut, y = ' + str(yplot)
             
        except:
            print("Plot_Fields_Horizontal No idea what happened here. File doesn't exist")
            return
    else:
        ax1 = fig.add_subplot(1,4,1); 
        ax2 = fig.add_subplot(1,4,2); 
        ax3 = fig.add_subplot(1,4,3); 
        ax4 = fig.add_subplot(1,4,4); 
    
    # if SPs['Do_SavePlots']:
    #     plt.savefig(str(SPs['iavg']) +'_'+str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
    #                 str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'Hor.png')
    # plt.show()

    
    fig.suptitle(s_title + '\n' + s1)    
    fig.tight_layout()




def Plot_Fields_Vertical(fname, fig):
    #val1,val2,val3,val4,title1,title2,title3,title4,SPs
    """
    This function takes the values from an .npy file saved in the tmpplot folder
    and uses them to plot vertical fields of the particles in the GUI frame.
    The folder is assumed to exist.
    The files are created by the corresponding function in Lib_plots.py that is called as the simulation is running. 
    This function is called only from the GUI interface.
    """
    
    #print('In Plot_FieldsV', fname)
    fig.clear()
    s_title='Vertical cut, z =' 
    s1 = ''
    
    iParams={}
    
    if fname!=' ':
        try:        
            iParams = np.load(fname, allow_pickle='TRUE').item()      
            val1 = iParams['dr']
            val2 = iParams['ux']
            val3 = iParams['uy']
            val4 = iParams['uz']
            title1 = iParams['title1']
            title2 = iParams['title2'] 
            title3 = iParams['title3'] 
            title4=  iParams['title4'] 
            SPs = iParams['SPs'] 
            s1 = str(SPs['iavg']+1) + ', ' + str(SPs['n']) + ', ' + str("{:.2f}".format(SPs[SPs['Par3str']])) + ', ' + str("{:.2f}".format(SPs[SPs['Par2str']])) + ', ' + str("{:.2f}".format(SPs[SPs['Par1str']]))
        
            zmid = int(SPs['nz']/2)
            X,Y = np.meshgrid(np.arange(SPs['nx']), np.arange(SPs['ny']))
            if title1 != '':
                ax1 = fig.add_subplot(1,4,1);
                sc = ax1.pcolormesh(X,Y, np.transpose(val1[:,:,zmid].real)) 
                ax1.set_title(title1); 
                ax1.set_xlabel('x'); 
                ax1.set_ylabel('y'); 
                ax1.set_xticks([]); 
                ax1.set_yticks([]); 
                fig.colorbar(sc)
            if title2 != '':
                ax2 = fig.add_subplot(1,4,2); 
                sc = ax2.pcolormesh(X,Y, np.transpose(val2[:,:,zmid].real)) 
                ax2.set_title(title2); 
                ax2.set_xlabel('x'); 
                ax2.set_ylabel('y'); 
                ax2.set_xticks([]); 
                ax2.set_yticks([]); 
                fig.colorbar(sc)
            if title3 != '':
                ax3 = fig.add_subplot(1,4,3); 
                sc = ax3.pcolormesh(X,Y, np.transpose(val3[:,:,zmid].real))
                ax3.set_title(title3); 
                ax3.set_xlabel('x'); 
                ax3.set_ylabel('y'); 
                ax3.set_xticks([]); 
                ax3.set_yticks([]); 
                fig.colorbar(sc)
            if title4 != '':
                ax4 = fig.add_subplot(1,4,4); 
                sc = ax4.pcolormesh(X,Y, np.transpose(val4[:,:,zmid].real))
                ax4.set_title(title4); 
                ax4.set_xlabel('x'); 
                ax4.set_ylabel('y'); 
                ax4.set_xticks([]); 
                ax4.set_yticks([]); 
                fig.colorbar(sc)
            s_title='Vertical cut, z = ' + str(zmid)
            
        except:
            print("Plot_Fields_Vertical No idea what happened here. File doesn't exist")
            return
    else:
        ax1 = fig.add_subplot(1,4,1); 
        ax2 = fig.add_subplot(1,4,2); 
        ax3 = fig.add_subplot(1,4,3); 
        ax4 = fig.add_subplot(1,4,4); 

    fig.suptitle(s_title + '\n' + s1)
    fig.tight_layout()
    # if SPs['Do_SavePlots']:
    #     plt.savefig(str(SPs['iavg'])+'_'+str(SPs['iPar1'])+'_'+str(SPs['iPar2'])+'_'+\
    #         str(SPs['iPar3'])+'_'+str(SPs['iovt'])+'Ver.png')



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