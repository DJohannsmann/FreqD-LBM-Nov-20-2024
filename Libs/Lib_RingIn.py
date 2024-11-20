import numpy as np
import time 
from scipy.ndimage import gaussian_filter
from Libs import Lib_General         as General
from Libs import Lib_Soft            as Soft
from Libs import Lib_OscBnd          as OscBnd
from Libs import Lib_Handle_Top      as Handle_Top
from Libs import Lib_FitRI           as FitRI
from Libs import Lib_StreamCollide   as StrColl
from Libs import Lib_IO              as IO
from Libs import Lib_Plots_from_Main as Plots_from_Main 
from Libs import Lib_Plots_for_GUI  as Plots_from_GUI 

def RingIn(SPs,FracVolSph,OscBndPars,\
        tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym,rhos,Do_Ref):
    nSph            = np.int64(SPs['nSph']) 
    nx              = np.int64(SPs['nx']) 
    ny              = np.int64(SPs['ny']) 
    nz              = np.int64(SPs['nz']) 
    n               = np.float64(SPs['n']) 
    om              = np.float64(SPs['om']) 
    tauInvBulk      = np.complex128(SPs['tauInvBulk'])
    UpdateMotionFac = np.float64(SPs['UpdateMotionFac'])
    SphPoss         = np.float64(SPs['SphPoss'])
    dimensions      = SPs['dimensions']
    Do_Plot_RingIns = SPs['Do_Plot_RingIns']
    if Do_Ref : Increased_Precision_Fac = 0.1
    else      : Increased_Precision_Fac = 1
    if Do_Ref : 
        ZBulk      = SPs['ZBulk']
        Lambda_TRT = SPs['Lambda_TRT']
        tauInv_Ref = 1./tauInvBulk
        if Lambda_TRT != 0 :  
            tauInv_Asym_Ref = (4.-2.*tauInv_Ref)/(2.-tauInv_Ref+4*Lambda_TRT*tauInv_Ref)
        else : tauInv_Asym_Ref = tauInv_Ref
        one_m_tauInv_m_Iom_Ref      = 1 - tauInv_Ref - 1j*om
        one_m_tauInv_m_Iom_Asym_Ref = 1 - tauInv_Ref - 1j*om

        one_m_tauInv_m_Iom_Ref      += 5./2.*om**2
        one_m_tauInv_m_Iom_Asym_Ref += 5./2.*om**2
        
    
    nd,cxs,cys,czs,ibars,wi,i_ups,i_notups,i_downs,i_notdowns = General.ReadStencil(SPs['dimensions'])      

    nuBulk = (1./tauInvBulk-0.5)/3.
    if dimensions == 1: MatricesTop = np.nan
    if dimensions == 2: MatricesTop = Handle_Top.Calc_MatricesTop_2D(nx,nuBulk,om)
    if dimensions == 3: MatricesTop = Handle_Top.Calc_MatricesTop_3D(nx,nz,nuBulk,om)
    FitInterval   = np.max([int(ny**2*0.05),10])
    PrintInterval = int(ny**2 * SPs['PrintIntervalFac'])
    PrintInterval = np.max([FitInterval,PrintInterval]) 
    FxCont_tot = 0
    if dimensions == 1: h = np.zeros((   ny,   nd),dtype = np.complex128)
    if dimensions == 2: h = np.zeros((nx,ny,   nd),dtype = np.complex128)
    if dimensions == 3: h = np.zeros((nx,ny,nz,nd),dtype = np.complex128)
    h1 = np.zeros((ny,nd),dtype = np.complex128)

    tbytRI_Extrapols = np.ones( 1000 * ny**2                 )*np.nan
    Dfcbyn_Extrapols = np.ones( 1000 * ny**2   ,dtype=complex)*np.nan
    tbytRI_RIs       = np.ones( 1000 * ny**2                 )*np.nan
    Dfcbyn_RIs       = np.ones( 1000 * ny**2   ,dtype=complex)*np.nan
    MotionPars_RI    = np.ones((1000 * ny**2,6),dtype=complex)*np.nan
    AuxPars_RI       = np.ones((1000 * ny**2,6),dtype=complex)*np.nan
    MotionParTitles = ['','','','','','']
    AuxParTitles    = ['','','','','','']

    time0=time.time(); countFits = 0; step = 0; Converged = False;
    if not Do_Ref and SPs['Do_OscBnd'] : 
        iSs,nLs,nLstot,xLs,yLs,zLs,OutsideLBMDomains,InParticles,i_BCs,PoiLs,\
            qs,uxLs,uyLs,uzLs,OscBndAmps,SphRespPars,UpdateMotionFac,\
            OscBndLocked,OscBndLockedTo,nSph,RSph,ySphbyR,rhoSph,iSiL_Lists = \
            OscBnd.Extract_OscBndPars(OscBndPars)
    B = Soft.compute_B(nx, ny, nz, int(SPs['nSph']), SPs['RSph'], SphPoss[0], SphPoss[1], SphPoss[2])
    while not Converged: 
        if Do_Ref : 
            h1,Fx_on_Wall = StrColl.FD_LBM_Step_Ref(h1,ny,nd,cxs,cys,czs,wi,ibars,\
                    i_ups,i_notups,i_downs,i_notdowns,ZBulk,om,\
                    tauInv_Ref,tauInv_Asym_Ref,one_m_tauInv_m_Iom_Ref,one_m_tauInv_m_Iom_Asym_Ref)
            MotionPars = np.ones(6,dtype = complex)*np.nan            
            AuxPars    = np.ones(6,dtype = complex)*np.nan            
        if not Do_Ref: 
            if SPs['Do_OscBnd']  : 
                if SPs['dimensions'] == 3: 
                    h,Fx_on_Wall,FxLs,FyLs,FzLs = \
                        StrColl.FreqDLBMStep_OscBnd_3D(h,nx,ny,nz,nd,cxs,cys,czs,wi,ibars,MatricesTop,\
                            nLstot,OutsideLBMDomains,i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
                            tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym)
                    uxLs,uyLs,uzLs,FxCont_tot,OscBndAmps,MotionPars,MotionParTitles,\
                        AuxPars,AuxParTitles = \
                        OscBnd.Update_Motion_3D(nLs,nLstot,iSs,xLs,yLs,zLs,uxLs,uyLs,uzLs,\
                            FxLs,FyLs,FzLs,nSph,om,\
                            OscBndAmps,SphPoss,SphRespPars,\
                            UpdateMotionFac,OscBndLocked,OscBndLockedTo,iSiL_Lists,n,SPs)
                if SPs['dimensions'] == 2: 
                    h,Fx_on_Wall,FxLs,FyLs = \
                        StrColl.FreqDLBMStep_OscBnd_2D(h,nx,ny,nd,cxs,cys,czs,wi,ibars,MatricesTop,\
                            nLstot,OutsideLBMDomains,i_BCs,PoiLs,qs,uxLs,uyLs,uzLs,\
                            tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym)
                    uxLs,uyLs,FxLiq,MotionPars,MotionParTitles,AuxPars,AuxParTitles = \
                        OscBnd.Update_Motion_2D(nLstot,FxLs)
            if not SPs['Do_OscBnd']  : 
                if SPs['ProblemType'] == 'SoftParticles' : 
                    h,Fx_on_Wall = \
                        StrColl.FreqDLBMStep_SoftPt_3D(h,nx,ny,nz,\
                            nd,cxs,cys,czs,wi,ibars,i_ups,i_notups,i_downs,i_notdowns,MatricesTop,\
                            tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym)
                if SPs['ProblemType'] == 'FilmResonance' : 
                    h,Fx_on_Wall = \
                        StrColl.FreqDLBMStep_SoftPt_1D(h,nx,ny,nz,nd,cxs,cys,czs,wi,ibars,\
                                i_ups,i_notups,i_downs,i_notdowns,om,MatricesTop,\
                                tauInvs,tauInvs_Asym,one_m_tauInvs_m_Iom,one_m_tauInvs_m_Iom_Asym)
                MotionPars = np.ones(6,dtype = complex)*np.nan            
                AuxPars    = np.ones(6,dtype = complex)*np.nan            
        DfcbynRaw = General.Calc_Dfcbyn(SPs,Fx_on_Wall,FxCont_tot,Do_Ref)
        if step >= len(tbytRI_RIs) - 1: 
            tbytRI_RIs       = np.append(tbytRI_RIs      ,np.ones( 1000 * ny**2                  )*np.nan)
            Dfcbyn_RIs       = np.append(Dfcbyn_RIs      ,np.ones( 1000 * ny**2   ,dtype= complex)*np.nan)
            tbytRI_Extrapols = np.append(tbytRI_Extrapols,np.ones( 1000 * ny**2                  )*np.nan)
            Dfcbyn_Extrapols = np.append(Dfcbyn_Extrapols,np.ones( 1000 * ny**2   ,dtype= complex)*np.nan)
            MotionPars_RI    = np.append(MotionPars_RI   ,np.ones((1000 * ny**2,6),dtype= complex)*np.nan)
            AuxPars_RI       = np.append(AuxPars_RI      ,np.ones((1000 * ny**2,6),dtype= complex)*np.nan)

        tbytRI_RIs[step] = float(step)/ny**2
        if Do_Ref : Dfcbyn_RIs[step] = DfcbynRaw
        else      : Dfcbyn_RIs[step] = DfcbynRaw - SPs['Dfcbyn_Ref']
        
        # with open(SPs['folder']+'\\errors1.txt', "a") as f:
        #     print('step               : ', step, file=f) 
        #     print('MotionPars         : ', MotionPars, file=f) 
        #     print('MotionPars_RI[step]: ', MotionPars_RI[step], file=f)             

        MotionPars_RI[step] = MotionPars 
     
        if step%FitInterval == 0 and step >= 10 : 
            if SPs['SigSmoothDfcbynsFac'] > 0 and step > 2*ny**2: 
                sig = step*SPs['SigSmoothDfcbynsFac']   
                i_ini =      int(np.max([10*sig,step/3.]))
                i_fin = step-int(10*sig)
                tbytRI_RIs4Fit = tbytRI_RIs[i_ini:i_fin]
                Dfcbyn_RIs4Fit = gaussian_filter(Dfcbyn_RIs,sigma=sig)
                Dfcbyn_RIs4Fit = Dfcbyn_RIs4Fit[i_ini:i_fin]
            else : 
                i_ini = int(step/3.)
                tbytRI_RIs4Fit = tbytRI_RIs[i_ini:step]
                Dfcbyn_RIs4Fit = Dfcbyn_RIs[i_ini:step]     
            try :     
                Dfcbyn_Extrapol,StdErr_Extrapol,amplitude,om_complex,Dfcbyn_RI_Fit = \
                    FitRI.Fit_RI(tbytRI_RIs4Fit,Dfcbyn_RIs4Fit)
            except : 
                Dfcbyn_Extrapol = np.nan
                amplitude       = np.nan
                om_complex      = np.nan
                Dfcbyn_RI_Fit   = np.nan
            countFits += 1    
            tbytRI_Extrapols[countFits] = step/ny**2
            Dfcbyn_Extrapols[countFits] = Dfcbyn_Extrapol
            
        if step%PrintInterval == 0 and step > 10:
            DriftFitResults40perc,DriftFitResults20perc = \
                FitRI.Calc_DriftFitResults(SPs,tbytRI_Extrapols,Dfcbyn_Extrapols,\
                    amplitude,om_complex,countFits)
            if not Do_Ref : 
                if Do_Plot_RingIns :  
                    if SPs['Do_from_GUI'] :  
                        Plots_from_GUI.Plot_RI(tbytRI_RIs4Fit, Dfcbyn_RIs4Fit, Dfcbyn_RI_Fit,tbytRI_Extrapols, Dfcbyn_Extrapols,countFits, DriftFitResults40perc, DriftFitResults20perc, SPs)      
                    else :  
                        Plots_from_Main.Plot_RI(tbytRI_RIs4Fit,Dfcbyn_RIs4Fit,Dfcbyn_RI_Fit,tbytRI_Extrapols,Dfcbyn_Extrapols,countFits)                    
                        
                if SPs['Do_Plot_MotionPars'] and SPs['ProblemType'] == 'StiffParticles':                                                          
                    if SPs['Do_from_GUI'] :  
                        Plots_from_GUI.Plot_MotionPars_RI( tbytRI_Extrapols[int(countFits/3.):countFits], MotionPars_RI[int(countFits/3.):countFits],MotionParTitles,SPs)
                    else :  
                        Plots_from_Main.Plot_MotionPars_RI(tbytRI_Extrapols[int(countFits/3.):countFits], MotionPars_RI[int(countFits/3.):countFits],MotionParTitles,SPs)

                print('t/t_RI',np.round(step/ny**2,2),\
                      'Drift40%',np.round(np.abs(DriftFitResults40perc),3),\
                      'Drift20%',np.round(np.abs(DriftFitResults20perc),3),)
            if np.abs(DriftFitResults40perc.real) < SPs['TargetSlopeFitResults']*Increased_Precision_Fac and \
               np.abs(DriftFitResults40perc.imag) < SPs['TargetSlopeFitResults']*Increased_Precision_Fac and \
               np.abs(DriftFitResults20perc.real) < SPs['TargetSlopeFitResults']*Increased_Precision_Fac and \
               np.abs(DriftFitResults20perc.imag) < SPs['TargetSlopeFitResults']*Increased_Precision_Fac : 
                   Converged = True; 
                   CompTimeMins = np.round((time.time()-time0)/60,2); 
                   SPs['ProblemFlag']     = 0
                   SPs['steps']           = step
                   SPs['tbytRI']          = step/ny**2
                   SPs['Dfcbyn_Extrapol'] = Dfcbyn_Extrapol
                   SPs['Dfratio']         = Dfcbyn_Extrapol.imag / (-Dfcbyn_Extrapol.real)
                   SPs['CompTimeMins'] = CompTimeMins
            if np.abs(Dfcbyn_Extrapol) > 1e7 or step/ny**2> SPs['MaxtbytRI'] : 
                CompTimeMins = np.round((time.time()-time0)/60,2);
                print('np.abs(Dfcbyn_Fit) > 1e7 or tbytRI > MaxtbytRI',SPs['MaxtbytRI'])
                #UPDATED
                if SPs['Do_from_GUI']: Plots_from_GUI.SimError(SPs)
                SPs['ProblemFlag']  = 1
                dr = np.ones(((nx,ny,nz)))*np.nan
                ux = np.ones(((nx,ny,nz)))*np.nan
                uy = np.ones(((nx,ny,nz)))*np.nan
                uz = np.ones(((nx,ny,nz)))*np.nan
                break
        step += 1    
        
    if Converged : 
        if Do_Ref :         
            SPs['Dfcbyn_Ref'] = Dfcbyn_Extrapol
            SPs['StdErr_Ref'] = StdErr_Extrapol
            SPs['Dfratio_Ref'] = Dfcbyn_Extrapol.imag/(-Dfcbyn_Extrapol.real)
            dr = np.ones(ny)*np.nan
            ux = np.sum([h]*cxs,axis=1)
            uy = np.ones(ny)*np.nan
            uz = np.ones(ny)*np.nan
        if not Do_Ref : 
            if SPs['ProblemType'] == 'SoftParticles' : 
                dr,ux,uy,uz = Soft.Calc_dr_ux_uy_uz_SoftPt_3D(h,cxs,cys,czs,nx,ny,nz)
                # MotionPars,MotionParTitles,AuxPars,AuxParTitles = \
                    # Soft.Calc_MotionPars_3D(nx,ny,nz,nd,cxs,cys,czs,i_ups,i_notups,ibars,wi,h,\
                    #                  tauInvs,FracVolSph,SPs,SphPoss)
                MotionPars,MotionParTitles,AuxPars,AuxParTitles = \
                    Soft.Calc_MotionPars_3D_UU(B,nx,ny,nz,nd,cxs,cys,czs,i_ups,i_notups,ibars,wi,h,\
                                     tauInvs,FracVolSph,SPs,SphPoss)
                MotionPars_RI[step] = MotionPars
                AuxPars_RI[   step] = AuxPars
            if SPs['ProblemType'] in ['StiffParticles','Roughness','SFA']: 
                dr,ux,uy,uz = \
                    OscBnd.Calc_dr_ux_uy_uz_OscBnd(h,cxs,cys,czs,\
                        OutsideLBMDomains,InParticles,OscBndAmps,SPs)
            if SPs['ProblemType'] == 'FilmResonance': 
                dr = np.zeros(ny)
                ux = np.ones(ny,dtype=np.complex128)*np.nan
                for y in range(ny): ux[y] = np.sum(h[y]*cxs)
                uy = np.zeros(ny)
                uz = np.zeros(ny)
        
            if SPs['ProblemType'] in ['SoftParticles','StiffParticles','Roughness','SFA'] :
                if SPs['Do_from_GUI'] :  
                    Plots_from_GUI.Plot_Fields_Horizontal(dr,ux,uy,uz,'Re($\Delta \\rho$)', 'Re(u$_{\mathrm{x}}$)','Re(u$_{\mathrm{y}}$)', 'Re(u$_{\mathrm{z}}$)',SPs,int(SPs['RSph']))
                    Plots_from_GUI.Plot_Fields_Vertical(dr,ux,uy,uz, 'Re($\Delta \\rho$)', 'Re(u$_{\mathrm{x}}$)','Re(u$_{\mathrm{y}}$)','Re(u$_{\mathrm{z}}$)',SPs)
                else : 
                    Plots_from_Main.Plot_Fields_Horizontal(dr,ux,uy,uz,'Re($\Delta \\rho$)', 'Re(u$_{\mathrm{x}}$)','Re(u$_{\mathrm{y}}$)', 'Re(u$_{\mathrm{z}}$)',SPs,int(SPs['RSph']))
                    Plots_from_Main.Plot_Fields_Vertical(dr,ux,uy,uz, 'Re($\Delta \\rho$)', 'Re(u$_{\mathrm{x}}$)','Re(u$_{\mathrm{y}}$)','Re(u$_{\mathrm{z}}$)',SPs)
                    


                # Plots.Plot_Top(ux,uy,uz,uz,\
                #     'Re(u$_{\mathrm{x}}$)',\
                #     'Re(u$_{\mathrm{y}}$)',\
                #     'Re(u$_{\mathrm{z}}$)','',SPs)


            if SPs['ProblemType'] == 'FilmResonance' :
                pass #Plots.Plot_DisplacementField_1D(h,SPs)

            MotionParsDict,AuxDict = IO.Make_MotionParsDict_AusParsDict(\
                MotionPars_RI[step-1],MotionParTitles,\
                AuxPars_RI[   step-1],AuxParTitles)
            print(SPs['iavg'],SPs['iPar1'],SPs['iPar2'],SPs['iPar3'],SPs['iovt'],\
                'Dfcbyn',np.round(SPs['Dfcbyn_Extrapol'],3),'CompTimeMins',SPs['CompTimeMins'])
            IO.Write_Config(SPs); IO.Save(SPs,MotionParsDict,AuxDict)    
    return   