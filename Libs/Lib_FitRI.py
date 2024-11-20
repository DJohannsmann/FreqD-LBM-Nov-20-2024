import numpy as np
import lmfit

def Fit_RI(tbytRI_RI,Dfcbyn_RI):
    nPars = 6
    ParNames = ['offset_Re','offset_Im','ampl_Re','ampl_Im','omc_Re','omc_Im']
    def Spiral(t,SpiralPars):
        offset     = SpiralPars[0]+1j*SpiralPars[1]
        amplitude  = SpiralPars[2]+1j*SpiralPars[3]
        om_complex = SpiralPars[4]+1j*SpiralPars[5]
        return offset+amplitude*np.exp(1j*om_complex*t)

    def residuals_lmfit(Pars):
        SpiralPars = np.ones(6)*np.nan
        for iPar in range(nPars): SpiralPars[iPar] = Pars[ParNames[iPar]].value
        Dfcbyn_RI_stacked = np.hstack((Dfcbyn_RI.real, Dfcbyn_RI.imag))
        Dfcbyn_RI_model = Spiral(tbytRI_RI,SpiralPars)
        Dfcbyn_RI_model_stacked = np.hstack((Dfcbyn_RI_model.real,Dfcbyn_RI_model.imag))
        return Dfcbyn_RI_stacked-Dfcbyn_RI_model_stacked

    def configure_Minimizer():
        p = lmfit.Parameters()
        GuessPars = np.array([Dfcbyn_RI[-1].real,
                              Dfcbyn_RI[-1].imag,
                              Dfcbyn_RI[ 0].real-Dfcbyn_RI[-1].real,
                              Dfcbyn_RI[ 0].imag-Dfcbyn_RI[-1].imag,
                              1/tbytRI_RI[-1], 
                              1/tbytRI_RI[-1]])
        for iPar in range(nPars): p.add(name=ParNames[iPar],value=GuessPars[iPar])
        mini = lmfit.Minimizer(residuals_lmfit,p,nan_policy = 'omit')
        return mini
    FitVals = np.ones(nPars)*np.nan
    StdErrs = np.ones(nPars)*np.nan
    mini = configure_Minimizer()
    result = mini.minimize(method='least_squares')
    for iPar in range(nPars):
        FitVals[iPar] = result.params[ParNames[iPar]].value
        StdErrs[iPar] = result.params[ParNames[iPar]].stderr
    Dfcbyn_Extrapol = FitVals[0]+1j*FitVals[1]
    StdErr_Extrapol = StdErrs[0]+1j*StdErrs[1]
    amplitude  = FitVals[2]+1j*FitVals[3]
    om_complex = FitVals[4]+1j*FitVals[5]
    Dfcbyn_RI_Fit = Spiral(tbytRI_RI,FitVals)
    return Dfcbyn_Extrapol,StdErr_Extrapol,amplitude,om_complex,Dfcbyn_RI_Fit

def Calc_DriftFitResults(SPs,tbytRIFits,DfcbynFits,amplitude,om_complex,countFits):
    if countFits > 5 and np.abs(amplitude) > 1e-22 and np.abs(om_complex.imag) < 1e12 :
        LowLim = int(countFits*0.6)
        AvgtbytRI  = np.nanmean(tbytRIFits[LowLim:countFits])
        AvgDfcbyn  = np.nanmean(DfcbynFits[LowLim:countFits])
        DriftFitResults40perc = np.nanmean((DfcbynFits[LowLim:countFits] - AvgDfcbyn)*\
                                      (tbytRIFits[LowLim:countFits] - AvgtbytRI)) / \
                           np.nanmean((tbytRIFits[LowLim:countFits] - AvgtbytRI)**2)
        LowLim = int(countFits*0.8)
        AvgtbytRI  = np.nanmean(tbytRIFits[LowLim:countFits])
        AvgDfcbyn  = np.nanmean(DfcbynFits[LowLim:])
        DriftFitResults20perc = np.nanmean((DfcbynFits[LowLim:countFits] - AvgDfcbyn)*\
                                      (tbytRIFits[LowLim:countFits] - AvgtbytRI)) / \
                           np.nanmean((tbytRIFits[LowLim:countFits] - AvgtbytRI)**2)
    else : DriftFitResults40perc = np.nan; DriftFitResults20perc = np.nan;
    return DriftFitResults40perc,DriftFitResults20perc