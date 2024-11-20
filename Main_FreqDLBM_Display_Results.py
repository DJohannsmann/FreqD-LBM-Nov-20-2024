import numpy as np
import matplotlib.pyplot as plt
import configparser

Colors = ['r','g','b','fuchsia','navy','orange','olive','cyan','darksalmon']
fname = r'flbm\FreqD-LBM-Output_2024-11-14-03_49_55.txt'
OvtLabels = ['15 MHz','25','35','45']
data = np.loadtxt(fname,skiprows=1); print('data.shape',data.shape)
config = configparser.ConfigParser(); config.read(fname[:-4]+'.cfg')
config.optionxform = str

SPs = {}
for section in config.sections():
    for option in config.options(section): SPs[option]=config.get(section,option)

print(SPs['par1str'],SPs['par2str'],SPs['par3str'])
navg  = int(SPs['navg' ])
nPar1 = int(SPs['npar1'])
nPar2 = int(SPs['npar2'])
nPar3 = int(SPs['npar3'])
novt  = int(SPs['novt'])


Dfbyns_all = np.ones((navg,nPar1,nPar2,nPar3,novt))*np.nan
DGbyns_all = np.ones((navg,nPar1,nPar2,nPar3,novt))*np.nan

Pars1 = np.ones((nPar1))*np.nan
Pars2 = np.ones((nPar2))*np.nan
Pars3 = np.ones((nPar3))*np.nan
ns    = np.ones((novt))*np.nan

if SPs['par1str'] == 'CovTarget' : SPs['par1str'] = 'coverage'
if SPs['par2str'] == 'CovTarget' : SPs['par2str'] = 'coverage'
if SPs['par3str'] == 'CovTarget' : SPs['par3str'] = 'coverage'

for il in range(len(data)):
    iavg  = int(data[il,0])
    iPar1 = int(data[il,1])
    iPar2 = int(data[il,2])
    iPar3 = int(data[il,3])
    iovt  = int(data[il,4])
    if SPs['par1str'] == 'coverage' : Pars1[iPar1]  = data[il,16]
    else : Pars1[iPar1]  = data[il,5] 
    if SPs['par2str'] == 'coverage' : Pars2[iPar2]  = data[il,16]
    else : Pars2[iPar2]  = data[il,6] 
    if SPs['par3str'] == 'coverage' : Pars3[iPar3]  = data[il,16]
    else : Pars3[iPar2]  = data[il,7] 
    ns[iovt]      = data[il,8]
    Dfbyns_all[iavg,iPar1,iPar2,iPar3,iovt]  = data[il,9] 
    DGbyns_all[iavg,iPar1,iPar2,iPar3,iovt]  = data[il,10]

Dfbyns_avg = np.nanmean(Dfbyns_all,axis=0)    
DGbyns_avg = np.nanmean(DGbyns_all,axis=0)    
Dfbyns_std = np.nanstd( Dfbyns_all,axis=0)    
DGbyns_std = np.nanstd( DGbyns_all,axis=0)    


plt.rcParams["figure.figsize"] = (5,3.);
plt.rcParams["font.size"] = 11;
for iPar1,Par1 in enumerate(Pars1):
    for iPar2,Par2 in enumerate(Pars2):
        plt.subplot(2,2,1)
        for iovt in range(novt): 
            plt.errorbar(Pars3,Dfbyns_avg[iPar1,iPar2,:,iovt],Dfbyns_std[iPar1,iPar2,:,iovt],\
            label = OvtLabels[iovt],color = Colors[iovt])
        plt.legend(fontsize = 7,loc = 'upper right')
        plt.ylabel('$\Delta f/n$ [Hz]'); 
        plt.xticks([])
    
        plt.subplot(2,2,3)
        for iovt in range(novt): 
            plt.errorbar(Pars3,DGbyns_avg[iPar1,iPar2,:,iovt],DGbyns_std[iPar1,iPar2,:,iovt],color = Colors[iovt])
        plt.ylabel('$\Delta \Gamma /n$ [Hz]'); 
        plt.xlabel(SPs['par3str'])
    
        plt.subplot(2,2,2)
        for iovt in range(novt): 
            Df_ratio_stds = ((Dfbyns_std[iPar1,iPar2,:,iovt]/Dfbyns_avg[iPar1,iPar2,:,iovt])**2 + (DGbyns_std[iPar1,iPar2,:,iovt]/DGbyns_avg[iPar1,iPar2,:,iovt])**2)**0.5*np.abs(DGbyns_avg[iPar1,iPar2,:,iovt]/Dfbyns_avg[iPar1,iPar2,:,iovt])
            plt.errorbar(-Dfbyns_avg[iPar1,iPar2,:,iovt],\
                -DGbyns_avg[iPar1,iPar2,:,iovt]/Dfbyns_avg[iPar1,iPar2,:,iovt],\
                  Df_ratio_stds,color = Colors[iovt])
        plt.ylabel('$\Delta \Gamma /(-\Delta f)$'); 
        plt.xticks([])
    
        plt.subplot(2,2,4)
        for iovt in range(novt): 
            ff_ratio_stds = ((Dfbyns_std[iPar1,iPar2,:,iovt]/Dfbyns_avg[iPar1,iPar2,:,iovt])**2 + \
                             (Dfbyns_std[iPar1,iPar2,:,0]/Dfbyns_avg[iPar1,iPar2,:,0])**2)**0.5*\
                              Dfbyns_avg[iPar1,iPar2,:,iovt]/Dfbyns_avg[iPar1,iPar2,:,0]
            plt.errorbar(-Dfbyns_avg[iPar1,iPar2,:,iovt],\
                1-Dfbyns_avg[iPar1,iPar2,:,iovt]/Dfbyns_avg[iPar1,iPar2,:,0],\
                ff_ratio_stds,color = Colors[iovt])

        plt.xlabel('$-\Delta f/n$ [Hz]'); 
        plt.ylabel('1 $-$ ( $(\Delta f/n)_{n}$ /\n $\Delta f/n_{low}$ )'); 
        

        plt.suptitle(SPs['par1str'] + ': ' + str(Par1) + '\n' +
                     SPs['par2str'] + ': ' + str(Par2),fontsize = 10)

        plt.tight_layout(); 
        plt.savefig(fname[:-4]+'.png'); 
        plt.show();   
