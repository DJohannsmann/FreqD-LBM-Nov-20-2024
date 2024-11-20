import numpy as np
from numba import jit

@jit(nopython=True)
def DisLTDMin3D(nx,ny,nz,x1,x2,y1,y2,z1,z2,DMin):
    r1 = ((x1-x2   )**2+(y1-y2)**2+(z1-z2   )**2)**0.5
    r2 = ((x1-x2-nx)**2+(y1-y2)**2+(z1-z2   )**2)**0.5
    r3 = ((x1-x2+nx)**2+(y1-y2)**2+(z1-z2   )**2)**0.5
    r4 = ((x1-x2   )**2+(y1-y2)**2+(z1-z2-nz)**2)**0.5
    r5 = ((x1-x2   )**2+(y1-y2)**2+(z1-z2+nz)**2)**0.5
    r6 = ((x1-x2-nx)**2+(y1-y2)**2+(z1-z2-nz)**2)**0.5
    r7 = ((x1-x2-nx)**2+(y1-y2)**2+(z1-z2+nz)**2)**0.5
    r8 = ((x1-x2+nx)**2+(y1-y2)**2+(z1-z2-nz)**2)**0.5
    r9 = ((x1-x2+nx)**2+(y1-y2)**2+(z1-z2+nz)**2)**0.5
    if (r1<=DMin)or(r2<=DMin)or(r3<=DMin)or(r4<=DMin)or(r5<=DMin)\
     or(r6<=DMin)or(r7<=DMin)or(r8<=DMin)or(r9<=DMin):
        lessthanDMin = True
    else: lessthanDMin = False    
    return lessthanDMin
    
@jit(nopython=True)
def DisLTDMin_inPlane(nx,nz,x1,x2,z1,z2,DMin):
    r1 = ((x1-x2   )**2+(z1-z2   )**2)**0.5
    r2 = ((x1-x2-nx)**2+(z1-z2   )**2)**0.5
    r3 = ((x1-x2+nx)**2+(z1-z2   )**2)**0.5
    r4 = ((x1-x2   )**2+(z1-z2-nz)**2)**0.5
    r5 = ((x1-x2   )**2+(z1-z2+nz)**2)**0.5
    r6 = ((x1-x2-nx)**2+(z1-z2-nz)**2)**0.5
    r7 = ((x1-x2-nx)**2+(z1-z2+nz)**2)**0.5
    r8 = ((x1-x2+nx)**2+(z1-z2-nz)**2)**0.5
    r9 = ((x1-x2+nx)**2+(z1-z2+nz)**2)**0.5
    if (r1<=DMin)or(r2<=DMin)or(r3<=DMin)or(r4<=DMin)or(r5<=DMin)\
     or(r6<=DMin)or(r7<=DMin)or(r8<=DMin)or(r9<=DMin):
        lessthanDMin = True
    else: lessthanDMin = False    
    return lessthanDMin

@jit(nopython=True)
def Check_for_Overlap(nx,nz,xSphs,zSphs,RSph,nSph,Gap_P2P):
    Overlap = False
    for i in range(nSph):
        for j in range(nSph):
            if i != j:
                if DisLTDMin_inPlane(nx,nz,xSphs[i],xSphs[j],zSphs[i],zSphs[j],2*RSph+Gap_P2P):
                  Overlap = True 
    return Overlap            

def Set_SphPoss_Random(nx,ny,nz,nSph,RSph,ySphbyR,Gap_P2P):
    xSphs = np.ones(nSph,dtype=np.float64)*np.nan
    ySphs = np.ones(nSph,dtype=np.float64)*np.nan
    zSphs = np.ones(nSph,dtype=np.float64)*np.nan
    ic_max = 1000000
    Overlap = True; ic = 0 
    if nSph>1:     
        while Overlap and ic<ic_max:
            for iS in range(nSph):
                xSphs[iS] = np.random.rand()*(nx-1)
                ySphs[iS] = RSph*ySphbyR
                zSphs[iS] = np.random.rand()*(nz-1)
            Overlap = Check_for_Overlap(nx,nz,xSphs,zSphs,RSph,nSph,Gap_P2P)
            ic += 1
            if ic == ic_max-1: 
                print('Set Shere Pos failed, coverage too high?')
                SphPoss = np.ones((3,nSph))*np.nan

        # print('# attempts positioning',ic)
    else: 
        xSphs[0] = (nx-1)/2.
        ySphs[0] = RSph*ySphbyR
        zSphs[0] = (nz-1)/2.
    SphPoss = np.array([xSphs,ySphs,zSphs])
    return SphPoss

def Set_SphPoss(SPs):
    nx = SPs['nx']
    ny = SPs['ny']
    nz = SPs['nz']
    nSph = SPs['nSph']
    RSph = SPs['RSph']
    ySphbyR = SPs['ySphbyR']
    Gap_P2P = SPs['Gap_P2P']
    SphPoss = Set_SphPoss_Random(nx,ny,nz,nSph,RSph,ySphbyR,Gap_P2P)
    return SphPoss


