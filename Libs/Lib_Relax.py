import numpy as np
from numba import jit

@jit(nopython = True)
def Relax_3D(hl,cxs,cys,czs,wi,ibars,
             tauInvl,tauInvl_Asym,one_m_tauInvl_m_Iom,one_m_tauInvl_m_Iom_Asym):
    ilist = [0,1,3,5,7,9,11,13,15,17]
    drl = np.sum(hl    )
    uxl = np.sum(hl*cxs)
    uyl = np.sum(hl*cys)
    uzl = np.sum(hl*czs)
    cdotu = cxs*uxl + cys*uyl + czs*uzl
    heq = wi*(drl + 3.*cdotu)
    for ibar in ilist:
        i = ibars[ibar]
        hSymm   = (hl[ ibar]+hl[ i])/2.
        hAsym   = (hl [ibar]-hl[ i])/2.
        heqSymm = (heq[ibar]+heq[i])/2.
        heqAsym = (heq[ibar]-heq[i])/2.
        hcolSymm = one_m_tauInvl_m_Iom      * hSymm + tauInvl      * heqSymm
        hcolAsym = one_m_tauInvl_m_Iom_Asym * hAsym + tauInvl_Asym * heqAsym
        hl[ibar] = hcolSymm + hcolAsym
        hl[i   ] = hcolSymm - hcolAsym
    return hl

### simpler, never called
@jit(nopython = True)
def Relax_3D_BGK(hl,cxs,cys,czs,wi,tauInvl,one_m_tauInvl_m_Iom):
    drl = np.sum(hl    )
    uxl = np.sum(hl*cxs)
    uyl = np.sum(hl*cys)
    uzl = np.sum(hl*czs)
    cdotu = cxs*uxl + cys*uyl + czs*uzl
    heq = wi*(drl + 3.*cdotu)
    hl = one_m_tauInvl_m_Iom * hl + tauInvl * heq
    return hl

@jit(nopython = True)
def Relax_2D(hl,cxs,cys,wi,ibars,
             tauInvl,tauInvl_Asym,one_m_tauInvl_m_Iom,one_m_tauInvl_m_Iom_Asym):
    ilist = [0,1,2,5,6]
    drl = np.sum(hl    )
    uxl = np.sum(hl*cxs)
    uyl = np.sum(hl*cys)
    cdotu = cxs*uxl + cys*uyl
    heq = wi*(drl + 3.*cdotu)
    for ibar in ilist:
        i = ibars[ibar]
        hSymm   = (hl[ ibar]+hl[ i])/2.
        hAsym   = (hl [ibar]-hl[ i])/2.
        heqSymm = (heq[ibar]+heq[i])/2.
        heqAsym = (heq[ibar]-heq[i])/2.
        hcolSymm = one_m_tauInvl_m_Iom      * hSymm + tauInvl      * heqSymm
        hcolAsym = one_m_tauInvl_m_Iom_Asym * hAsym + tauInvl_Asym * heqAsym
        hl[ibar] = hcolSymm + hcolAsym
        hl[i   ] = hcolSymm - hcolAsym
    return hl

### simpler, never called
@jit(nopython = True)
def Relax_2D_BGK(hl,cxs,cys,wi,tauInvl,one_m_tauInvl_m_Iom):
    drl = np.sum(hl    )
    uxl = np.sum(hl*cxs)
    uyl = np.sum(hl*cys)
    cdotu = cxs*uxl + cys*uyl
    heq = wi*(drl + 3.*cdotu)
    hl = one_m_tauInvl_m_Iom * hl + tauInvl * heq
    return hl

