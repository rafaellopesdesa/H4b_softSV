ELECTRON = 11
POSITRON = -ELECTRON
EMINUS = ELECTRON
EPLUS = POSITRON
MUON = 13
ANTIMUON = -MUON
TAU = 15
ANTITAU = -TAU

NU_E = 12
NU_EBAR = -NU_E
NU_MU = 14
NU_MUBAR = -NU_MU
NU_TAU = 16
NU_TAUBAR = -NU_TAU

PHOTON = 22
GAMMA = PHOTON
GLUON = 21
WPLUSBOSON = 24
WMINUSBOSON = -WPLUSBOSON
WPLUS = WPLUSBOSON
WMINUS = WMINUSBOSON
Z0BOSON = 23
ZBOSON = Z0BOSON
Z0 = Z0BOSON
HIGGSBOSON = 25
HIGGS = HIGGSBOSON

DQUARK = 1
UQUARK = 2
SQUARK = 3
CQUARK = 4
BQUARK = 5
TQUARK = 6

PROTON = 2212
ANTIPROTON = -PROTON
PBAR = ANTIPROTON
NEUTRON = 2112
ANTINEUTRON = -NEUTRON

PI0 = 111
PIPLUS = 211
PIMINUS = -PIPLUS
K0L = 130
K0S = 310
KPLUS = 321
KMINUS = -KPLUS
ETA = 221
ETAPRIME = 331
PHI = 333
OMEGA = 223

ETAC = 441
JPSI = 443
PSI2S = 100443

D0 = 421
DPLUS = 411
DMINUS = -DPLUS
DSPLUS = 431
DSMINUS = -DSPLUS
ETAB = 551
UPSILON1S = 553
UPSILON2S = 100553
UPSILON3S = 200553
UPSILON4S = 300553

B0 = 511
BPLUS = 521
BMINUS = -BPLUS
B0S = 531
BCPLUS = 541
BCMINUS = -BCPLUS

LAMBDA = 3122
SIGMA0 = 3212
SIGMAPLUS = 3222
SIGMAMINUS = 3112
LAMBDACPLUS = 4122
LAMBDACMINUS = 4122
LAMBDAB = 5122
XI0 = 3322
XIMINUS = 3312
XIPLUS = -XIMINUS
OMEGAMINUS = 3334
OMEGAPLUS = -OMEGAMINUS

REGGEON = 110
POMERON = 990
ODDERON = 9990
GRAVITON = 39
NEUTRALINO1 = 1000022
GRAVITINO = 1000039
GLUINO = 1000021

Location = {'nj':1, 'nq3':2, 'nq2':3, 'nq1':4, 'nl':5, 'nr':6, 'n':7, 'n8':8, 'n9':9, 'n10':10}

def abspid(pid):
    return abs(pid)

def extraBits(pid):
    return abspid(pid)//10000000

def digit(loc, pid):
    numerator = 10**(Location[loc]-1)
    return (abspid(pid)//numerator) % 10

def fundamentalID(pid):
    
        if (extraBits(pid) > 0):
            return 0;
        if (digit('nq2',pid) == 0 and digit('nq1',pid) == 0):
            return abspid(pid)%10000
        elif (abspid(pid) <= 100): 
            return abspid(pid)
        else:
            return 0

def isNucleus(pid):
    if (abspid(pid) == 2212):
        return True

    if ((digit('n10',pid) == 1 ) and ( digit('n9',pid) == 0 )) :
        if( (abspid(pid)//10)%1000 >= (abspid(pid)//10000)%1000 ):
            return True
        return False

def Z(pid):

    if (abspid(pid) == 2212):
        return 1
    if( isNucleus(pid) ):
        return (abspid(pid)//10000)%1000;
    return 0;


def A(pid):

    if( abspid(pid) == 2212 ):
        return 1
    if( isNucleus(pid) ):
        return (abspid(pid)//10)%1000;
    return 0

def Lambda(pid):

    if( abspid(pid) == 2212 ):
        return 0
    if( isNucleus(pid) ):
        return digit('n8',pid)
    return 0

def isMeson(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( abspid(pid) <= 100 ):
        return False
    if( fundamentalID(pid) <= 100 and fundamentalID(pid) > 0 ):
        return False
    aid = abspid(pid);
    if( aid == 130 or aid == 310 or aid == 210 ):
        return True
    if( aid == 150 or aid == 350 or aid == 510 or aid == 530 ):
        return True
    if( pid == 110 or pid == 990 or pid == 9990 ):
        return True
    if(    digit('nj',pid) > 0 and digit('nq3',pid) > 0 and digit('nq2',pid) > 0 and digit('nq1',pid) == 0 ):
        if( digit('nq3',pid) == digit('nq2',pid) and pid < 0 ):
            return False
        else:
            return True
    return False

def isBaryon(pid):
    
    if( extraBits(pid) > 0 ):
        return False
    if( abspid(pid) <= 100 ):
        return False
    if( fundamentalID(pid) <= 100 and fundamentalID(pid) > 0 ):
        return False
    if( abspid(pid) == 2110 or abspid(pid) == 2210 ):
        return True
    if(    digit('nj',pid) > 0 and digit('nq3',pid) > 0 and digit('nq2',pid) > 0 and digit('nq1',pid) > 0 ):
        return True
    return False

def isDiQuark(pid):
    
    if( extraBits(pid) > 0 ):
        return False    
    if( abspid(pid) <= 100 ):
        return False
    if( fundamentalID(pid) <= 100 and fundamentalID(pid) > 0 ):
        return False
    if(    digit('nj',pid) > 0 and digit('nq3',pid) == 0 and digit('nq2',pid) > 0 and digit('nq1',pid) > 0 ):
        return True
    return False


def isPentaquark(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( digit('n',pid) != 9 ):
        return False
    if( digit('nr',pid) == 9 or digit('nr',pid) == 0 ):
        return False
    if( digit('nj',pid) == 9 or digit('nl',pid) == 0 ):
        return False
    if( digit('nq1',pid) == 0 ):
        return False
    if( digit('nq2',pid) == 0 ):
        return False
    if( digit('nq3',pid) == 0 ):
        return False
    if( digit('nj',pid) == 0 ):
        return False
    if( digit('nq2',pid) > digit('nq1',pid) ):
        return False
    if( digit('nq1',pid) > digit('nl',pid) ):
        return False
    if( digit('nl',pid) > digit('nr',pid) ):
        return False
    return True

def isHadron(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( isMeson(pid) ):
        return True
    if( isBaryon(pid) ):
        return True
    if( isPentaquark(pid) ):
        return True
    return False

def isLepton(pid):

    if( extraBits(pid) > 0 ):
        return False
    if( fundamentalID(pid) >= 11 and fundamentalID(pid) <= 18 ):
        return True
    return False

def isSUSY(pid):

    if( extraBits(pid) > 0 ):
        return False
    if( digit('n',pid) != 1 and digit('n',pid) != 2 ):
        return False
    if( digit('nr',pid) != 0 ):
        return False
    if( fundamentalID(pid) == 0 ):
        return False
    return True

def isRhadron(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( digit('n',pid) != 1 ):
        return False
    if( digit('nr',pid) != 0 ):
        return False
    if( isSUSY(pid) ):
        return False
    if( digit('nq2',pid) == 0 ):
        return False
    if( digit('nq3',pid) == 0 ):
        return False
    if( digit('nj',pid) == 0 ):
        return False
    return True

def hasUp(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( fundamentalID(pid) > 0 ):
        return False
    if( digit('nq3',pid) == 2 or digit('nq2',pid) == 2 or digit('nq1',pid) == 2 ):
        return True
    return False

def hasDown(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( fundamentalID(pid) > 0 ):
        return False
    if( digit('nq3',pid) == 1 or digit('nq2',pid) == 1 or digit('nq1',pid) == 1 ):
        return True
    return False

def hasStrange(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( fundamentalID(pid) > 0 ):
        return False
    if( digit('nq3',pid) == 3 or digit('nq2',pid) == 3 or digit('nq1',pid) == 3 ):
        return True
    return False

def hasCharm(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( fundamentalID(pid) > 0 ):
        return False
    if( digit('nq3',pid) == 4 or digit('nq2',pid) == 4 or digit('nq1',pid) == 4 ):
        return True
    return False

def hasBottom(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( fundamentalID(pid) > 0 ):
        return False
    if( digit('nq3',pid) == 5 or digit('nq2',pid) == 5 or digit('nq1',pid) == 5 ):
        return True
    return False

def hasTop(pid):
    if( extraBits(pid) > 0 ):
        return False
    if( fundamentalID(pid) > 0 ):
        return False
    if( digit('nq3',pid) == 6 or digit('nq2',pid) == 6 or digit('nq1',pid) == 6 ):
        return True
    return False


def isValid(pid):
    if( extraBits(pid) > 0 ):
        if( isNucleus(pid) ):
            return True
        return False
    if( isSUSY(pid) ):
        return True
    if( isRhadron(pid) ):
        return True
    if( isMeson(pid) ):
        return True
    if( isBaryon(pid) ):
        return True
    if( isDiQuark(pid) ):
        return True
    if( fundamentalID(pid) > 0 ):
        if(pid > 0 ):
            return True
        else:
            return False
    if( isPentaquark(pid) ):
        return True
    return False

def jSpin(pid):
    if( fundamentalID(pid) > 0 ):
        fund = fundamentalID(pid)
        if( fund > 0 and fund < 7 ):
            return 2
        if( fund == 9 ):
            return 3
        if( fund > 10 and fund < 17 ):
            return 2
        if( fund > 20 and fund < 25 ):
            return 3
        return 0
    elif( extraBits(pid) > 0 ):
        return 0
    return abspid(pid)%10

def sSpin(pid):
    if( not isMeson(pid) ):
        return 0
    if( digit('n',pid) == 9 ):
        return 0
    inl = digit('nl',pid)
    js = digit('nj',pid)

    if( inl == 0 and js >= 3 ):
        return 1
    elif ( inl == 0  and js == 1 ):
        return 0
    elif ( inl == 1  and js >= 3 ):
        return 0
    elif ( inl == 2  and js >= 3 ):
        return 1
    elif ( inl == 1  and js == 1 ):
        return 1
    elif ( inl == 3  and js >= 3 ):
        return 1
    return 0

def lSpin(pid):
    if( not isMeson(pid) ):
        return 0
    if( digit('n',pid) == 9 ):
        return 0
    inl = digit('nl',pid)
    js = digit('nj',pid)
    if( inl == 0 and js == 3 ):
        return 0;
    elif ( inl == 0 and js == 5 ):
        return 1
    elif ( inl == 0 and js == 7 ):
        return 2
    elif ( inl == 0 and js == 9 ):
        return 3
    elif ( inl == 0  and js == 1 ):
        return 0
    elif ( inl == 1  and js == 3 ):
        return 1
    elif ( inl == 1  and js == 5 ):
        return 2
    elif ( inl == 1  and js == 7 ):
        return 3
    elif ( inl == 1  and js == 9 ):
        return 4
    elif ( inl == 2  and js == 3 ):
        return 1
    elif ( inl == 2  and js == 5 ):
        return 2
    elif ( inl == 2  and js == 7 ):
        return 3
    elif ( inl == 2  and js == 9 ):
        return 4
    elif ( inl == 1  and js == 1 ):
        return 1
    elif ( inl == 3  and js == 3 ):
        return 2
    elif ( inl == 3  and js == 5 ):
        return 3
    elif ( inl == 3  and js == 7 ):
        return 4
    elif ( inl == 3  and js == 9 ):
        return 5
    return 0

def threeCharge(pid):
    charge=0
    ch100 = [ -1, 2,-1, 2,-1, 2,-1, 2, 0, 0,
              -3, 0,-3, 0,-3, 0,-3, 0, 0, 0,
              0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 3, 0, 0, 3, 0, 0, 0,
              0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 6, 3, 6, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    q1 = digit('nq1',pid);
    q2 = digit('nq2',pid)
    q3 = digit('nq3',pid)
    ida = abspid(pid)
    sid = fundamentalID(pid)
    if( ida == 0 or extraBits(pid) > 0 ):
        return 0
    elif ( sid > 0 and sid <= 100 ):
        charge = ch100[sid-1]
        if(ida==1000017 or ida==1000018):
            charge = 0
        if(ida==1000034 or ida==1000052):
            charge = 0
        if(ida==1000053 or ida==1000054):
            charge = 0
        if(ida==5100061 or ida==5100062):
            charge = 6
    elif( digit('nj',pid) == 0 ):
        return 0
    elif( isMeson(pid) ):
        if( q2 == 3 or q2 == 5 ):
            charge = ch100[q3-1] - ch100[q2-1]
        else:
            charge = ch100[q2-1] - ch100[q3-1];
                
    elif( isDiQuark(pid) ):
        charge = ch100[q2-1] + ch100[q1-1]
    elif( isBaryon(pid) ):
        charge = ch100[q3-1] + ch100[q2-1] + ch100[q1-1]
    else:
        return 0

    if( charge == 0 ):
        return 0
    elif( pid < 0 ):
        charge = -charge
        
    return charge
    


def charge(pid):
    return threeCharge(pid)/3.0

def isCharged(pid):
    return threeCharge(pid) != 0

def isNeutral(pid):
    return threeCharge(pid) == 0

def isQuark(pid):
    return (abs(pid) >= 1 and abs(pid)<7)

def isParton(pid):
    return (pid == GLUON or isQuark(pid))

def isPhoton(pid):
    return (pid == PHOTON)

def isElectron(pid):
    return (abs(pid) == ELECTRON)

def isMuon(pid):
    return (abs(pid) == MUON)

def isTau(pid):
    return (abs(pid) == TAU)

def isChLepton(pid):
    apid = abs(pid)
    return (apid == 11 or apid == 13 or apid == 15)

def isNeutrino(pid):
    apid = abs(pid)
    return (apid == 12 or apid == 14 or apid == 16)

def isWplus(pid):
    return (pid == WPLUSBOSON)

def isWminus(pid):
    return (pid == WMINUSBOSON)

def isW(pid):
    return (abs(pid) == WPLUSBOSON)
        
def isZ(pid):
    return (pid == Z0BOSON)

def isHiggs(pid):
    return (pid == HIGGSBOSON or pid == 26)
  
def isBSMHiggs(pid):
    return (pid == 35 or pid == 36)    

def isTop(pid):
    return (abs(pid) == 6)

def isLightMeson(pid):
    return (isMeson(pid) and not (hasCharm(pid) or hasBottom(pid)))

def isLightBaryon(pid):
    return (isBaryon(pid) and not (hasCharm(pid) or hasBottom(pid)))

def isLightHadron(pid):
    return (isHadron(pid) and not (hasCharm(pid) or hasBottom(pid)))

def isHeavyMeson(pid):
    return (isMeson(pid) and (hasCharm(pid) or hasBottom(pid)))

def isHeavyBaryon(pid):
    return (isBaryon(pid) and (hasCharm(pid) or hasBottom(pid)))

def isHeavyHadron(pid):
    return (isHadron(pid) and (hasCharm(pid) or hasBottom(pid)))

def isBottomMeson(pid):
    return (hasBottom(pid) and isMeson(pid))

def isBottomBaryon(pid):
    return (hasBottom(pid) and isBaryon(pid))

def isBottomHadron(pid):
    return (hasBottom(pid) and isHadron(pid))

def isCharmMeson(pid):
    return (isMeson(pid) and hasCharm(pid))

def isCharmBaryon(pid):
    return (isBaryon(pid) and hasCharm(pid))

def isCharmHadron(pid):
    return (isHadron(pid) and hasCharm(pid))

def isStrangeMeson(pid):
    return (isMeson(pid) and hasStrange(pid))

def isStrangeBaryon(pid):
    return (isBaryon(pid) and hasStrange(pid))

def isStrangeHadron(pid):
    return (isHadron(pid) and hasStrange(pid))

def isGenSpecific(pid):
    return (pid >= 80 and pid < 101)

def isResonance(pid):
    return (isW(pid) or isZ(pid) or isHiggs(pid) or isTop(pid))

def isTransportable(pid):
    return (isPhoton(pid) or isHadron(pid) or isLepton(pid))
