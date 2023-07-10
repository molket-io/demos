## A module includes the fundamental constants, units, and conversion factors 
#  needed in atomic and molecular physics/theoretical and computational chemistry.
#

module constants

    const cm1 = 0.000004556335252750422::Float64 # used to convert the potential from Eh to cm^-1
                                                 # reference: NIST ....
    const Ang = 1.889726133921251660652::Float64 # used to convert from Ang to bohrs 
                                                 # reference: 
    const mC = 21874.6617573460462154::Float64   # mass of Oxygen atom in atomic units
                                                 # mass                in amu 
    const mO = 29156.94559967216628138::Float64
    #const pi=3.1415926535897932384626433832795029::Float64
    const NA     = 6.02214199e23::Float64          # 1998 CODATA Avogadro value from NIST
    const Ang    = 1/0.5291772083::Float64         # 1998 CODATA
    const eV     = 1/27.2113834::Float64           # 1998 CODATA (was 1/27.2113961)
    const cm1    = 1/219474.63137098::Float64      # 1998 CODATA (was = 1/219474.63068)
    const amu    = 1/5.485799110e-4::Float64       # 1998 CODATA
    # amu = 1/5.4857990946e-4      // 2010 CODATA (electron mass in u)
    kcal   = 1/627.50956::Float64
    c0     = 299792458.::Float64              # m/s, by definition
    coulomb= 1/1.602176462e-19::Float64      # 1998 CODATA, elementary charge e=1.6e-19 C
    second = 1/2.418884326500e-17::Float64   # 1/(E_h/hbar) from NIST "non SI units"
    c0au   = c0*1e+10*Ang/second::Float64    # speed of light in au:
    joule  = 1/4.35974381e-18::Float64       # E_h, NIST
    kboltz = 1.3806503e-23*joule::Float64    # NIST J/K * E_h/J = E_h/K
    kHz    = cm1*1e1/c0::Float64
    MHz    = cm1*1e4/c0::Float64
    GHz    = cm1*1e7/c0::Float64
    THz    = cm1*1e10/c0::Float64
    Vcm1   = joule/(1e8*Ang*coulomb)::Float64        # 1 V/cm = 3.112D-09 E_h/(e a_0); (V = J/C)
    tesla  = joule*second/(coulomb*1e20*Ang^2)::Float64      # magnetic induction in SI
                                                    # from F = q*v*B
    gauss  = 1e-4*tesla::Float64                             # Gaussian units
    g_e    = 2.0023193043718::Float64                        # Lande g factor for electron
                                                    # (gyromagnetic ratio)
    g_p    = 5.585694713::Float64                            # NIST 2010 CODATA
    gamma_p= 2.675222005e8/(second*tesla)::Float64
    # frequency to wavenumber: 1MHz = 1e4/c0 cm-1
    #  masses of most abundant isotopes from NIST 2003
    #  http://physics.nist.gov/PhysRefData/Compositions/index.html
    mH    = 1.0078250321*amu::Float64        # 99.9885 % natural abundance
    mD    = 2.0141017780*amu::Float64        #  0.0115 % natural abundance
    mHe   = 4.0026032497*amu::Float64
    mHe3  = 3.0160293097*amu::Float64
    mLi   = 7.0160040*amu::Float64
    mLi6  = 6.0151223*amu::Float64
    mBe   = 9.0121821*amu::Float64
    mB    = 11.0093055*amu::Float64
    mB10  = 10.0129370*amu::Float64
    mC    = 12*amu::Float64
    mC13  = 13.0033548378*amu::Float64
    mN    = 14.0030740052*amu::Float64       # 99.632 % natural abundance
    mN15  = 15.0001088984*amu::Float64       #  0.368 % natural abundance
    mO    = 15.9949146221*amu::Float64
    mO17  = 16.99913150*amu::Float64
    mO18  = 17.9991604*amu::Float64
    mF    = 18.99840320*amu::Float64
    mNe   = 19.9924401759*amu::Float64
    mNe21 = 20.99384674*amu::Float64
    mNe22 = 21.99138551*amu::Float64
    mNa   = 22.98976967*amu::Float64
    mMg   = 23.98504190*amu::Float64
    mAl   = 26.98153844*amu::Float64
    mSi   = 27.9769265327*amu::Float64
    mP    = 30.97376151*amu::Float64
    mS    = 31.97207069*amu::Float64
    mCl   = 34.96885271*amu::Float64 # 75.78 abundance
    mCl37 = 36.96590260*amu::Float64 # 24.22 abundance
    mAr   = 39.962383123*amu::Float64
    mK    = 38.9637069*amu::Float64
    mCa   = 39.9625912*amu::Float64
    mI    = 126.904468*amu::Float64
    mBr79 = 78.9183376*amu::Float64  # abundance 50.69
    mBr   = mBr79::Float64
    mBr81 = 80.916291*amu::Float64   # abundance 49.31
    
    mXe   = 131.293*amu::Float64     # standard atomic weight
    mSc45   = 44.9559102*amu::Float64        # NIST
    
    
    # astrostuff
    deg         = 2.0*pi/360.0::Float64
    arcminute   = deg/60.0::Float64
    arcsecond   = arcminute/60.0::Float64
    mas         = 1e-3*arcsecond::Float64    # mas = milli-arc-second
    muas        = 1e-6*arcsecond::Float64    # muas = micro-arc-second
    
    #
    mP          = 1.007276466812*amu::Float64 # NIST 2010 CODATA, proton mass in u
    muN         = 1/(2*mP)::Float64           # nuclear magneton
    
    
end # end of module
