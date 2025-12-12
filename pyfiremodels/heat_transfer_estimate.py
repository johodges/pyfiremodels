import numpy as np

def getJHcolors():
    colors = np.array([[0, 65, 101],
                       [229, 114, 0],
                       [136, 139, 141],
                       [170, 39, 44],
                       [119, 197, 213],
                       [161, 216, 132],
                       [255, 200, 69],
                       [101, 0, 65],
                       [0, 229, 114],
                       [141, 136, 139],
                       [44, 170, 39],
                       [213, 119, 197],
                       [132, 161, 216],
                       [69, 255, 200],
                       [65, 101, 0],
                       [114, 0, 229]], dtype=np.float32)
    colors = colors/255
    return colors

def calculateCeilingJetTemperature(r, Q, H, Ti=293):
    # SFPE Handbook chapter 14
    if r > 0.18*H:
        # SFPE Handbook chapter 14 - Eq 14.3
        T = 5.38 * ( ( Q**(2/3) ) / ( H ** (5/3) ) ) / ( (r/H)**(2/3) ) + Ti
    else:
        # SFPE Handbook chapter 14 - Eq 14.2
        T = 16.9*(Q**(2/3))/(H**(5/3))+Ti
    return T

def calculateCeilingJetVelocity(r, Q, H):
    # SFPE Handbook chapter 14
    if r > 0.15*H:
        # SFPE Handbook chapter 14 - Eq 14.5
        U = 0.197 * ( (Q/H)**(1/3) ) / ( (r/H)**(5/6) )
    else:
        # SFPE Handbook chapter 14 - Eq 14.4
        U = 0.947 * ( (Q/H)**(1/3) )
    U = max([U, 0.1])
    return U

def calculateVirtualSource(Q, D):
    # Heskestad
    z0 = 1.02*D + 0.083 * (Q**0.4)
    return z0

def calculateSprinklerTemperature(r, H, t, Ti, RTI, alpha, Qmax, p=2, g=9.81):
    Q = alpha * (t ** p)
    Q = min([Q, Qmax])
    cp = getAirCpFromTempK(Ti) / 1000 # kJ/kg-C
    rho_0 = getAirRhoFromTempK(Ti)
    u = calculateCeilingJetVelocity(r, Q, H)
    Tg = calculateCeilingJetTemperature(r, Q, H, Ti)
    dT = Tg - Ti
    D = 0.146 + 0.242 * r / H
    t2sf = 0.861 * ( 1 + (r/H) )
    A = g/(cp*Ti*rho_0)
    t2s = t / (A ** (-1/(3 + p))*alpha**(-1 / (3 + p)) * H**(4 / (3 + p)))
    
    if t2s < t2sf:
        dT2s = 0
        u2s = 0
        Y = 0
        Td = Ti
    else:
        u_over_u2s = (A ** (1 / (3 + p))) * (alpha ** (1 / (3 + p)) ) * (H ** ((p-1)/(3+p)))
        dT_over_dT2s = (Ti/g)*(A ** (2 / (3 + p))) * (alpha ** (2 / (3 + p)) ) * (H ** (-(5-p)/(3+p)))
        dT2s = ( (t2s-t2sf) / (0.146 + 0.242 * (r/H)) ) ** (4/3)
        u2s_over_dt2sr = 0.59 * ( (r/H) **(-0.63) )
        Y = (3/4) * (u_over_u2s ** 0.5) * ((u2s_over_dt2sr) ** 0.5) * (dT2s / RTI) * (t / t2s) * D
        Td = dT_over_dT2s * dT2s * (1 - ((1 - np.exp(-Y)) / Y) ) + Ti
        #dT2s = ((t2s-t2sf) / (0.146 + 0.242 * (r / H)) ) ** (4/3)
        #u2s = 0.59 * ( (r / H) ** (-0.63) ) * (dT2s ** 0.5)
        #Y = (3/4) * ( (u / u2s)**0.5 ) * ( ( u2s / (dT2s ** 0.5))**0.5 ) * (dT2s / RTI) * (t/t2s) * D
        #Td = (dT / dT2s) * dT2s * (( 1 - ((1 - np.exp(-Y)) / Y))) + Ti
    return Td

def calculateSprinklerActivationTime(r, H, Ti, RTI, alpha, Qmax, Tactivation):
    t = 0.0
    Td = Ti
    while Td < Tactivation:
        t = t + 0.5
        Td_old = Td
        Td = calculateSprinklerTemperature(r, H, t, Ti, RTI, alpha, Qmax)
        if (t > 3600):
            return -1
    return t

def calculateSprinklerActivationTimeSpecifiedProfile(r, H, Ti, RTI, t1, Q1, Tactivation, dt=0.1):
    ts = np.linspace(0, np.max(t1), int(np.max(t1)/dt + 1))
    Ts = np.zeros_like(ts) + Ti
    Qs = np.interp(ts, t1, Q1)
    us = np.zeros_like(ts)
    Tgs = np.zeros_like(ts)
    for i in range(1, len(ts)):
        Q = Qs[i-1]
        u = calculateCeilingJetVelocity(r, Q, H)
        Tg = calculateCeilingJetTemperature(r, Q, H, Ti)
        dTdt = (u**0.5)*(Tg-Ts[i-1]) / RTI
        Ts[i] = Ts[i-1] + dTdt*dt
        us[i] = u
        Tgs[i] = Tg
    
    inds = np.where(Ts > Tactivation)[0]
    if len(inds) > 0:
        tActivation = ts[inds[0]]
    else:
        tActivation = np.nan
    return ts, Ts, tActivation, us, Tgs


def getGrowthRates(speed, peak=None, time=None, p=2):
    if speed == "ultrafast":
        (peak, time) = (1000, 75)
    elif speed == "fast":
        (peak, time) = (1000, 150)
    elif speed == "medium":
        (peak, time) = (1000, 300)
    elif speed == "slow":
        (peak, time) = (1000, 600)
    alpha = peak / (time ** p)
    return alpha

def calculateT(Q, z, Tamb, g=9.81):
    # Drysdale 2011 Table 4.2
    try:
        z[0]
        z = np.array(z)
    except:
        z = np.array([z])
    heightFactor = z/(Q**(2/5))
    T = np.zeros_like(z)
    for i, z0 in enumerate(z):
        if heightFactor[i] <0.08: #z <= H*0.6:
            # Continuous flaming region
            k = 6.8 # m(1/2)/s
            eta = 0.5
            C = 0.9
        elif heightFactor[i] <= 0.2: #z <= H*1.2:
            # Intermittent flaming region
            k = 1.9 # m/kW(1/5)s
            eta = 0
            C = 0.9
        else:
            # Plume above flame
            k = 1.1 # m(4/3)/kW(1/3)s
            eta = -1/3
            C = 0.9
        DT = ((k/C)**2)*(heightFactor[i]**(2*eta-1))*(Tamb/(2*g))
        T[i] = Tamb + DT
    return T

def calculateV(Q, z, Tamb, g=9.81):
    try:
        z[0]
        z = np.array(z)
    except:
        z = np.array([z])
    # Drysdale 2011 Table 4.2
    heightFactor = z/(Q**(2/5))
    u = np.zeros_like(z)
    for i, z0 in enumerate(z):
        if heightFactor[i] <0.08: #z <= H*0.6:
            # Continuous flaming region
            k = 6.8 # m(1/2)/s
            eta = 0.5
        elif heightFactor[i] <= 0.2: #z <= H*1.2:
            # Intermittent flaming region
            k = 1.9 # m/kW(1/5)s
            eta = 0
        else:
            # Plume above flame
            k = 1.1 # m(4/3)/kW(1/3)s
            eta = -1/3
        u[i] = k*Q**(1/5)*(z[i]/Q**(2/5))**eta
    return u

def getAirCpFromTempK(Temp):
    # J/kg-K
    Ts = [100, 150, 200, 250, 300,
          350, 400, 450, 500, 550,
          600, 650, 700, 750, 800,
          850, 900, 950, 1000, 1100,
          1200, 1300, 1400, 1500, 1600,
          1700, 1800, 1900, 2000, 2100,
          2200, 2300, 2400, 2500, 3000]
    cps = [1.032, 1.012, 1.007, 1.006, 1.007,
           1.009, 1.014, 1.021, 1.030, 1.040,
           1.051, 1.063, 1.075, 1.087, 1.099,
           1.110, 1.121, 1.131, 1.141, 1.159,
           1.175, 1.189, 1.207, 1.230, 1.248,
           1.267, 1.286, 1.307, 1.337, 1.372,
           1.417, 1.478, 1.558, 1.665, 2.726]
    cp = np.interp(Temp, Ts, cps)*(10**3)
    return cp

def getAirRhoFromTempK(Temp):
    # kg/m3
    Ts = [100, 150, 200, 250, 300,
          350, 400, 450, 500, 550,
          600, 650, 700, 750, 800,
          850, 900, 950, 1000, 1100,
          1200, 1300, 1400, 1500, 1600,
          1700, 1800, 1900, 2000, 2100,
          2200, 2300, 2400, 2500, 3000]
    rhos = [3.5562, 2.3364, 1.7458, 1.3947, 1.1614,
            0.9950, 0.8711, 0.7740, 0.6964, 0.6329,
            0.5804, 0.5356, 0.4975, 0.4643, 0.4354,
            0.4097, 0.3868, 0.3666, 0.3482, 0.3166,
            0.2902, 0.2679, 0.2488, 0.2322, 0.2177,
            0.2049, 0.1935, 0.1833, 0.1741, 0.1658,
            0.1582, 0.1513, 0.1448, 0.1389, 0.1135]
    rho = np.interp(Temp, Ts, rhos)
    return rho

def getAirNuFromTempK(Temp):
    # m2/s
    Ts = [100, 150, 200, 250, 300,
          350, 400, 450, 500, 550,
          600, 650, 700, 750, 800,
          850, 900, 950, 1000, 1100,
          1200, 1300, 1400, 1500, 1600,
          1700, 1800, 1900, 2000, 2100,
          2200, 2300, 2400, 2500, 3000]
    nus = [2.000,4.426,7.590,11.440,15.890,
           20.920,26.410,32.390,38.790,45.570,
           52.690,60.210,68.100,76.370,84.930,
           93.800,102.900,112.200,121.900,141.800,
           162.900,185.100,213.000,240.000,268.000,
           298.000,329.000,362.000,396.000,431.000,
           468.000,506.000,547.000,589.000,841.000]
    nu = np.interp(Temp, Ts, nus)*(10**-6)
    return nu

def getAirPrFromTempK(Temp):
    # m2/s
    Ts = [100, 150, 200, 250, 300,
          350, 400, 450, 500, 550,
          600, 650, 700, 750, 800,
          850, 900, 950, 1000, 1100,
          1200, 1300, 1400, 1500, 1600,
          1700, 1800, 1900, 2000, 2100,
          2200, 2300, 2400, 2500, 3000]
    Prs = [0.786,0.758,0.737,0.720,0.707,
           0.700,0.690,0.686,0.684,0.683,
           0.685,0.690,0.695,0.702,0.709,
           0.716,0.720,0.723,0.726,0.728,
           0.728,0.719,0.703,0.685,0.688,
           0.685,0.683,0.677,0.672,0.667,
           0.655,0.647,0.630,0.613,0.536]
    Pr = np.interp(Temp, Ts, Prs)
    return Pr

def getAirKFromTempK(Temp):
    # m2/s
    Ts = [100, 150, 200, 250, 300,
          350, 400, 450, 500, 550,
          600, 650, 700, 750, 800,
          850, 900, 950, 1000, 1100,
          1200, 1300, 1400, 1500, 1600,
          1700, 1800, 1900, 2000, 2100,
          2200, 2300, 2400, 2500, 3000]
    Ks = [9.34,13.80,18.10,22.30,26.30,
          30.00,33.80,37.30,40.7,43.9,
          46.9,49.7,52.4,54.9,57.3,
          59.6,62.0,64.3,66.7,71.5,
          76.3,82.0,91.0,100.0,106.0,
          113.0,120.0,128.0,137.0,147.0,
          160.0,175.0,196.0,222.0,486.0]
    Ks = np.interp(Temp, Ts, Ks)*10**-3
    return Ks

def getAirMuFromTempK(T):
    muV = np.array([18.13,21.74,25.73,29.28,32.5,35.47,38.25,40.85,43.32,45.66,47.88,50.01,52.06,54.04,55.96,57.82,59.61])*10**-6
    muT = np.array([20,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600])
    mu = np.interp(T-273,muT,muV)
    return mu

def calculateLf(Q, D):
    Lf = 0.235*Q**0.4 - 1.02*D
    return Lf

def calculateVolume(Qc):
    # Orloff 1982 (quoted as 0.1-0.7m diameter only)
    V = Qc/1200
    return V

def calculateFrequency(D,g=9.81):
    # Drysdale 2011 4.31
    f = (0.50)*((g/D)**0.5)
    return f


def calculateVirtualSource(Q,D):
    # Drysdale 2011 4.26
    #z0 = (-1.02+0.083*(Q**(2/5))/D)*D
    z0 = -1.02*D+0.083*(Q**(2/5))
    return z0


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    colors = getJHcolors()
    
    
    Q = 1000#1.5*7000 #13400*3 #35000
    Tamb = 20 + 273.15
    zs = np.linspace(0, 2, 101)
    T = calculateT(Q, zs, Tamb, g=9.81)
    V = calculateV(Q, zs, Tamb, g=9.81)
    
    Pr = getAirPrFromTempK(T)
    cp = getAirCpFromTempK(T)
    k = getAirKFromTempK(T)
    rho = getAirRhoFromTempK(T)
    nu = getAirNuFromTempK(T)
    mu = getAirMuFromTempK(T)
    alpha = k/(rho*cp)
    g = 9.81
    L = 1
    T1 = np.array(T)
    #T1[:] = 1500
    Ra = (2*g*abs(T1-Tamb)*(L**3))/((T1+Tamb)*alpha*nu)
    Re = rho*V*L/mu
    
    Nu_Ra = (0.825 + 0.324*Ra**(1/6))**2
    Nu_Re = (0.0296*Re**0.8)*(Pr**(1/3))
    
    h_Ra = (k/L)*Nu_Ra
    h_Re = (k/L)*Nu_Re
    
    lw = 3
    fs = 16
    plt.figure(figsize=(6,8))
    plt.plot(h_Ra, zs, '-', linewidth=lw, label='Natural', color=colors[0])
    plt.plot(h_Re, zs, '-', linewidth=lw, label='Forced', color=colors[1])
    
    plt.grid()
    plt.xlabel("Heat Transfer Coefficient (W/m$^{2}$K)", fontsize=fs)
    plt.ylabel("Height Above Fire (m)", fontsize=fs)
    plt.legend(fontsize=fs) #, bbox_to_anchor=(0.9, 0.8))
    plt.xlim(0, 14)
    plt.xticks([0, 2, 4, 6, 8, 10, 12, 14])
    plt.ylim(0, 2)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    plt.tick_params(labelsize=fs)
    #plt.title("HRR = %0.1f MW"%(Q/1000), fontsize=fs)
    plt.tight_layout()
    
    
    