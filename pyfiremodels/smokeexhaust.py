import numpy as np

def maximumVolumetricFlow(d, Ts, Tamb, dwall, width, height, units='imperial'):
    Di = 2*width*height/(width + height)
    d_over_Di = d / Di
    
    if d_over_Di < 2:
        print("Error d_over_Di must be greater than 2, but was %0.1f"%(d_over_Di))
    else:
        print("d_over_Di = %0.1f"%(d_over_Di))
    if dwall > 2*Di:
        gamma = 1.0
    else:
        gamma = 0.5
    
    if units == 'imperial':
        prefactor = 452
    elif units == 'metric':
        prefactor = 4.16
    v_max = prefactor*gamma*(d**(5/2)) * (((Ts - Tamb) / Tamb)**(1/2))
    return v_max

def minimumSeparationDistance(V, units='imperial'):
    if units == 'imperial':
        prefactor = 0.065
    elif units == 'metric':
        prefactor = 0.9
    S_min = prefactor * (V ** (1/2))
    return S_min

if __name__ == "__main__":
    d = 30
    Ts = 80 + 491.67 # F to R
    Tamb = 70 + 491.67 # F to R
    dwall = 0
    width = 5
    height = 5
    
    v_max = maximumVolumetricFlow(d, Ts, Tamb, dwall, width, height, units='imperial')
    print("Max cfm = %0.1f ft3/min"%(v_max))
    
    n_vents = 2
    
    S_min = minimumSeparationDistance(35000)
    
    print(S_min)
    
    
    velocity = v_max / (n_vents * width * height)
    
    print("Velocity = %0.1f ft/min"%(velocity))
    
    