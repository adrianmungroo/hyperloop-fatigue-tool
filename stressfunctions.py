import numpy as np
def distributed_load(rho,radius,thickness): 
    '''Calculates the weight load per unit length of the tube itself'''
    if rho == 0:
        return 1e-32
    q = rho*np.pi*(((radius+thickness)**2) - (radius**2)) * 9.81  # N/m
    return q

def constant_bending_moment(rho,radius,thickness,dx):
    """egg"""
    return distributed_load(rho,radius,thickness)*dx**2/8

def pod_moment(rho,radius,thickness,dx,pod_mass):
    '''Finds the bending moment with the pod assumed to be in the center of the tube'''
    return pod_mass*9.81*dx/4  # Nm

def sig_theta(internal_p,amb_p,radius,thickness):
    '''Finds the hoop stress in the tube'''
    if (2*radius/thickness) > 20:
        return (internal_p-amb_p)*radius/thickness  # Pa
    sig_theta = axial_pressure_stress(internal_p,amb_p,radius,thickness) - (radius**2 * (radius+thickness)**2 * (amb_p-internal_p)/(radius**2 * ((radius+thickness)**2 - radius**2)))
    return sig_theta

def inertia(radius,thickness):
    '''Finds the moment of inertia of the tube'''
    return np.pi/4*((radius+thickness)**4 - radius**4)

def thermal_stress(alpha,mod_elasticity,temperature_change,radIndex):
    '''Finds the stress due to the change in temperature'''
    if 0 in {alpha,mod_elasticity,temperature_change,radIndex}:
        return 1e-32
    return alpha*mod_elasticity*temperature_change*radIndex

def bending_stress(radius,thickness,moment):
    '''Finds the stress due to bending moment of the pod passing through the tube'''
    if radius == 0:
        return 1e-32
    if inertia(radius,thickness) == 0:
        return 1e32
    return moment*(radius+thickness)/inertia(radius,thickness)

def axial_pressure_stress(internal_p,amb_p,radius,thickness):
    '''finds the axial stress caused by pressure'''
    #return (internal_p-amb_p)*radius/(2*thickness)
    if thickness == 0:
        return 1e32
    if (2*radius / thickness) > 20:
        return (internal_p-amb_p)*radius/(2*thickness)
    return (internal_p*radius**2 - amb_p*(radius+thickness)**2)/((radius+thickness)**2 - radius**2)

# def sig_axial(internal_p,amb_p,rho,radius,thickness,dx,pod_mass,alpha,mod_elasticity,temperature_change):
#     " Finds the total axial stress according to the component stresses"
#     sig_axial = axial_pressure_stress(internal_p,amb_p,radius,thickness) + bending_stress(rho,radius,thickness,dx,pod_mass) + thermal_stress(alpha,mod_elasticity,temperature_change)
#     return sig_axial

def von_mises(self): # find where this formula comes from? Not used anymore but kept as legacy code
    sigma_e = np.sqrt((self.sig_theta()**2+self.sig_axial()**2+(self.sig_axial()-self.sig_theta())**2)/2)
    return sigma_e

def convection(radius,thickness,dx,speed,mass,tubetemp,podtemp):
    '''Function that gives the resultant increase in temperature from convection as the pod enters the slice of tube'''
    ri = radius
    ro = ri+thickness
    pr = 0.7  #prandtl number
    beta = 0.00285 # K^-1
    visc = 20.92*1e-6 # m^2 / s
    k = 0.03 #W/mK
    alpha = 29.9*1e-6 #m^2 / s
    To = tubetemp #35 celcius. ambient temperature of the tube
    Ti = podtemp #100 celcius temperature of the pod
    length = dx #meters. Length of pod
    c = 490 # joules / kilogram * kelvin
    time = dx/speed #second. Length of time that the pod is present to transfer heat
    g = 9.81
    L_c = (2*(np.log(ro/ri))**(4/3))/((ri**(-3/5)+ro**(-3/5))**(5/3)) #meters
    Ra = g*beta*(Ti-To)*L_c**3.0 / (visc*alpha) #assumes thta Ti - To > 0. 
    if Ra < 0:
        keff = -1*k*0.74*(pr/(0.861+pr))**(1/4) * (-1*Ra)**(1/4) #W/mK
    else:
        keff = k*0.74*(pr/(0.861+pr))**(1/4) * Ra**(1/4) #W/mK
    q_dot = 2*np.pi*L_c*keff*(Ti-To)/(np.log(ro/ri)) #W/m
    heat_power = q_dot * length
    dT = heat_power*time/(c*mass)
    return dT

def radiation(radius,thickness,dx,speed,mass,tubetemp,podtemp):
    '''Function that gives the resultant increase in temperature from radiation as the pod enters the slice of tube'''
    L = dx #meters. Length of pod at any point in time
    ri = radius #m. Radius of pod
    ro = ri+thickness #m. Radius of tube
    epi = 0.07 #polished steel
    epo = 0.2 #outer tube steel emmisivity
    Ti = podtemp +273 #Temperature of pod
    To = tubetemp +273 #Temperature of tube
    c = 490 # joules / kilogram * kelvin
    time = dx/speed #second. Length of time that the pod is present to transfer heat
    A = np.pi*2*ri*L #Surface area of pod approximated as cylinder
    sigma = 5.67*1e-8 #W / (m^2 * K^4)
    q12 = sigma*A*(Ti**4 - To**4)/(1/epi + ri*(1-epo)/(ro*epo))
    dT = q12*time/(c*mass)
    return dT

def cycles2Fail(alternating_stress,mean_stress,S_ut,S_e,fatigue_stength_fraction):
    '''Finds the cycles to failure for the given stress'''

    safety_factor = goldman_safety_factor(alternating_stress, mean_stress, S_e, S_ut)

    # if safety_factor > 1:
    #     return np.inf
    # elif safety_factor == 0:
    #     return 0
    # elif safety_factor < 0:
    #     return None
    # else:
    #     sigma_rev = alternating_stress / (1-mean_stress/S_ut)  # finding reversible stress approximation
    #     a = ((fatigue_stength_fraction*S_ut)**2)/S_e
    #     b = -(1/3)*np.log10(fatigue_stength_fraction*S_ut/S_e)
    #     N = (sigma_rev/a)**(1/b)
    #     return N

    sigma_rev = alternating_stress / (1-mean_stress/S_ut)  # finding reversible stress approximation
    a = ((fatigue_stength_fraction*S_ut)**2)/S_e
    b = -(1/3)*np.log10(fatigue_stength_fraction*S_ut/S_e)
    N = (sigma_rev/a)**(1/b)
    return N

def bucklingCheck(force,radius,thickness,modElas,length):
    ''' Checks to see if the tube will buckle under compressive loads '''
    critical = 4*(3.14159**2)*modElas*inertia(radius,thickness)/(length**2)
    if force < -1*critical:
        return True
    return False

def fatigueStengthFraction(sut):
    '''calculates the fatigue strength fraction based on page 285 in Shigleys. Do not extrapolate too far'''
    sut = sut*0.00014503773/1000 #convert to kpsi
    return -5.501*10**(-8)*sut**3 + 2.887*10**-5*sut**2 - 5.535*10**-3*sut + 1.164

def surfaceFactor(sut,steel):
    ''' return the endurance stress modifying factor relating to the surface quality of the material.
        Assumes that the steel is cold rolled'''
    if steel == 'cold':
        return 4.51*(sut**(-0.265))
    elif steel == 'hot':
        return 57.7*(sut**(-0.718))

def temperatureFactor(T):
    ''' return the endurance stress modifying factor relating to the surface quality of the material.
        Assumes that the steel is cold rolled'''
    return 0.975 + (0.432E-3)*T - (0.115E-5)*T**2 + (0.104E-8)*T**3 - (0.595E-12)*T**4

def pod_stress_series(max_pod_stress,time_spent_in_tube,period,increments_per_second):
    total = 86400
    time = np.linspace(0,total,total*increments_per_second)

    signal = []
    for t in time:
        if  t%period <= time_spent_in_tube: 
           signal.append(max_pod_stress*np.sin(2*np.pi*(t%period)/(2*time_spent_in_tube)))
        else:
            signal.append(0)
    
    return signal

def goldman_safety_factor(alternating_stress, mean_stress, S_e, S_ut):
    return 1/(alternating_stress/S_e + mean_stress/S_ut)