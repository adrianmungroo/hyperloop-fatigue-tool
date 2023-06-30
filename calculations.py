import fractions
import stressfunctions as func
import numpy as np
from matplotlib import pyplot as plt
import rainflow
import pandas as pd

increments_per_second = 10 # 1 for seconds, 1000 for milliseconds

def failure(yearTempSwing,dayTempSwing,dx,internalP,ambientP,youngMod,density,ultimate,tubethick,tuberadius,podmass,period,speed,podtemp,tubetemp,radIndex,alpha,steelType):

    amb_p = ambientP # Ambient pressure outside the tube [Pa]

    internal_p = internalP # Pressure inside of the tube [Pa]

    thickness = tubethick/1000 # [m] Input as millimeters

    radius = tuberadius/100  # [m] Input as centimeters

    support_separation = dx  # This is the distance between any two supports in the system 100 [m]

    rho = density  # [kg/m3]

    alpha = alpha * 10**-6  # 11.6 This is the coefficient of thermal expansion https://www.engineeringtoolbox.com/linear-expansion-coefficients-d_95.html [m/(mK)] 

    mod_elasticity = youngMod * 10 ** 9  # The young's modulus of the material. This can be changed to simulate other materials [Pa] Input in Giga Pascals

    S_ut = ultimate * 10**6  # The ultimate stress of the material. This can be changed to simulate other materials [Pa] Input in Mega Pascals

    # emissivity = 1.0 #this is not used yet

    pod_mass = podmass  # [kg]

    fatigue_stength_fraction = func.fatigueStengthFraction(S_ut)  # Figure 6-18 Page 285 in Shigleys

    # damping_factor = 0.5 # A crude estimate of how much axial stress dampers will fractionally decrease the stress. Not used yet

    temperature_swing_year = yearTempSwing # Estimate of how much the temperature goes above and below its manufactured temperature during a year

    temperature_swing_day = dayTempSwing # Estimate of how much the temperature goes above and below its manufactured temperature during a day

    pod_speed = speed # Pod speed [m/s]

    mass = support_separation*func.distributed_load(rho,radius,thickness)/9.81 # Mass of tube [kg]

    pod_period = period*60 # The interval of time before the next pod arrives [mins]

    ##                                                                               CALCULATIONS
     
    # We find what will be a constant stress. A stress that will not change as long as the system variables stay the same
    # In this case this is the hoop stress, and axial pressure stress. Notice that these are in different directions. Will probably separate them

    constant = [func.axial_pressure_stress(internal_p, amb_p, radius,thickness) + func.bending_stress(radius, thickness, func.constant_bending_moment(rho,radius,thickness,support_separation)),  
                func.sig_theta(internal_p,amb_p,radius,thickness)]
    
    average_stress = constant[0] # The stress that acts at any given time  

    max_yearly_temperature_stress = func.thermal_stress(alpha, mod_elasticity, temperature_swing_year, radIndex=1) # Finds the maximum stress. Basically the alternating stress

    rms_yearly_temp_stress = np.sqrt(2)*max_yearly_temperature_stress # We will be assuming that yearly temp acts as a rms constant

    average_stress += rms_yearly_temp_stress # Adding onto the average stress for our assumption

    max_daily_temperature_stress = func.thermal_stress(alpha, mod_elasticity, temperature_swing_day, radIndex)

    min_daily_temperature_stress = -1*func.thermal_stress(alpha, mod_elasticity, temperature_swing_day, 1) # radIndex = 1 since it's at night

    daily_temperature_alternating_stress = (max_daily_temperature_stress - min_daily_temperature_stress)/2 

    daily_average_stress = (max_daily_temperature_stress + min_daily_temperature_stress)/2 # can now draw a sine wave since we have alternating and average

    # Here we will attempt to model the stress changes as the pod enters the tube segment. First we look at how temperature increases as the pod enters

    temperature_increase = func.convection(radius, thickness, support_separation, pod_speed, mass, tubetemp, podtemp) + func.radiation(radius, thickness, support_separation, pod_speed, mass, tubetemp, podtemp) 
    
    pod_thermal_stress_on_tube = func.thermal_stress(alpha, mod_elasticity, temperature_increase, 1)

    # Now we find the pod induced stress as the pod is in the tube. This is a sum of its bending stress and temperature increase stress

    max_pod_stress = func.bending_stress(radius,thickness,func.pod_moment(rho,radius,thickness,support_separation,pod_mass)) + pod_thermal_stress_on_tube

    min_pod_stress = 0 # When the pod isn't there, there's no stress

    time_spent_in_tube = support_separation / pod_speed

    ## Setting up data to plot a full stress history
    day_in_seconds = 86400 
    time = np.linspace(0,day_in_seconds,day_in_seconds*increments_per_second)


    # print("Generating stress signals...")
    pod_stress_signal = func.pod_stress_series(max_pod_stress,time_spent_in_tube,pod_period,increments_per_second)
    tube_stress_signal = daily_temperature_alternating_stress*np.sin(2*np.pi*time/day_in_seconds) + daily_average_stress
    stress_signal = pod_stress_signal + tube_stress_signal + average_stress

    # plt.plot(time,stress_signal)
    # plt.show()

    Sm = 0.75*S_ut #setting the y-intercept for the start of the high cycle loading. Originally 0.75. 8.715021

    ka = func.surfaceFactor(S_ut/1e6, steelType) #surface factor ka

    kc = 0.85 #loading factor

    T = tubetemp*(9/5) + 32

    kd = func.temperatureFactor(T) #temperature factor

    ke = 0.814 #reliability factor. Assumes a 99% reliability
    
    kf = 0.5 #miscellaneous factor. was 0.8

    modifier = ka*kc*kd*ke*kf

    S_e = S_ut*modifier #setting the endurance stress of the S-N curve. We can modify further, maybe multiply by half?

    # print("Starting rainflow analysis...")

    safety_factors = []
    cycles_to_fail = []
    rng_list = []
    mean_list = []
    # count_list = []

    accumulated_damage = 0
    c = 0

    for rng, mean, count, _, _ in rainflow.extract_cycles(stress_signal): 
        safety_factors.append(func.goldman_safety_factor(rng, mean, S_e, S_ut))
        cycles = func.cycles2Fail(rng,mean,S_ut,S_e,fatigue_stength_fraction)
        cycles_to_fail.append(cycles)
        accumulated_damage += count/cycles
        # c += 1
        rng_list.append(rng)
        mean_list.append(mean)
        # count_list.append(count)
    
    average_stress_range = np.mean(rng_list)
    average_mean_stress = np.mean(mean_list)
    average_safety_factor = round(func.goldman_safety_factor(average_stress_range, average_mean_stress, S_e, S_ut),5)
    
    # print(c)

    # data  = {'Stress Range' : rng_list,
    #          'Mean Stress' : mean_list,
    #          'Stress Count' : count_list,
    #          'Safety Factors' : safety_factors,
    #          'Cycles to Fail': cycles_to_fail}

    # df = pd.DataFrame(data)
    # df.to_excel(r'C:\Users\adria\Dropbox\Academic\Hyperloop\Python Code\successful\v22\export_dataframe.xlsx', index = False, header=True)

    if accumulated_damage < 1e-10:
        num_days = np.inf
    else:
        num_days = round(1/accumulated_damage,2)

    # mean_safety_factor = round(np.mean(safety_factors),3)
    mean_safety_factor = average_safety_factor

    return (mean_safety_factor, accumulated_damage, num_days)


# msf, acum_dmg, days = failure(30, 15, 100, 101, 101325, 200, 7860, 600, 30, 150, 5000, 1, 150, 200, 20, 1.5, 11.6, 'cold') #above ground
msf, acum_dmg, days = failure(10, 10,   150    ,101, 101325, 200, 7860, 600, 30,    150     ,5000, 1 ,150, 200, 10, 1.5, 11.6, 'cold') #underground
#update_(yearTempSwing,dayTempSwing,dxx,internalP,ambientP,youngMod,density,ultimate,tubethick,tuberadius,podmass,period,speed,podtemp,tubetemp,radIndex,alpha,steelType):

print(f"This system's daily stresses have a mean safety factor of {msf}. Once cycle accumulates damage of {acum_dmg/1440} so it will exist for {1/(acum_dmg/1440)} cycles")
