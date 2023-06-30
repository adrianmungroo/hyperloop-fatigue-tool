from cmath import nan
import numpy as np
import matplotlib.pyplot as plt
import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
import plotly.graph_objs as go
import stressfunctions as func
import calculations as calc_main
import rainflow

## stylesheet here 

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

#layout of online app here

app = dash.Dash(__name__, external_stylesheets= external_stylesheets)
server = app.server

app.layout = html.Div([
    html.Div([
        html.Div([
             dcc.Graph(id = 'goldman-graph'),       
        ], className="row"),
        html.Div([
             dcc.Graph(id = 'sn-curve'),       
        ], className="row"),
    ], className = "six columns"),

    html.Div([
        html.Div([
            html.Div([ #column for tube design variables

                html.H5('Tube Material Properties'),

                html.Label('Density (kg/m^3)'),
                dcc.Input(id = 'density', value=7860, type = 'number'),

                html.Label('Modulus of Elasticity (GPa)'),
                dcc.Input(id = 'youngMod', value=200, type = 'number'),

                html.Label('Ultimate Tensile Strength (MPa)'),
                dcc.Input(id = 'ultimate', value=600, type = 'number'),

                html.Label('Thermal Expansion Coefficient (m/mK * 10^-6)'),
                dcc.Input(id = 'alpha', value=11.6, type = 'number'),

                html.Label('Type of steel used'),
                dcc.Dropdown(
                    id='steelType',
                    options=[
                        {'label': 'Cold Rolled', 'value': 'cold'},
                        {'label': 'Hot Rolled', 'value': 'hot'},
                    ],
                    value='cold'
                ),

            ],className = 'three columns'),

            html.Div([ #column for tube design variables
                html.H5('Tube Design Variables'),

                html.Label('Thickness (mm)'),
                dcc.Input(id = 'tubethick', value=30, type = 'number'),

                html.Label('Radius (cm)'),
                dcc.Input(id = 'tuberadius', value=150, type = 'number'),

                html.Label('Internal Pressure (Pa)'),
                dcc.Input(id = 'internalP', value=101, type = 'number'),

                html.Label('Tube resting temperature (Celcius)'),
                dcc.Input(id = 'tubetemp', value=20, type = 'number'),

                html.Label('Distance between piers (m)'),
                dcc.Input(id = 'dx', value=100, type = 'number'),

            ],className = 'three columns'),

            html.Div([ #column for environmental variables
                html.H5('Environment Variables'),

                html.Label('Yearly Temperature Swing (Celcius)'),
                dcc.Input(id = 'year-temp-slider', value=12, type = 'number'),

                html.Label('Daily Temperature Swing (Celcius)'),
                dcc.Input(id = 'day-temp-slider', value=5, type = 'number'),

                html.Label('Ambient Pressure (Pa)'),
                dcc.Input(id = 'ambientP', value=101325, type = 'number'),

                html.Label('Radiation Index'),
                dcc.Input(id = 'rad', value=1.5, type = 'number'),

            ],className = 'three columns'),

            html.Div([ #column for pod variables
                html.H5('Pod Variables'),

                html.Label('Mass (kg)'),
                dcc.Input(id = 'podmass', value=5000, type = 'number'), #originally 2578.21903108
                
                html.Label('Speed (m/s)'),
                dcc.Input(id = 'speed', value=150, type = 'number'),

                html.Label('Surface Temperature (Celcius)'),
                dcc.Input(id = 'podtemp', value=200, type = 'number'),

                html.Label('Operation interval (minutes)'),
                dcc.Input(id = 'period', value=1, type = 'number'),  

            ],className = 'three columns'),

        ], className="row"),
        html.H4(['Cycles to failure'],style={'textAlign': 'center'}),
        html.Div(
            id='failure-text', className = "row"
        ),
    ], className="six columns"),
])


@app.callback(
    [
        Output('sn-curve','figure'),
        Output('goldman-graph','figure'),
        Output('failure-text','children'),
    ],
    [
        Input('year-temp-slider','value'),
        Input('day-temp-slider','value'),
        Input('dx','value'),
        Input('internalP','value'),
        Input('ambientP','value'),
        Input('youngMod','value'),
        Input('density','value'),
        Input('ultimate','value'),
        Input('tubethick','value'),
        Input('tuberadius','value'),
        Input('podmass','value'),
        Input('period','value'),
        Input('speed','value'),
        Input('podtemp','value'),
        Input('tubetemp','value'),
        Input('rad','value'),
        Input('alpha','value'),
        Input('steelType','value')
    ]
)


def update_(yearTempSwing,dayTempSwing,dxx,internalP,ambientP,youngMod,density,ultimate,tubethick,tuberadius,podmass,period,speed,podtemp,tubetemp,radIndex,alpha,steelType):
    increments_per_second = 10
    amb_p = ambientP # Ambient pressure outside the tube [Pa]

    internal_p = internalP # Pressure inside of the tube [Pa]

    thickness = tubethick/1000 # [m] Input as millimeters

    radius = tuberadius/100  # [m] Input as centimeters

    support_separation = dxx  # This is the distance between any two supports in the system 100 [m]

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

    for rng, mean, count, _, _ in rainflow.extract_cycles(stress_signal): 
        safety_factors.append(func.goldman_safety_factor(rng, mean, S_e, S_ut))
        cycles = func.cycles2Fail(rng,mean,S_ut,S_e,fatigue_stength_fraction)
        cycles_to_fail.append(cycles)
        accumulated_damage += count/cycles
        rng_list.append(rng)
        mean_list.append(mean)
        # count_list.append(count)
    
    average_stress_range = np.mean(rng_list)
    average_mean_stress = np.mean(mean_list)
    average_safety_factor = round(func.goldman_safety_factor(average_stress_range, average_mean_stress, S_e, S_ut),5)
    

    if accumulated_damage < 1e-10:
        num_days = np.inf
    else:
        num_days = round(1/accumulated_damage,2)

    # mean_safety_factor = round(np.mean(safety_factors),3)
    mean_safety_factor = average_safety_factor

    failure = f"Mean safety factor = {round(mean_safety_factor,3)}, Daily accumulated Damage = {round(accumulated_damage,6)}, Days till failure = {round(num_days,3)}"
    if str(num_days) == "nan":
        failure = "The structure will fail transiently, will not withstand even 1 cycle."


    sn_y = [np.log10(Sm),np.log10(S_e),np.log10(S_e)]
    sn_x = [10**3,10**6,10**8]

    sn_curve = {'data':[
        {'x': sn_x, 'y':sn_y,'type' :'line', 'name': 'Max Allowable Stress for Given Cycle'},
        {'x': sn_x, 'y': np.ones(len(sn_x))*np.log10(average_stress_range), 'type': 'scatter', 'name': 'Mean System Alternating Stress'}
    ],'layout': dict(
        xaxis = {'type': 'log','title': '# of Cycles'},
        yaxis = {'type': 'log','title': 'Log of Stress (log10(MPa))'},
        transition = {'duration' : 750},
        title = "S-N Curve"
    )}

    goldman_graph = {'data':[
            {'x': [0,S_ut],'y':[S_e,0],'type' :'line','name': 'Goldman Line'},
            {'x': [average_mean_stress],'y':[average_stress_range,0],'name': 'Mean safety factor point'}
        ],'layout': dict(
            xaxis = {'type': 'linear','title': 'Average Stress (Pa)'},
            yaxis = {'type': 'linear','title': 'Alternating Stress (Pa)'},
            transition = {'duration' : 1000},
            title = "Goldman Graph"
        )}

    return sn_curve, goldman_graph, failure,

## run app server here
 
if __name__ == '__main__':
    app.run_server(debug=True,dev_tools_ui=False,dev_tools_props_check=False)
    server = app.server