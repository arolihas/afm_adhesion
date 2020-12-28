from math import pi, sqrt, pow, sin
import numpy as np
import csv 
import pandas as pd 
import progressbar # pip install progressbar if this does not work

wave_points = 1e6
spring_constant = 0.5 # N/m units
resonant_frequency = 25000 
quality_factor = 3
temperature = 300 

mass = spring_constant / pow((2 * pi * resonant_frequency), 2)
b = sqrt(spring_constant * mass) / quality_factor
sample_points_period = 500
dt = 1 / (sample_points_period * resonant_frequency)
dt = 5e-7
total_time = dt * wave_points

equilibrium_begin_pos = 5e-9
equilibrium_end_pos = 3e-10
equilibrium_distance = equilibrium_end_pos - equilibrium_begin_pos
equilibrium_pos = equilibrium_begin_pos

x = np.linspace(equilibrium_end_pos, equilibrium_begin_pos, 5000)
simulation_force = []
for i in x:
    val = 1e-11/i
    sim_value = 100 * pow(val, 6) - 0.0006 * pow(val, 3)
    simulation_force.append(sim_value)

driving_frequency = resonant_frequency
driving_amplitude = 0 # at resonance 5e-9
driving_fraction = driving_frequency / resonant_frequency
driving_force = driving_amplitude * spring_constant * driving_fraction * sqrt(pow((1/driving_fraction - driving_fraction), 2) + pow(1/quality_factor, 2))

velocity = 0
acceleration = 0
position_i = equilibrium_pos

eqb_pos = [0] * int(wave_points)
force_profile = [0] * int(wave_points)
positions = [0] * int(wave_points)
scaled_position = [0] * int(wave_points)
impacts = [0] * int(wave_points/sample_points_period)

velocity_turning = 0
oversampling = 10
oversampling_loop = 0
i = j = 0


with progressbar.ProgressBar(max_value=wave_points) as bar:

    while i < (wave_points - 1):
        initial_velocity = velocity

        min_index = np.argmin(abs(x - position_i))
        forcetrace_temp = simulation_force[min_index]
        
        equilibrium_pos = equilibrium_begin_pos + (equilibrium_distance * j/(wave_points * oversampling))
        noise = np.random.normal(0, sqrt(temperature * 1.38e-23 * spring_constant / pi * quality_factor * dt * resonant_frequency))
        F_drive = driving_force * sin(2 * pi * j * dt * driving_frequency)
        acceleration = 1 / mass * (-1 * spring_constant * (position_i - equilibrium_pos) - b * velocity + F_drive + forcetrace_temp + noise)
        
        velocity += (acceleration * dt)
        position_i += (velocity * dt)
        j += 1

        if (velocity*initial_velocity < 0 and velocity > 0):
            velocity_turning += 1
            # impacts[velocity_turning] = forcetrace_temp

        if (oversampling_loop == oversampling):
            i += 1
            positions[i] = position_i
            eqb_pos[i] = equilibrium_pos
            force_profile[i] = forcetrace_temp
            oversampling_loop = 0
        oversampling_loop += 1 
        bar.update(i)

for p in range(1,len(positions)):
    scaled_position[p]  = positions[p] - (equilibrium_distance * p/wave_points + equilibrium_begin_pos)

eqb_pos = np.asarray([p * 1e9 for p in eqb_pos])
scaled_position = np.asarray([sp * 1e9 for sp in scaled_position])

eqb_pos_filename = 'pos_k0.5_T{}_DFS.csv'.format(temperature)
scaled_position_filename = 'defl_k0.5_T{}_DFS.csv'.format(temperature)

pd.DataFrame(eqb_pos).to_csv(eqb_pos_filename, header=None)
pd.DataFrame(scaled_position).to_csv(scaled_position_filename, header=None)