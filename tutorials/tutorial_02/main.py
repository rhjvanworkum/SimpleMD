import numpy as np
import json

class Particle():
    def __init__(self):
        self.x = np.zeros(3)
        self.v = np.zeros(3)
        self.a = np.zeros(3)

# number of unit cells used on the x, y and z axis
n_cells = 3
n_unit_cells = np.array([n_cells, n_cells, n_cells])

# number of particles in the simulation
n_particles = n_unit_cells.prod()
particles = [Particle() for i in range(n_particles)]


density = 0.8
# cell_size is the array containing the length of each unit cell
cell_size = 1 / np.sqrt(density)
# region is the array containing the spatial dimenions of the region simulated
region = n_unit_cells * cell_size

dt = 0.005
n_steps = 1000
temperature = 1

r_cut = np.power(2, (1/6))

trajectory = np.zeros((n_steps, n_particles, 3))

## initializing the positions
i = 0
n_x, n_y, n_z = n_unit_cells[0], n_unit_cells[1], n_unit_cells[2]
for x in range(n_x):
    for y in range(n_y):
        for z in range(n_z):
            # put the position in the middle of the unit cell
            position = np.array([x + 0.5, y + 0.5, z + 0.5])
            # scale the position with the cell size dimensions
            position *= cell_size
            # center the position around the origin (0, 0, 0)
            position -= 0.5 * region

            particles[i].x = position

            i += 1

# initializing the velocities
vel_magnitude = np.sqrt(3 * (1 - 1 / n_particles) * temperature)
v_sum = 0

for particle in particles:
    particle.v = np.random.random_sample(3)
    particle.v *= vel_magnitude
    v_sum += particle.v

for particle in particles:
    particle.v -= (v_sum / n_particles)

# the accelerations can stay zero

# now we start the simulation
n = 0
while (n < n_steps):

    for particle in particles:
        print(particle.x)

    # compute the forces acting on the system
    compute_forces()

    # update velocity and position, apply boundary conditions
    i = 0
    for particle in particles:
        particle.v += particle.a * dt
        particle.x += particle.v * dt

        wrap_around(particle.x)

        trajectory[n][i] = particle.x

        i += 1

    n += 1

# export trajectory
with open('test.json', 'w') as f:
    json.dump(trajectory.tolist(), f)