# SimpleMD
A simple Molecular Dynamics project, including the LJ soft fluid model and TIP-4P water model

### How to use
1. Go into main.py
2. Change the Settings object to the settings you prefer:
  - N -> dimensionality of the space - keep this one 3
  - density -> density of the system in dimensionless units
  - temperature -> temperature of the system in dimensionless units
  - step_limit -> the amount of timesteps
  - step_avg -> the interval of steps after you which a properties are evaluated and a summary is printed
  - step_adjust_temp -> the interval of steps after which the particle velocities are scaled down to the right temperature
  - step_equilibrium -> currently not used
  - show_progress -> prints out the progress of the calculation
  - show_summary -> prints out the summary of the properties of the system
3. Choose a system:
  - gas_particles: a Lennard-Jones soft fluid simulation 
  - TIP4P: a TIP-4P simulation of water molecules
4. And Integrator & Interaction Method:
  - gas_particles is compatible with both LeapFrog and VerletVelocity Intergrator and pair_by_pair & cell_division Interaction methods
  - TIP4P is currently only working with LeapFrogIntegrator and pair_by_pair interaction method
5. Add a Reporter, if you would want to visualize the results afterwards
  - args for the reporter are: Reporter object, amount of steps after which to report, the model(either "LJ" or "TIP4P) and the output file location
6. Run main.py

### How to visualize
The results of the simulation will get stored as a numpy array in a .json file, which can be visualized afterwards on a webpage on my personal site: (coming soon).
