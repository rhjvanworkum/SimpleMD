from models.gas_particles import gas_particles
from models.TIP4P import TIP4P
from objects.Settings import Settings
from integrators.LeapFrogIntegrator import LeapFrogIntegrator
from integrators.PredictorCorrectorIntegrator import PredictorCorrectorIntegrator
from integrators.VerletVelocityIntegrator import VerletVelocityIntegrator
from integrators.RungeKuttaIntegrator import RungeKuttaIntegrator
from interaction_methods.pair_by_pair import ComputeForcesPBP
from reporter import Reporter

settings = Settings(
    N=3,
    density=0.8,
    temperature=2,
    delta_t=0.005,
    step_limit=1000,
    step_avg=100,
    step_adjust_temp=10,
    step_equilibrium=20,
    show_progress=False,
    show_summary=True,
)

system = gas_particles(settings, LeapFrogIntegrator, ComputeForcesPBP)

# add reporter
system.add_reporter(Reporter, 1, "LJ", "output/test_3d_3_chain.json")

# run the program
system.simulate()