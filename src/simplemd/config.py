"""Simulation configuration."""


class Settings:
    """Container for the parameters of a simulation run.

    Args:
        N: Number of unit cells per dimension (the system is ``N**3`` cells).
        density: Number density in reduced (dimensionless) units.
        temperature: Target temperature in reduced units.
        delta_t: Integration time step.
        step_limit: Total number of integration steps to run.
        step_avg: Interval (in steps) over which properties are averaged/printed.
        step_adjust_temp: Interval at which velocities are rescaled to ``temperature``.
        step_equilibrium: Reserved for future use (currently unused).
        show_progress: Print percentage progress while running.
        show_summary: Print a property summary every ``step_avg`` steps.
        seed: Optional RNG seed. When set, the run is fully reproducible. When
            ``None`` (the default) the global RNGs are left untouched, preserving
            the original non-deterministic behavior.
    """

    def __init__(
        self,
        N=5,
        density=0.9,
        temperature=1,
        delta_t=0.005,
        step_limit=10000,
        step_avg=100,
        step_adjust_temp=10,
        step_equilibrium=20,
        show_progress=False,
        show_summary=True,
        seed=None,
    ):
        self.N = N
        self.density = density
        self.temperature = temperature
        self.delta_t = delta_t

        self.step_limit = step_limit
        self.step_avg = step_avg
        self.step_adjust_temp = step_adjust_temp
        self.step_equilibrium = step_equilibrium

        self.show_progress = show_progress
        self.show_summary = show_summary

        self.seed = seed
