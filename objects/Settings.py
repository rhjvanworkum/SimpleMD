class Settings():
    def __init__(self, N=5, density=0.9, temperature=1, delta_t=0.005, step_limit=10000, step_avg=100,
                 step_adjust_temp=10, step_equilibrium=20, show_progress=False, show_summary=True):
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