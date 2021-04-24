from models.gas_particles import gas_particles
from models.TIP4P import TIP4P

from objects.Settings import Settings
from reporter import Reporter

from integrators.LeapFrogIntegrator import LeapFrogIntegrator
from integrators.VerletVelocityIntegrator import VerletVelocityIntegrator

from interaction_methods.pair_by_pair import ComputeForcesPBP
from interaction_methods.cell_division import ComputeForcesCD

from flask import Flask, request
import json
import numpy as np

app = Flask(__name__)

system_dict = {
    'LJ': gas_particles,
    'TIP4P': TIP4P
}

integrator_dict = {
    'leapfrog': LeapFrogIntegrator,
    'verlet': VerletVelocityIntegrator
}

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

@app.route('/api/run_job', methods=['POST'])
def run_job():

    data = request.json

    settings = Settings(
        N=data['N'],
        density=data['density'],
        temperature=data['temperature'],
        delta_t=data['delta_t'],
        step_limit=['steps'],
        step_avg=100,
        step_adjust_temp=10,
        step_equilibrium=20,
        show_progress=False,
        show_summary=True,
    )

    if data['system'] == 'LJ':
        system = system_dict[data['system']](settings, integrator_dict[data['integrator']], ComputeForcesCD)
    else:
        system = system_dict[data['system']](settings, integrator_dict[data['integrator']], ComputeForcesPBP)

    # add reporter
    system.add_reporter(Reporter, 1, data['system'])

    # run the program
    system.simulate()

    return json.dumps({'output': system.reporter.traj}, cls=NumpyEncoder)

if __name__ == '__main__':
    app.run(port=5000, debug=True, use_reloader=False)