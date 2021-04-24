import numpy as np
import json


class Reporter:

    def __init__(self, parent, step_report_avg, model, output_file=None):
        self.system = parent
        self.step_report_avg = step_report_avg
        self.model = model

        if (model == "LJ"):
            self.traj = np.zeros((int(self.system.step_limit / step_report_avg), self.system.n_mol, self.system.n_dim))
        elif (model == "TIP4P"):
            self.traj = np.zeros((int(self.system.step_limit / step_report_avg), self.system.n_mol * 3, self.system.n_dim))

        self.i = 0

        self.output_file = output_file

    def step(self):
        if self.system.step_count % self.step_report_avg == 0:
            if (self.model == "LJ"):
                self.report()
            elif (self.model == "TIP4P"):
                self.report_water()

            self.i += 1

    def report(self):
        for n in range(self.system.n_mol):
            self.traj[self.i][n] = self.system.mol[n].r

    def report_water(self):
        i = 0
        for n in range(self.system.n_mol):
            for j in range(self.system.n_sites):
                if j != 1:
                    self.traj[self.i][i] = self.system.sites[self.system.n_sites * n + j].r
                    i += 1

    def export(self):
        with open(self.output_file, 'w') as f:
            json.dump(self.traj.tolist(), f)