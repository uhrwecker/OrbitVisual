from env.orbital import Orbital
import numpy as np

import itertools
import copy

class OrbitalEnvironment():
    def __init__(self, masses, r_list, v_list):
        self.states = []
        for r, v in zip(r_list, v_list):
            self.states.append((r, v))
        self.masses = masses
        self.permutated_masses = list(itertools.permutations(self.masses))
        self.k = [6e-11 for ind in range(len(masses))]
        self.orbitals = []
        if len(masses) > 1:
            for mass_list in self.permutated_masses:
                self.orbitals.append(Orbital(self.k[:-1], masses=list(mass_list[1:])))
        else:
            self.orbitals.append(Orbital(self.k))

    def add_orbit(self, r, v, mass):
        for orbit in self.orbitals:
            orbit.masses.append(mass)
            orbit.k.append(6e-11)
        self.orbitals.append(Orbital(self.k[0], masses=self.masses))
        self.masses.append(mass)
        self.k.append(6e-11)
        self.states.append((r, v))

    def step(self):
        new_states = []
        for n, state in enumerate(self.states):
            r = state[0]
            v = state[1]

            rest_states = copy.deepcopy(np.array(self.states).tolist())
            rest_states.remove(np.array(state).tolist())
            rest_states = [state[0]+state[1] for state in rest_states]
            
            new_states.append(self.orbitals[n].cowell_multibody(r, v, rest_states))

        self.states = new_states#

        return self.states
                

def main():
    o = OrbitalEnvironment([1e24, 1e24], [[100, 0, 0], [0, 0, 0]], [[1, 0, 0], [0, 0, 0]])
    o.step()
    print(o.step())

if __name__ == '__main__':
    main()
