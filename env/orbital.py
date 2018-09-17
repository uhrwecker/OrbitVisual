import matplotlib.pyplot as pl
from scipy.integrate import ode
import numpy as np

# r** = G*mj (rj - ri) / rij^3

class Orbital():
    def __init__(self, k, tof=1, nbsteps=1000, rtol=1e-10, masses=[1000, 100]):
        self.k = k
        self.tof = tof
        self.nbsteps = nbsteps
        self.rtol = rtol
        self.masses = masses

    def cowell(self, r0, v0, ad=None):
        """Propagates orbit using Cowell's formulation.
        Parameters
        ----------
        k : float
            Gravitational constant of main attractor (km^3 / s^2).
        r0 : array
            Initial position (km).
        v0 : array
            Initial velocity (km).
        ad : function(t0, u, k), optional
             Non Keplerian acceleration (km/s2), default to None.
        tof : float
            Time of flight (s).
        nsteps : int, optional
            Maximum number of internal steps, default to 1000.
        rtol : float, optional
            Maximum relative error permitted, default to 1e-10.
        Raises
        ------
        RuntimeError
            If the algorithm didn't converge.
        Notes
        -----
        This method uses a Dormand & Prince method of order 8(5,3) available
        in the ``scipy.integrate.ode`` module.
        """
        k = self.k
        nsteps = self.nbsteps
        tof = self.tof
        rtol = self.rtol
        x, y, z = r0
        vx, vy, vz = v0
        u0 = np.array([x, y, z, vx, vy, vz])

        # Set the non Keplerian acceleration
        if ad is None:
            ad = lambda t0, u, k: (0, 0, 0)

        # Set the integrator
        rr = ode(self.func_twobody).set_integrator('dop853', rtol=rtol, nsteps=nsteps)
        rr.set_initial_value(u0)  # Initial time equal to 0.0
        rr.set_f_params(k, ad)  # Parameters of the integration

        # Make integration step
        rr.integrate(tof)

        if rr.successful():
            r, v = rr.y[:3], rr.y[3:]
        else:
            raise RuntimeError("Integration failed")

        return r, v

    def cowell_multibody(self, r0, v0, other_r):
        k = self.k
        nsteps = self.nbsteps
        tof = self.tof
        rtol = self.rtol
        masses = self.masses
        
        x, y, z = r0
        vx, vy, vz = v0
        u0 = np.array([x, y, z, vx, vy, vz])

        rr = ode(self.func_multibody).set_integrator('dop853', rtol=rtol, nsteps=nsteps)
        rr.set_initial_value(u0)
        rr.set_f_params(k, masses, other_r)

        rr.integrate(tof)

        if rr.successful():
            r, v = rr.y[:3], rr.y[3:]
        else:
            raise RuntimeError("Integration failed")

        return r, v

    @staticmethod
    def func_multibody(t0, u, k_list, m_list, other_u_list):
        x, y, z, vx, vy, vz = u

        sumx = 0
        sumy = 0
        sumz = 0
        for index in range(0, len(other_u_list)):
            distance = ((x-other_u_list[index][0])**2 +
                        (y-other_u_list[index][1])**2 +
                        (z-other_u_list[index][2])**2)**1.5
            sumx += (k_list[index]*m_list[index]*(other_u_list[index][0]-x)) / distance**3
            sumy += (k_list[index]*m_list[index]*(other_u_list[index][1]-y)) / distance**3
            sumz += (k_list[index]*m_list[index]*(other_u_list[index][2]-z)) / distance**3

        du = np.array([
            vx,
            vy,
            vz,
            sumx,
            sumy,
            sumz
        ])
        return du

    @staticmethod
    def func_twobody(t0, u, k, ad):
        """Differential equation for the initial value two body problem.
        This function follows Cowell's formulation.
        Parameters
        ----------
        t0 : float
            Time.
        u : ndarray
            Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
        k : float
            Standard gravitational parameter.
        ad : function(t0, u, k)
             Non Keplerian acceleration (km/s2).
        """
        ax, ay, az = ad(t0, u, k)

        x, y, z, vx, vy, vz = u
        r3 = (x**2 + y**2 + z**2)**1.5

        du = np.array([
            vx,
            vy,
            vz,
            -k * x / r3 + ax,
            -k * y / r3 + ay,
            -k * z / r3 + az
        ])
        return du


def main():
    o1 = Orbital([1000000], masses=[100000000000])
    o2 = Orbital([100000], masses=[100000000000])

    u1 = [100, 0, 0, 0.1, 0.1, 0]
    u2 = [-50, 0, 0, 0.1, 0.1, 0]

    r1, v1 = o1.cowell_multibody(u1[0:3], u1[3:], [u2])
    r2, v2 = o2.cowell_multibody(u2[:3], u2[3:], [r1+v1])

    print(r1, v1)
    print(r2, v2)

    x1 = [r1[0]]
    y1 = [r1[1]]
    x2 = [r2[0]]
    y2 = [r2[1]]

    for i in range(0, 12):
        r1, v1 = o1.cowell_multibody(r1, v1, [r2+v2])
        r2, v2 = o2.cowell_multibody(r2, v2, [r1+v1])
        print(r1, r2)
        print(v1, v2)

        x1.append(r1[0])
        x2.append(r2[0])
        y1.append(r1[1])
        y2.append(r2[1])


    pl.plot(x1, y1)
    pl.plot(x2, y2)
    pl.grid()
    pl.show()

if __name__ == '__main__':
    main()

