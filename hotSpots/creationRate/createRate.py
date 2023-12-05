
import numpy as np

mu = 1.2566*10**(-6)
mp = 1.6726*10**(-27)

B0 = 800.0
n0 = 1.0*10**(27)
Lz = 80.0*10**(-6)


def dPhi_dt(Er):
    vA = B0/np.sqrt(mu*n0*mp)
    print(vA)
    return Lz*Er*B0*vA/1000

def Er(dPhi_dt):
    vA = B0/np.sqrt(mu*n0*mp)
    return dPhi_dt/(Lz*B0*vA)*1000


if __name__ == "__main__":
    print(dPhi_dt(0.2))
    print(Er(6.0))
