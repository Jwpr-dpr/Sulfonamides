from ase.atoms import Atoms
from ase.build import molecule
from ase.visualize import view
from gpaw import GPAW, PW, FermiDirac
from ase.io import write, read
from ase.optimize import BFGS
from ase.units import Bohr
from gpaw import GPAW,Mixer
from ase.calculators.dftd3 import DFTD3

cosas = ['BPHQ', 'BPSQ', 'BPTQ', 'SMZ']

for coso in cosas:
    label = coso
    st = read(label+'.cif')
    st.center(vacuum=10)
    st.pbc =[False,False,False]
    st.write(label+'.xyz')
