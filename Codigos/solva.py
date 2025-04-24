import gpaw.solvation as solv
from ase.io import read
from gpaw import GPAW, FermiDirac, Mixer
from dftd4.ase import DFTD4
from ase.calculators.mixing import SumCalculator
from ase.parallel import paropen
import copy

def get_ghost_setups(indice, target='ghost'):
    """
    Crea un diccionario de configuraciones para átomos fantasmas (ghost) o reales (paw).

    Args:
        indice (int): Índice límite para la configuración.
        target (str): 'ghost' para átomos fantasmas, 'paw' para átomos reales.

    Returns:
        dict: Diccionario de configuraciones.
    """
    my_setup = {'default': 'ghost' if target == 'ghost' else 'paw'}
    for i in range(indice+1):
        my_setup[i] = 'paw' if target == 'ghost' else 'ghost'
    return my_setup

def calculate_energies(input_dir='relax-d4-lcao-2/', output_file='results_ads_eners_spinPol-solva.csv'):
    """
    Calcula las energías de adsorción para diferentes configuraciones.

    Args:
        input_dir (str): Directorio donde se encuentran los archivos XYZ relajados.
        output_file (str): Nombre del archivo CSV de salida.
    """

    positions = ['P51', 'P52', 'P53', 'P61', 'P62', 'P63']
    puntos = {'BPSQ': 131, 'BPHQ': 118, 'BPTQ': 123}  # Diccionario para mapear puntos a índices

    calc_dic = {'mode': 'lcao',
                'xc': 'PBE',
                'h': 0.2,
                'basis': 'dzp',
                'spinpol' : True,
		'hund' : True,
                'parallel': {'sl_auto': True, 'augment_grids': True},
                'kpts': {'size': (1, 1, 1), 'gamma': True},
                'occupations': FermiDirac(0.01)}

    with paropen(output_file, "w") as fout:
        fout.write("label,Etotal-0,Etotal-1,Emol-0,Emol-1,Esurf-0,Esurf-1\n")

        for punto, indice_ in puntos.items():
            for pos in positions:
                label = f"{punto}-SMZ-{pos}"
                st = read(f'{input_dir}{label}-relaxed-lcao-lcao-spinPol.xyz')
                st.center(vacuum=10)
                st.pbc = [False, False, False]

                # Cálculo de la energía de la molécula
                st_mol = st.copy()
                calc_mol = solv.SolvationGPAW(**calc_dic,
                                setups=get_ghost_setups(indice_, target='ghost'),
                                txt=f'logs/{label}-etotal-mol-solva.log',
                                **solv.get_HW14_water_kwargs())
                st_mol.calc = DFTD4(method='PBE').add_calculator(calc_mol)
                st_mol.get_potential_energy()
                Etotal_mol = [c.get_potential_energy() for c in st_mol.calc.calcs]

                # Cálculo de la energía de la superficie
                st_surf = st.copy()
                calc_surf = GPAW(**calc_dic,
                                 setups=get_ghost_setups(indice_, target='paw'),
                                 txt=f'logs/{label}-etotal-surf-solva.log')
                st_surf.calc = DFTD4(method='PBE').add_calculator(calc_surf)
                st_surf.get_potential_energy()
                Etotal_surf = [c.get_potential_energy() for c in st_surf.calc.calcs]

                # Cálculo de la energía total
                st_total = st.copy()
                calc_total = GPAW(**calc_dic,
                                  txt=f'logs/{label}-etotal-solva.log')
                st_total.calc = DFTD4(method='PBE').add_calculator(calc_total)
                st_total.get_potential_energy()
                Etotal = [c.get_potential_energy() for c in st_total.calc.calcs]

                # Escritura de resultados
                fout.write(f"{label},{Etotal[0]},{Etotal[1]},{Etotal_mol[0]},{Etotal_mol[1]},{Etotal_surf[0]},{Etotal_surf[1]}\n")

                # Limpieza de memoria
                del calc_mol, calc_surf, calc_total, st_mol, st_surf, st_total

if __name__ == "__main__":
    calculate_energies(input_dir="relax-d4-lcao-spinPol/")
