from ase.io import read
from gpaw import GPAW, FermiDirac, Mixer
from dftd4.ase import DFTD4
from ase.calculators.mixing import SumCalculator
from ase.parallel import paropen
import copy
import gc


def calculate_energies(input_dir='xyz-relax/', output_file='results_ads_eners_spinPol.csv'):
    """
    Calcula las energías de adsorción para diferentes configuraciones.

    Args:
        input_dir (str): Directorio donde se encuentran los archivos XYZ relajados.
        output_file (str): Nombre del archivo CSV de salida.
    """

    calc_dic = {'mode': 'lcao',
                'xc': 'PBE',
                'h': 0.2,
                'basis': 'dzp',
                'spinpol' : True,
                'hund' : True,
                'mixer': Mixer(0.07, 5, 100),
                'parallel': {'sl_auto': True, 'augment_grids': True},
                'kpts': {'size': (1, 1, 1), 'gamma': True},
                'occupations': FermiDirac(0.01)}

    with paropen(output_file, "w") as fout:
        fout.write("label,Etotal-0,Etotal-1,Emol-0,Emol-1,Esurf-0,Esurf-1\n")

        for punto in ['BPSQ', 'BPHQ', 'BPTQ','SMZ']:
 
            label = f"{punto}-relaxed-lcao"
            st = read(f'{input_dir}{label}.xyz')
            #st.center(vacuum=10)
            st.pbc = [False, False, False]
            calc_gpaw = GPAW(**calc_dic,
                    txt=f'logs/{label}-etotal-mol.log')
            st.calc = DFTD4(method='PBE').add_calculator(calc_gpaw)
            st.get_potential_energy()
            Etotal_ = [c.get_potential_energy() for c in st.calc.calcs]

            # Escritura de resultados
            fout.write(f"{label},{Etotal_[0]},{Etotal_[1]}\n")
            # Limpieza de memoria
            del calc_gpaw
            del st
            gc.collect()
if __name__ == "__main__":
    calculate_energies()

