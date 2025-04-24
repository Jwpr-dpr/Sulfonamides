from ase.atoms import Atoms
from ase.build import molecule
from ase.visualize import view
from gpaw import GPAW, PW, FermiDirac
from ase.io import write, read
from ase.optimize import BFGS
from ase.units import Bohr
from gpaw import GPAW, Mixer
from dftd4.ase import DFTD4
from ase.calculators.mixing import SumCalculator
import gpaw.solvation as solv
import os
import gc

def run_calculation(cif_dir="relax-d4-lcao-2", log_dir="logs", relax_dir="relax-d4-lcao-spinPol-solva"):
    """
    Realiza cálculos de relajación estructural para una serie de moléculas.

    Args:
        cif_dir (str): Directorio donde se encuentran los archivos XYZ de entrada.
        log_dir (str): Directorio para guardar los archivos de registro de GPAW.
        relax_dir (str): Directorio para guardar los archivos XYZ relajados.
    """

    # Crea los directorios si no existen
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(relax_dir, exist_ok=True)

    # Define las etiquetas de las moléculas
    qds = ['BPHQ', 'BPSQ', 'BPTQ']
    positions = ['P51', 'P52', 'P53', 'P61', 'P62', 'P63']
    # Itera sobre todas las combinaciones de etiquetas
    for qd in qds:
        for pos in positions:
            label = qd+'-SMZ-'+pos+'-relaxed-lcao' 
            cif_file = os.path.join(cif_dir, f"{label}.xyz")

            # Verifica si el archivo XYZ existe
            if not os.path.exists(cif_file):
               print(f"Warning: CIF file not found: {cif_file}")
               continue

            # Carga la estructura desde el archivo XYZ
            st = read(cif_file)

            # Centra la estructura y establece condiciones de contorno no periódicas
            st.center(vacuum=10)
            st.pbc = [False, False, False]

            # Define los nombres de los archivos de registro y salida
            log_file = os.path.join(log_dir, f"{label}-lcao-spinPol-relax-solva.log")
            xyz_file = os.path.join(relax_dir, f"{label}-lcao-spinPol.xyz")

            # Configura el cálculo GPAW
            calc = solv.SolvationGPAW(
                mode='lcao',  # Diferencias finitas
                xc='PBE',  # Funcional de intercambio-correlación PBE
                h=0.2,  # Espaciado de la malla
                basis='dzp',  # Base de doble zeta polarizada
                nbands ='110%',
	        spinpol = True,
                hund = True,
                parallel={'sl_auto': True, 'augment_grids': True,'use_elpa':False},  # Paralelización automática
                kpts={'size': (1, 1, 1), 'gamma': True},  # Punto gamma
                occupations=FermiDirac(0.01),  # Ocupaciones de Fermi-Dirac (temperatura electrónica)
                txt=log_file,  # Archivo de registro
                **solv.get_HW14_water_kwargs()
                )

            # Combina GPAW con DFTD4 para la corrección de dispersión
            d4 = SumCalculator([DFTD4(method='PBE'), calc])
            st.calc = d4

            # Realiza la relajación estructural usando el optimizador BFGS
            relax = BFGS(st)
            relax.run(fmax=0.04)  # Fuerza máxima permitida

            # Guarda la estructura relajada en un archivo XYZ
            st.write(xyz_file)

            # Limpia los objetos y libera memoria
            del calc
            gc.collect()

if __name__ == "__main__":
    run_calculation(cif_dir="relax-d4-lcao-2", log_dir="logs", relax_dir="relax-d4-lcao-spinPol-solva")
