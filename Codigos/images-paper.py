import os
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.io.pov import get_bondpairs, set_high_bondorder_pairs
import fire
import tempfile
import shutil

def render_structure(input_file, output_image, rotation=None, canvas_width=1200, high_bondorder_pairs=None):
    """
    Renderiza un archivo CIF a una imagen usando POV-Ray en un directorio temporal seguro.

    Parámetros:
        input_file (str): Ruta al archivo CIF.
        output_image (str): Ruta de salida para la imagen PNG.
        rotation (str): Eje para rotar la vista desde (ej., 'y' para vista lateral, 'x', 'z').
        canvas_width (int): Ancho del lienzo renderizado.
        high_bondorder_pairs (set): Pares de átomos con enlaces dobles/triples.
    """
    atoms = read(input_file)
    atoms.cell = None
    atoms.pbc = False

    # Crear directorio temporal
    with tempfile.TemporaryDirectory() as tmpdir:
        pov_file = os.path.join(tmpdir, "temp.pov")
        temp_output = os.path.join(tmpdir, "temp.png")

        radii = [{'O': 1.0, 'S': 0.8, 'P': 0.6, 'C': 0.4, 'N': 0.3, 'H': 0.2}[at.symbol] for at in atoms]
        bondpairs = get_bondpairs(atoms, radius=1.1)
        if high_bondorder_pairs:
            bondpairs = set_high_bondorder_pairs(bondpairs, high_bondorder_pairs)

        write(pov_file, atoms, format='pov',
              radii=radii,
              povray_settings=dict(canvas_width=canvas_width, bondatoms=bondpairs))

        with open(pov_file, 'a') as f:
            f.write('\nbackground { color rgb <1, 1, 1> }\n')
            f.write('\n#declare Bond_radius = 0.1;\n')
            f.write('\n#declare Bond_color = rgb <0.5, 0.5, 0.5>;\n')

            # Añadir rotación de la cámara con el vector 'up' corregido
            if rotation == 'y': 
                f.write('\ncamera { location <0, -10, 35> look_at <0, 1, 0> up <0, 1, 0> rotate <90, 0, 0> }\n') #ajustar aquís
            elif rotation == 'z':
                f.write('\ncamera { location <0, 0, 25> look_at <0, 0, 0> up <0, 1, 0> }\n')
            else:
                f.write('\ncamera { location <0, 0, -10> look_at <0, 0, 0> up <0, 1, 0> }\n')

        # Cambiar directorio y renderizar
        old_cwd = os.getcwd()
        os.chdir(tmpdir)
        os.system(f"pvengine /RENDER temp.pov +W{canvas_width} +H{canvas_width} +Otemp.png -D -V /EXIT")
        os.chdir(old_cwd)

        # Copiar resultado a la ruta de salida
        shutil.copy(temp_output, output_image)


def main(cif_path: str, cif_relax_path: str, output_path: str):
    """
    Renderiza moléculas específicas desde dos planos (z e y) y genera un gráfico de cuadrícula.
    """
    molecules = [
        ("BPSQ", "P52"),
        ("BPHQ", "P52"),
        ("BPTQ", "P61")
    ]
    views = [("z", "z"), ("y", "y")]  # Usar "z" e "y" para las vistas

    original_images = []
    relaxed_images = []

    for mol, pos in molecules:
        for view_name, rotation in views:
            original_file = os.path.join(cif_path, f"{mol}-SMZ-{pos}.cif")
            relaxed_file = os.path.join(cif_relax_path, f"{mol}-SMZ-{pos}-relaxed-dftd3.cif")

            original_img = os.path.join(output_path, f"{mol}_{pos}_original_{view_name}.png")
            relaxed_img = os.path.join(output_path, f"{mol}_{pos}_relaxed_{view_name}.png")

            render_structure(original_file, original_img, rotation=rotation)
            render_structure(relaxed_file, relaxed_img, rotation=rotation)

            original_images.append(original_img)
            relaxed_images.append(relaxed_img)

    fig, axes = plt.subplots(2, 6, figsize=(18, 6))
    titles = [f"{mol}-{pos} {v}-view" for mol, pos in molecules for v, _ in views]

    for i, (orig_img, relax_img) in enumerate(zip(original_images, relaxed_images)):
        axes[0, i].imshow(plt.imread(orig_img))
        axes[0, i].axis('off')
        axes[0, i].set_title(f"{titles[i]}-startPoint")

        axes[1, i].imshow(plt.imread(relax_img))
        axes[1, i].axis('off')
        axes[1, i].set_title(f"{titles[i]}-endPoint")

    axes[1, 0].set_ylabel("Relajado")
    axes[0, 0].set_ylabel("Original")
    plt.tight_layout()
    output_plot = os.path.join(output_path, "grid_plot.png")
    plt.savefig(output_plot)
    print(f"Gráfico guardado: {output_plot}")

if __name__ == "__main__":
    fire.Fire(main)
