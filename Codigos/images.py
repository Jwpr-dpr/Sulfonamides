'''
Rutina para correrlo en bash
El propósito de hacerlo de esta forma es poder evitar la duplicación de archivos en el cluster o en el pc que se trabaje, en vez de ello, accedemos directamente 
a la carpeta que contenga los archivos y realizamos nuestra tarea. Esto por medio de la definición de parámetros con ayuda de la librería fire. En este caso, 
los parámetros se pasan como argumentos al script en la línea de comandos, similar a una visión de C/C++, haciéndolos más fácilmente reutilizables.
Args:
    mol: string, nombre del punto cuántico (BPSQ, BPTQ o BPHQ)
    cif_path: string, ruta al directorio con los CIF originales (starting point)
    cif_relax_path: string, ruta al directorio con los CIFs de la estructuras relajadas (ending point)
    output_path: string, ruta al directorio donde se almacenará la imagen
Example:
    python images.py --mol="BPTQ" --cif_path="\Users\Dell\Dropbox\Sulfodamides\Resultados-actualizados\cif" --cif_relax_path="\Users\Dell\Dropbox\Sulfodamides\results\relaxed-DFTD3" --output_path="\Users\Dell\Dropbox\Sulfodamides"
'''
import os
import sys
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.visualize import view
from ase.io.pov import get_bondpairs, set_high_bondorder_pairs
import fire 


def render_structure(input_file, output_image,rotation=None, canvas_width=1200, high_bondorder_pairs=None):
    """
    Render a CIF file to an image with enhanced quality, bonds, and white background.
    """
    input_file = os.path.normpath(input_file)
    output_image = os.path.normpath(output_image)
    atoms = read(input_file)
    atoms.cell = None  # Removes the cell information
    pov_file = os.path.normpath(input_file.replace(".cif", ".pov"))
    # Export the structure to POV-Ray format
    atoms.write(pov_file)

    radii = [{'O': 1.0, 'S': 0.8, 'P': 0.6, 'C': 0.4, 'N': 0.3, 'H': 0.2}[at.symbol] for at in atoms]

    # Get bond pairs and set high bond orders if provided
    bondpairs = get_bondpairs(atoms, radius=1.1)
    if high_bondorder_pairs:
        bondpairs = set_high_bondorder_pairs(bondpairs, high_bondorder_pairs)

    # Write POV-Ray file with enhanced settings
    write(pov_file, atoms, format='pov',
          radii=radii,
          povray_settings=dict(canvas_width=canvas_width, bondatoms=bondpairs))


    # Enhance the POV-Ray rendering file by adding a white background and bonds
    with open(pov_file, 'a') as f:
        f.write('\nbackground { color rgb <1, 1, 1> }\n')  # White background
        f.write('\n#declare Bond_radius = 0.1;\n')  # Adjust bond radius if necessary
        f.write('\n#declare Bond_color = rgb <0.5, 0.5, 0.5>;\n')  # Bond color

    # Render the image using POV-Ray with high resolution
    os.system(f"pvengine /RENDER {pov_file} +W1600 +H1200 +O{output_image} -D -V /EXIT")


def main(mol:str, cif_path:str, cif_relax_path:str,output_path:str):
    """
    Main routine to create and save the molecules grid plot.
    """
    pos_list = ['P51', 'P52', 'P53', 'P61', 'P62', 'P63']
    original_images = []
    relaxed_images = []

    for pos in pos_list:
        # Define filenames
        original_file = os.path.normpath(f"{cif_path}/{mol}-SMZ-{pos}.cif")
        relaxed_file = os.path.normpath(f"{cif_relax_path}/{mol}-SMZ-{pos}-relaxed-dftd3.cif")
        
        print(f"Original file: {original_file}")
        print(f"Relaxed file: {relaxed_file}")

        # Render images
        original_image = original_file.replace(".cif", ".png")
        relaxed_image = relaxed_file.replace(".cif", ".png")
        
        render_structure(original_file, original_image)
        render_structure(relaxed_file, relaxed_image)
        
        original_images.append(original_image)
        relaxed_images.append(relaxed_image)
    
    
    fig, axes = plt.subplots(2, 6, figsize=(18, 6))
    for i, (orig_img, relax_img) in enumerate(zip(original_images, relaxed_images)):
        orig_img_data = plt.imread(orig_img)
        relax_img_data = plt.imread(relax_img)

        axes[0, i].imshow(orig_img_data)
        axes[0, i].axis('off')
        axes[0, i].set_title(f"{pos_list[i]} Starting Point")

        axes[1, i].imshow(relax_img_data)
        axes[1, i].axis('off')
        axes[1, i].set_title(f"{pos_list[i]} Ending Point")
    
    # Save the final plot
    plt.tight_layout()
    plt.savefig(f'{output_path}/{mol}Arrayplotdfdt.png')
    print(f"Grid plot saved at {output_path}")

if __name__ == "__main__":
    fire.Fire(main)