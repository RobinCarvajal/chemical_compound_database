################################################### IMPORTS ############################################################
from io import BytesIO  # func: get_molecule_image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# calling libraries
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import os


########################################################################################################################
#                                                UTILITIES                                                             #
########################################################################################################################

###################################### FUNCTION TO CREATE AN IMAGE FROM MOL ############################################

def get_molecule_image(mol):
    """
    Explanation: Generates an image_byte_array of a molecule
    :param mol: (rdkit.Chem.Mol)- The molecule.
    :return image_byte_array: The image in binary format.
    """
    # Generate 2D coordinates
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol)

    # Convert image to binary data
    img_byte_array = BytesIO()
    img.save(img_byte_array, format='PNG')
    img_byte_array = img_byte_array.getvalue()

    return img_byte_array


######################################## FUSED AROMATIC RINGS COUNT FUNCTION ###########################################
def count_fused_aromatic_rings(mol):
    """
    Count the number of fused aromatic rings in a molecule.

    :param mol: The molecule.
    :type mol: rdkit.Chem.Mol

    :return: Number of fused aromatic rings.
    :rtype: int
    """
    # Get the SSSR (smallest set of smallest rings)
    sssr = Chem.GetSymmSSSR(mol)

    # Initialize a set to store unique fused aromatic rings
    fused_aromatic_ring_indices = set()

    # Check each ring in the SSSR
    for i, ring in enumerate(sssr):
        # Check if the ring is aromatic
        if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring):
            # Check if any atom in the ring is also part of another aromatic ring
            for j, other_ring in enumerate(sssr):
                if i != j and any(atom in other_ring for atom in ring) and all(
                        mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in other_ring):
                    # This ring is fused with another aromatic ring
                    fused_aromatic_ring_indices.add(i)

    return len(fused_aromatic_ring_indices)

## Testing function

# Define SMILES strings and expected counts
structures = {
    'c1ccccc1': 0,
    'C1=CC=C2C=CC=CC2=C1': 2,
    'C1=CC=C2C=C3C=CC=CC3=CC2=C1': 3,
    'C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2': 4,
    'C1=CC=C2C=C3C=C4C=C5C=CC=CC5=CC4=CC3=CC2=C1': 5,
    'C1=CC=C2C=C3C(=CC2=C1)C=CC4=CC5=CC6=CC=CC=C6C=C5C=C43': 6,
    'C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7': 7,
    'C1=CC2=C3C4=C1C=CC5=CC6=C7C8=C(C=CC9=C8C1=C(C=C9)C=C(C3=C1C7=C54)C=C2)C=C6': 10
}

# Perform assertions
for smiles, expected_count in structures.items():
    assert count_fused_aromatic_rings(Chem.MolFromSmiles(smiles)) == expected_count

############################################### NEW FUNCTION ###########################################################


def find_directory(directory_name):
    """
    Change the current working directory to the specified directory name by traversing up the directory hierarchy.

    :param directory_name: The name of the directory to search for and change to.
    :type directory_name: str

    :return: None
    :rtype: None

    Prints:
        If the directory is found, prints the path of the new current working directory.
        If the directory is not found, prints a message indicating that the directory was not found.
    """
    # Get the current working directory
    current_path = os.getcwd()

    # Continue looping until we find the desired directory or reach the root directory
    while os.path.basename(current_path) != directory_name:
        # Get the parent directory path
        parent_path = os.path.dirname(current_path)

        # If the parent directory is the same as the current directory,
        # it means we've reached the root directory and the desired directory doesn't exist
        if parent_path == current_path:
            print(f"Directory '{directory_name}' not found.")
            return

        # Move up to the parent directory
        current_path = parent_path

    # Change the current working directory to the desired directory
    os.chdir(current_path)

    # Print the path of the new current working directory
    print(f"Changed directory to: {os.getcwd()}")

########################################################################################################################

def ask_csv_filename():
    """
    Opens a file dialog to select a CSV file and returns the selected file path.

    :return: The selected file path.
    :rtype: str
    """
    root = tk.Tk()
    root.withdraw()  # Hide the main Tkinter window

    user_input = filedialog.asksaveasfilename(initialdir="/",  # Initial directory to open the dialog
                                              title="Select CSV File",  # Title of the dialog window
                                              filetypes=(("CSV files", "*.csv"),    # File types filter for CSV files
                                                         ("All files", "*.*")))     # File types filter for all files

    return user_input  # Return the path selected by the user


########################################################################################################################

def ask_excel_filename():
    """
    Opens a file dialog to select an Excel file and returns the selected file path.

    :return: The selected file path.
    :rtype: str
    """
    root = tk.Tk()
    root.withdraw()  # Hide the main Tkinter window

    user_input = filedialog.asksaveasfilename(initialdir="/",  # Initial directory to open the dialog
                                              title="Select Excel File",  # Title of the dialog window
                                              filetypes=(("Excel files", "*.xls"),  # File types filter for Excel files
                                                         ("All files", "*.*")))     # File types filter for all files

    return user_input  # Return the path selected by the user


########################################################################################################################

def lists_to_csv(headings, data_lists):
    """
    Convert lists of data into a CSV file.

    :param headings: List of column headings.
    :type headings: list

    :param data_lists: List of lists containing the data.
    :type data_lists: list[list]

    :return: None
    :rtype: None
    """
    file_path = ask_csv_filename()
    # Create DataFrame
    df = pd.DataFrame(data_lists, columns=headings)

    # Export DataFrame to CSV
    try:
        df.to_csv(file_path, index=False)
        messagebox.showinfo("Export Successfully", f"CSV exported successfully. CSV does not contain images, they can be generated using the SMILES.")

    except Exception as e:
        messagebox.showinfo("Export Unsuccessful", f"An error occurred while exporting the DataFrame: {str(e)}")

    return None

########################################################################################################################

def lists_to_excel(headings, data_lists):
    """
    Export data_lists to an Excel file with the specified headings.

    :param headings: List containing column headings.
    :param data_lists: List of lists containing data to be exported.
    """
    file_path = ask_excel_filename()

    if not file_path:
        return

    # Create DataFrame
    df = pd.DataFrame(data_lists, columns=headings)

    # Export DataFrame to Excel
    try:
        # Writing to .xls format using xlwt engine
        df.to_excel(file_path, index=False, engine='xlsxwriter')
        messagebox.showinfo("Export Successfully", "Excel exported successfully.")

    except Exception as e:
        messagebox.showinfo("Export Unsuccessful", f"An error occurred while exporting the DataFrame: {str(e)}")
