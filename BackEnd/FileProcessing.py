# calling modules
from BackEnd.Objects import SDFLoader
from BackEnd.DatabaseManager import DatabaseManager
from BackEnd.Utils import find_directory

# calling libraries
import tkinter as tk
from tkinter import filedialog, messagebox
import os


def load_file() -> str:
    """
    Open a file dialog to select a file and return the selected file path.

    :return: The selected file path.
    :rtype: str
    """
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    file_path = filedialog.askopenfilename()

    return file_path


def loading_sdf(file_path: str) -> str:
    """
    Load data from an SDF file and create a database with the extracted information.

    :param file_path: The path to the SDF file to load.
    :type file_path: str

    :return: The path to the created database.
    :rtype: str
    """
    # Finding the base directory "DRUG_DISCOVERY"
    find_directory("DRUG_DISCOVERY")
    # Create the directory
    new_dir = "Databases"
    os.makedirs(new_dir, exist_ok=True)

    # database path
    db_name = "drug_database"
    db_path = os.path.join(new_dir, db_name + '.db')

    ## Remove the database file if it exists
    if os.path.exists(db_path):
        print(f"Removing existing database at {db_path}")
        os.remove(db_path)

    ## load the sdf file

    if file_path:
        sdf_file = file_path
        molecule_list = []
        sdf_info = SDFLoader(sdf_file)

        for molecule in sdf_info.get_data():
            molecule_info = [molecule.CdId, molecule.image, molecule.smiles, molecule.IUPAC_name, molecule.MW,
                             molecule.HBA, molecule.HBD, molecule.LogP, molecule.LogD, molecule.PSA,
                             molecule.Rings, molecule.FusedAromaticRings, molecule.RotatableBonds]
            molecule_list.append(molecule_info)
        # creating the database manager
        db_manager = DatabaseManager(db_path)
        # creating database
        db_manager.create_db()
        # insert data into the database
        db_manager.insert_db(molecule_list)

        # display message
        messagebox.showinfo("Database Created", f"Database created successfully.")

    return db_path


#######################

def sort_by_column(list_of_lists: list, index: int) -> list:
    """
    Sorts a list of lists based on the specified index within each sublist.

    :param list_of_lists: The list of lists to be sorted.
    :type list_of_lists: list

    :param index: The index based on which the sorting should be done within each sublist.
    :type index: int

    :return: The sorted list of lists.
    :rtype: list
    """
    return sorted(list_of_lists, key=lambda x: x[index])


def sort_by_direction(sorted_data: list, descending: bool = False) -> list:
    """
    Rearranges sorted data in either ascending or descending order.

    :param sorted_data: The sorted data to be rearranged.
    :type sorted_data: list

    :param descending: If True, sorts the data in descending order, otherwise in ascending order. Default is False.
    :type descending: bool

    :return: The sorted data in the specified order.
    :rtype: list
    """
    if descending:
        return sorted_data[::-1]
    else:
        return sorted_data
