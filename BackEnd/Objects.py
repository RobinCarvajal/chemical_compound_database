################################################### IMPORTS ############################################################

# packages
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import rdMolDescriptors
from dataclasses import dataclass
from tkinter import filedialog, messagebox, simpledialog

# modules
from BackEnd.Utils import count_fused_aromatic_rings, get_molecule_image


########################################################################################################################
#                                                   OBJECTS                                                            #
########################################################################################################################

@dataclass
class Molecule:
    CdId: str
    image: str
    smiles: str
    IUPAC_name: str
    MW: float
    HBA: int
    HBD: int
    LogP: float
    LogD: float
    PSA: float
    Rings: int
    FusedAromaticRings: int
    RotatableBonds: int

class SDFLoader:
    def __init__(self, FileName):
        """
        Initialize SDFLoader with the given file name.

        :param FileName: The name of the SDF file to load.
        """
        self._FileName = FileName

    @property
    def FileName(self):
        """
        Getter for the file name property.

        :return: The name of the SDF file.
        """
        return self._FileName

    def get_data(self) -> Molecule:
        """
        Generator function to parse molecules from the SDF file and yield Molecule objects.

        :return: A Molecule object representing each parsed molecule.
        """
        molecules_not_parsed = []  # To store names and SMILES of molecules that did not parse
        all_molecules_parsed = True  # Flag to track if all molecules parsed successfully

        # load the data from the sdf file
        with Chem.SDMolSupplier(self.FileName) as supplier:
            for mol in supplier:
                if mol:
                    try:
                        # Ro5 descriptors
                        smiles = Chem.MolToSmiles(mol)
                        IUPAC_name = mol.GetProp('Name')
                        CdId = mol.GetProp('CdId')
                        MW = Descriptors.MolWt(mol)
                        HBA = Descriptors.NOCount(mol)
                        HBD = Descriptors.NHOHCount(mol)
                        LogP = Descriptors.MolLogP(mol)
                        LogD = mol.GetProp('LogD')
                        Rings = rdMolDescriptors.CalcNumRings(mol)
                        FusedAromaticRings = count_fused_aromatic_rings(mol)
                        PSA = Descriptors.TPSA(mol)
                        RotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

                        # Generate structure image
                        # img_path = get_molecule_image(mol, CdId)
                        img_bytes = get_molecule_image(mol)

                        yield Molecule(CdId, img_bytes, smiles, IUPAC_name, MW, HBA, HBD, LogP, LogD, PSA, Rings,
                                       FusedAromaticRings, RotatableBonds)

                    except Exception as e:
                        # If an exception occurs, log the molecule name and SMILES
                        molecules_not_parsed.append((IUPAC_name, smiles))
                        all_molecules_parsed = False

        # Save names and SMILES of molecules that did not parse into a txt file
        with open("molecules_not_parsed.txt", "w") as file:
            for name, smiles in molecules_not_parsed:
                file.write(f"{name}\t{smiles}\n")

        # Display messagebox if any molecules could not be parsed
        if not molecules_not_parsed:
            messagebox.showinfo("All molecules parsed", "All molecules in the SDF file parsed correctly.")
        elif all_molecules_parsed:
            messagebox.showinfo("Errors while parsing",
                                "All molecules in the SDF file could not be parsed. Please check the input file.")
        else:
            messagebox.showinfo("Errors while parsing",
                                "Some molecules could not be parsed. Details are located in a .txt file.")
