def Lipinski(data_from_query: list) -> list:
    """
    Apply Lipinski's Rule of Five to filter molecules.

    :param data_from_query: A list of tuples representing molecule data fetched from the database.
    :type data_from_query: list

    :return: A list of tuples containing molecule data that pass Lipinski's Rule of Five criteria.
    :rtype: list
    """
    # Empty list to store molecules passing Lipinski's Rule of Five
    filtered_molecules = []

    # Constants representing the cutoff values for Lipinski's Rule of Five
    MW_cutoff = 500
    HBA_cutoff = 10
    HBD_cutoff = 5
    LogP_cutoff = 5

    # Iterate over each row of molecule data
    for row in data_from_query:
        # Unpack the row into individual variables for easier comparison
        ID, Image, Smiles, Molecule, MW, HBA, HBD, LogP, LogD, PSA, Rings, FusedAromaticRings, NoRotatableBonds = row

        # Check if the molecule passes Lipinski's Rule of Five
        if (MW <= MW_cutoff and
                HBA <= HBA_cutoff and
                HBD <= HBD_cutoff and
                LogP <= LogP_cutoff):
            # If the molecule passes, add it to the list of filtered molecules
            filtered_molecules.append(row)

    # Return the list of molecules passing Lipinski's Rule of Five
    return filtered_molecules


def LeadLikeness(data_from_query: list) -> list:
    """
    Apply lead-likeness filters to molecules.

    :param data_from_query: A list of tuples representing molecule data fetched from the database.
    :type data_from_query: list

    :return: A list of tuples containing molecule data that pass lead-likeness filters.
    :rtype: list
    """
    # Empty list to store molecules passing lead-likeness filters
    filtered_molecules = []

    # Constants representing the cutoff values for lead-likeness filters
    MW_cutoff = 450
    LogD_cutoff = 4
    Rings_cutoff = 4
    NoRotatableBonds_cutoff = 10
    HBA_cutoff = 8
    HBD_cutoff = 5

    # Iterate over each row of molecule data
    for row in data_from_query:
        # Unpack the row into individual variables for easier comparison
        ID, Image, Smiles, Molecule, MW, HBA, HBD, LogP, LogD, PSA, Rings, FusedAromaticRings, NoRotatableBonds = row

        # Check if the molecule passes lead-likeness filters
        if (MW <= MW_cutoff and
                (-LogD_cutoff <= LogD <= LogD_cutoff) and
                Rings <= Rings_cutoff and
                NoRotatableBonds <= NoRotatableBonds_cutoff and
                HBA <= HBA_cutoff and HBD <= HBD_cutoff):
            # If the molecule passes, add it to the list of filtered molecules
            filtered_molecules.append(row)

    # Return the list of molecules passing lead-likeness filters
    return filtered_molecules


def Bioavailability(data_from_query: list) -> list:
    """
    Apply bioavailability filters to molecules.

    :param data_from_query: A list of tuples representing molecule data fetched from the database.
    :type data_from_query: list

    :return: A list of tuples containing molecule data that pass bioavailability filters.
    :rtype: list
    """
    # Empty list to store molecules passing bioavailability filters
    filtered_molecules = []

    # Constants representing the cutoff values for bioavailability filters
    MW_cutoff = 500
    LogP_cutoff = 5
    HBD_cutoff = 5
    HBA_cutoff = 10
    NoRotatableBonds_cutoff = 10
    PSA_cutoff = 200
    FusedAromaticRings_cutoff = 5

    # Iterate over each row of molecule data
    for row in data_from_query:
        # Unpack the row into individual variables for easier comparison
        ID, Image, Smiles, Molecule, MW, HBA, HBD, LogP, LogD, PSA, Rings, FusedAromaticRings, NoRotatableBonds = row

        # Counter for parameters met
        parameters_met = 0

        # Check if parameters meet cutoff values
        if MW <= MW_cutoff:
            parameters_met += 1
        if LogP <= LogP_cutoff:
            parameters_met += 1
        if HBD <= HBD_cutoff:
            parameters_met += 1
        if HBA <= HBA_cutoff:
            parameters_met += 1
        if NoRotatableBonds <= NoRotatableBonds_cutoff:
            parameters_met += 1
        if PSA <= PSA_cutoff:
            parameters_met += 1
        if FusedAromaticRings <= FusedAromaticRings_cutoff:
            parameters_met += 1

        # Check if at least 6 parameters are met
        if parameters_met >= 6:
            # If the molecule passes, add it to the list of filtered molecules
            filtered_molecules.append(row)

    # Return the list of molecules passing bioavailability filters
    return filtered_molecules
