# chemical_compound_database

This shows the code of a graphical user interphase to explore and visualise a chemical compounds database. Build using base python and tkinter.
## Methods

### Modularity
The program is built modularly for both the back and front end.

### Data extraction and parsing
Data from the .sdf file of bioactive compounds was extracted using RDKit OSS in Python. The file contained a hundred compounds, some with predefined properties and others computed by RDKit. Data extraction was done via a custom generator class (`Objects.SDFLoader`), which created `Molecule` data class (`Objects.Molecule`) objects. These attributes were stored in lists for database insertion.

### Database Management
The `DatabaseManager` module sets up an SQLite database ("drug_database.db"). The database will always be created at this location and is erased when the program closes. It handles data insertion and querying using Python's SQLite3 package.

### Images generation
Images of the molecules were obtained using `rdkit.Chem.Draw`, transformed into bytes using `io.BytesIO`, and stored in the database as BLOB. Later, they were decoded and used in the GUI.

### Fused Aromatic Rings Count
The `count_fused_aromatic_rings` function uses RDKit to identify fused aromatic rings in a molecule by analyzing its smallest set of smallest rings (SSSR). It checks for aromaticity and detects fused rings by examining shared atoms. The count is determined from the unique fused ring indices.

### Filtering criteria
The database consistently retrieves uniform information, accessing all data from the main table with each query. However, filtering occurs upon data retrieval within the program, utilizing functions aligned with widely recognized criteria such as Lipinski's rules, lead-likeness, and bioavailability. This approach ensures that only relevant data is processed and displayed in the GUI.

## Properties Table

| Property                           | Lipinski’s | Lead-likeness | Bioavailability |
|------------------------------------|------------|---------------|-----------------|
| Molecular Weight                   | ≤ 500      | ≤ 450         | ≤ 500           |
| LogD                               | -          | -4 ≤ x ≤ 4    | ≤ 5             |
| Number of Hydrogen-bond Donors     | ≤ 5       | ≤ 5           | ≤ 5             |
| Number of Hydrogen-bond Acceptors   | ≤ 10      | ≤ 8           | ≤ 10            |
| LogP                               | ≤ 5       | -             | -               |
| Ring Count                         | -          | ≤ 4           | -               |
| Fused Aromatic Rings Count        | -          | -             | ≤ 5             |
| Polar Surface Area                 | -          | -             | ≤ 200           |
| Rotatable Bonds Count              | -          | ≤ 10          | -               |
| Lipinski’s                         | ≤ 10      | -             | ≤ 10            |

## Front-End Design
The application's front-end, built with *Tkinter*, employs a modular structure with distinct classes for various GUI elements. **MenuBar** manages the menu, **Tabs** handles tab navigation, and **Visualiser** houses **DataTable** for database visualisation. **DataTable** decodes byte data into images using `Image` and `ImageTK`. **App** merges these components for a cohesive interface, enhancing code organisation and usability. Each class has methods that enable data visualisation and interaction with other classes when running the GUI.

### Starting GUI
The GUI program is triggered by executing the `main.py` script in the command line in this fashion: `python3 main.py` in the *DRUG_DATABASE* directory.

## How to use the GUI
The database loading process starts with selecting an SDF file through the File menu's OpenSDF option, prompting the user to navigate directories and choose the file. Subsequently, the program initiates background processing to create the database, which is confirmed by a pop-up upon successful creation. Exiting the application is also accessible via the File menu.

Initially, the database appears in the visualiser frame, which can be navigated through horizontal and vertical scroll bars, presenting unfiltered and unsorted data. The Sort tab streamlines database management with Filter, Sort by, Order, and a Search box. Filters include Lipinski’s rule, Lead-likeness, and Bioavailability parameters. Sorting by numerical values is facilitated by "Sort by", applicable only to numerical columns, while "Order" permits ascending or descending arrangement. The “Search” box targets the molecule column.

The Export tab facilitates exporting the current table in CSV or XLS format, initiating a directory selection window upon activation. Additionally, the Help menu accesses local HTML documentation within the DRUG_DATABASE directory, offering guidance without internet dependency. The database is deleted upon program closure to mitigate storage issues and ensure future database creation without data conflicts.
