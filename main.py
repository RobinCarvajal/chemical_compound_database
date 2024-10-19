# Imports
import os
import tkinter as tk
from GUI_components.App import App


########################################################################################################################
#                                                      Running the App                                                 #
########################################################################################################################
def main():
    """
    Main function to run the drug discovery application.

    It creates the application window, initializes the database path, and removes any existing database file.

    :return: None
    :rtype: None
    """
    root = tk.Tk()
    app = App(root, "Drug Discovery", (1200, 600))
    root.mainloop()
    db_path = "Databases/drug_database.db"
    if os.path.exists(db_path):
        print(f"Removing existing database at {db_path}")
        os.remove(db_path)


if __name__ == "__main__":
    main()
