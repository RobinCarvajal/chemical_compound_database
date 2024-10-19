# Imports
import os
import tkinter as tk
import webbrowser
# BackEnd
from BackEnd.DatabaseManager import DatabaseManager
from BackEnd.FileProcessing import load_file, loading_sdf
from BackEnd.Utils import find_directory


class MenuBar(tk.Menu):
    def __init__(self, master, visualiser):
        """
        Initialize a MenuBar instance.

        :param master: The parent widget.
        :param visualiser: The visualiser instance.
        """
        super().__init__(master)

        self.visualiser = visualiser
        self.data_table = visualiser.data_table

        # File menu
        file_menu = tk.Menu(self, tearoff=False)
        file_menu.add_command(label="Open SDF", command=self.open_file)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.quit)
        self.add_cascade(label="File", menu=file_menu)

        # Help menu
        help_menu = tk.Menu(self, tearoff=False)
        help_menu.add_command(label="Manual", command=self.display_manual)
        self.add_cascade(label="Help", menu=help_menu)

    def open_file(self):
        """
        Open an SDF file and display its content in the data table.
        """
        file_path = load_file()
        db_path = loading_sdf(file_path)
        print("Opening file...")
        db = DatabaseManager(db_path)
        headings = db.header_db()
        data = db.query_db()

        self.data_table.draw_table(headings, data)

    @staticmethod
    def display_manual():
        """
        Open the manual in a web browser.
        """
        find_directory("DRUG_DISCOVERY")
        webbrowser.open("file://" + os.getcwd() + "/Documentation/documentation.html", new=2)
        print("Displaying manual...")

