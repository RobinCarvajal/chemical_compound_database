# Imports
from tkinter import ttk
# BackEnd
from GUI_components.DataTable import DataTable

class Visualiser(ttk.Frame):
    def __init__(self, master):
        """
        Initialize the Visualiser frame.

        :param master: The parent widget.
        """
        super().__init__(master)

        # self.grid(column=0, row=1, sticky='nsew')
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        style = ttk.Style()
        style.configure("White.TFrame", background="white")
        self.configure(style="White.TFrame")

        # Treeview
        self.data_table = DataTable(self)
        self.data_table.grid(column=0, row=0, padx=10, pady=10, sticky='nsew')
