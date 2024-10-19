# Imports
import tkinter as tk
from tkinter import ttk
# BackEnd
from BackEnd.DatabaseManager import DatabaseManager
from BackEnd.FileProcessing import sort_by_column, sort_by_direction
from BackEnd import Filters
from BackEnd.Utils import lists_to_csv, lists_to_excel


class Tabs(ttk.Frame):
    def __init__(self, master, visualiser):
        """
        Initialize the Tabs frame.

        :param master: The parent widget.
        :param visualiser: An instance of the visualiser.
        """
        super().__init__(master)

        self.columnconfigure(0, weight=1)
        # to display the data in treeview
        self.visualiser = visualiser
        self.data_table = visualiser.data_table

        ## For filters and sorting
        # variables
        self.sort_by_var = None
        self.reverse_sort_var = None
        self.filter_var = None
        # menus
        self.filter_menu = None
        self.sort_by_menu = None
        self.reverse_sort_menu = None
        # searchbox
        self.search_entry = None

        # Create widgets
        self.create_widgets()

    def create_widgets(self):
        """
        Create notebook(tabs) containing sorting and export tabs.

        """
        # Create notebook(tabs)
        notebook = ttk.Notebook(self)
        notebook.grid(column=0, row=0, sticky='nsew')

        # Sort Tab
        sort_tab = ttk.Frame(self)
        notebook.add(sort_tab, text="Sort")
        # Call the method to create sorting dropdown menu
        self.create_sorting_widget(sort_tab)

        # Export Tab
        export_tab = ttk.Frame(notebook)
        notebook.add(export_tab, text="Export")
        self.create_export_widget(export_tab)

    def create_sorting_widget(self, master):
        """
        Create sorting widget containing comboboxes and an entry box for sorting, filtering, and searching.

        :param master: The parent widget.
        """
        # List of headings for the sorting options
        headings = ["ID", "MW", "HBA", "HBD", "LogP", "LogD", "PSA", "Rings", "FusedAromaticRings", "RotatableBonds"]

        # String variables for storing selected values
        sort_by_var = tk.StringVar()
        sort_by_var.set("ID")

        reverse_sort_var = tk.StringVar()
        reverse_sort_var.set("ASC")

        filter_var = tk.StringVar()
        filter_var.set("No filter")

        # Create a frame to contain all sorting widgets
        sort_frame = ttk.Frame(master)
        sort_frame.grid(row=0, column=0, padx=10, pady=30, sticky="nsew")

        # Label and combobox for filtering
        filter_label = ttk.Label(sort_frame, text="Filter:")
        filter_label.grid(row=0, column=0, padx=(0, 5))

        self.filter_menu = ttk.Combobox(sort_frame, textvariable=filter_var,
                                        values=('No filter', 'Lipinski', 'Lead-Likeness', 'Bioavailability'))
        self.filter_menu.grid(row=0, column=1, padx=(0, 5))

        # Label and combobox for sorting by column
        sort_by_label = ttk.Label(sort_frame, text="Sort by:")
        sort_by_label.grid(row=0, column=2, padx=(0, 5))

        self.sort_by_menu = ttk.Combobox(sort_frame, textvariable=sort_by_var, values=headings)
        self.sort_by_menu.grid(row=0, column=3, padx=(0, 5))

        # Label and combobox for sorting order
        reverse_sort_label = ttk.Label(sort_frame, text="Order:")
        reverse_sort_label.grid(row=0, column=4, padx=(0, 5))

        self.reverse_sort_menu = ttk.Combobox(sort_frame, textvariable=reverse_sort_var, values=('ASC', 'DESC'))
        self.reverse_sort_menu.grid(row=0, column=5, padx=(0, 5))

        # Label and entry box for searching
        search_label = ttk.Label(sort_frame, text="Search:")
        search_label.grid(row=0, column=6, padx=(10, 5))

        self.search_entry = ttk.Entry(sort_frame, width=40)  # Adjust width here
        self.search_entry.grid(row=0, column=7, padx=(0, 5))
        self.search_entry.bind("<Return>", self.searchbox_filter)  # Bind Return event

        # Store variables for later use
        self.sort_by_var = sort_by_var
        self.reverse_sort_var = reverse_sort_var
        self.filter_var = filter_var

        # Bind combobox selection to sorting function
        self.filter_menu.bind("<<ComboboxSelected>>", self.click_filter)
        self.sort_by_menu.bind("<<ComboboxSelected>>", self.click_sort_column)
        self.reverse_sort_menu.bind("<<ComboboxSelected>>", self.click_sort_direction)

    def create_export_widget(self, master):
        """
        Create export widget containing buttons to export current table data as Excel or CSV.

        :param master: The parent widget.
        """
        reverse_sort_label = ttk.Label(master, text="Export current table as:")
        reverse_sort_label.grid(row=0, column=0, padx=10, pady=30)

        # Export button
        export_excel_button = ttk.Button(master, text="Excel", command=self.export_as_xls)
        export_excel_button.grid(column=1, row=0, padx=10, pady=30)

        # Export button
        export_csv_button = ttk.Button(master, text="CSV", command=self.export_as_csv)
        export_csv_button.grid(row=0, column=2, padx=10, pady=30)

    def click_filter(self, event=None):
        """
        Apply the selected filter to the data and redraw the table.

        :param event: The event that triggers the filter (default is None).
        :return: Filtered data based on the selected filter.
        """
        # Get the selected filter value
        selected_filter = self.filter_var.get()

        # Perform filtering based on the selected filter value
        if selected_filter == "No filter":
            # Display the original table
            unfiltered_data = self.no_filter()
            print("Not filtered...")
            return unfiltered_data
        elif selected_filter == "Lipinski":
            filtered_data = self.lipinski_filter()
            print("Filtering by Lipinski...")
            return filtered_data
        elif selected_filter == "Lead-Likeness":
            filtered_data = self.lead_likeness_filter()
            print("Filtering by Lead-Likeness...")
            return filtered_data
        elif selected_filter == "Bioavailability":
            filtered_data = self.bioavailability_filter()
            print("Filtering by Bioavailability...")
            return filtered_data

    def click_sort_column(self, event=None):
        """
        Sort the data based on the selected column and redraw the table.

        :param event: The event that triggers the sort column (default is None).
        :return: Sorted data based on the selected column.
        """
        # getting the value in the combobox
        selected_filter = self.sort_by_var.get()
        # getting the headings
        filepath = "Databases/drug_database.db"
        db = DatabaseManager(filepath)
        headings = db.header_db()
        # Call click_filter to get the filtered table
        filtered_data = self.click_filter()

        # Perform actions based on the selected filter
        filter_to_column = {
            "ID": 0,
            "MW": 4,
            "HBA": 5,
            "HBD": 6,
            "LogP": 7,
            "LogD": 8,
            "PSA": 9,
            "Rings": 10,
            "FusedAromaticRings": 11,
            "RotatableBonds": 12
        }

        if selected_filter in filter_to_column:
            column_index = filter_to_column[selected_filter]
            sorted_by_column = sort_by_column(filtered_data, column_index)
            print(f"Sorting by {selected_filter}...")
            self.data_table.draw_table(headings, sorted_by_column)
            return sorted_by_column

    def click_sort_direction(self, event=None):
        """
        Sort the data based on the selected sort direction (ASC or DESC) and redraw the table.

        :param event: The event that triggers the sort direction (default is None).
        :return: Sorted data based on the selected sort direction.
        """
        # Get the selected sort direction
        sort_direction = self.reverse_sort_var.get()
        # getting the headings
        filepath = "Databases/drug_database.db"
        db = DatabaseManager(filepath)
        headings = db.header_db()
        # Call click_sort_column to get the sorted table
        sorted_by_column = self.click_sort_column()

        # Perform actions based on the selected sort direction (ASC or DESC)
        if sort_direction == "ASC":
            # Sorting logic for ascending order
            sorted_by_direction = sort_by_direction(sorted_by_column, descending=False)
            print("Sorting in ascending order...")
            self.data_table.draw_table(headings, sorted_by_direction)
            return sorted_by_direction
        elif sort_direction == "DESC":
            # Sorting logic for descending order
            sorted_by_direction = sort_by_direction(sorted_by_column, descending=True)
            print("Sorting in descending order...")
            self.data_table.draw_table(headings, sorted_by_direction)
            return sorted_by_direction

    def searchbox_filter(self, event=None):
        """
        Filter the data based on the search query entered in the search box and redraw the table.

        :param event: The event that triggers the searchbox filter (default is None).
        :return: Filtered data based on the search query.
        """
        # getting the headings
        filepath = "Databases/drug_database.db"
        db = DatabaseManager(filepath)
        headings = db.header_db()

        # Getting the search query
        search_query = str(self.search_entry.get().strip())

        # Getting the data
        data = self.click_sort_direction()

        # Filtered data
        filtered_data = []

        # Only search in the "molecule" column (assuming it's the first column)
        for row in data:
            if search_query.lower() in row[3].lower():
                filtered_data.append(row)

        # Update the data table
        self.data_table.draw_table(headings, filtered_data)

        return filtered_data


    def export_as_csv(self):
        """
        Export the current table data as a CSV file.

        """
        # getting the headings
        filepath = "Databases/drug_database.db"
        db = DatabaseManager(filepath)
        headings = db.header_db()
        headings.pop(1)  # remove the image column
        data = self.searchbox_filter()
        # remove the image column
        new_data = []
        for sublist in data:
            sublist = list(sublist)
            del sublist[1]
            new_data.append(sublist)
        # exporting as csv
        lists_to_csv(headings, new_data)

    def export_as_xls(self):
        """
        Export the current table data as an Excel file.

        """
        # getting the headings
        filepath = "Databases/drug_database.db"
        db = DatabaseManager(filepath)
        headings = db.header_db()
        headings.pop(1)  # remove the image column
        data = self.searchbox_filter()
        # remove the image column
        new_data = []
        for sublist in data:
            sublist = list(sublist)
            del sublist[1]
            new_data.append(sublist)
        # exporting as xls
        lists_to_excel(headings, new_data)

    def no_filter(self):
        """
        Display the data without applying any filters and redraw the table.

        :return: Unfiltered data.
        """
        filepath = "Databases/drug_database.db"
        db = DatabaseManager(filepath)
        headings = db.header_db()
        data = db.query_db()
        self.data_table.draw_table(headings, data)
        return data

    def lipinski_filter(self):
        """
        Filter the data by Lipinski criteria and redraw the table.

        :return: Filtered data after applying the Lipinski filter.
        """
        filepath = "Databases/drug_database.db"
        db = DatabaseManager(filepath)
        headings = db.header_db()
        data = db.query_db()
        filtered_data = Filters.Lipinski(data)
        self.data_table.draw_table(headings, filtered_data)
        return filtered_data

    def lead_likeness_filter(self):
        """
        Filter the data by lead-likeness and redraw the table.

        :return: Filtered data after applying the lead-likeness filter.
        """
        filepath = "Databases/drug_database.db"
        db = DatabaseManager(filepath)
        headings = db.header_db()
        data = db.query_db()
        filtered_data = Filters.LeadLikeness(data)
        self.data_table.draw_table(headings, filtered_data)
        return filtered_data

    def bioavailability_filter(self):
        """
        Filter the data by bioavailability and redraw the table.

        :return: Filtered data after applying the bioavailability filter.
        """
        filepath = "Databases/drug_database.db"
        db = DatabaseManager(filepath)
        headings = db.header_db()
        data = db.query_db()
        filtered_data = Filters.Bioavailability(data)
        self.data_table.draw_table(headings, filtered_data)
        return filtered_data
