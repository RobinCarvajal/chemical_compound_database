# BackEnd
from GUI_components.MenuBar import MenuBar
from GUI_components.Tabs import Tabs
from GUI_components.Visualiser import Visualiser

class App:
    # Main setup
    def __init__(self, master, title, size):
        """
        Initialize the application.

        :param master: The master Tkinter window.
        :type master: tkinter.Tk
        :param title: The title of the application window.
        :type title: str
        :param size: The size of the application window in pixels (width, height).
        :type size: tuple
        """
        super().__init__()  # Initialize the superclass
        self.master = master  # Set the master window
        self.master.title(title)  # Set the title of the window
        self.master.geometry(f'{size[0]}x{size[1]}')  # Set the geometry of the window
        self.master.minsize(size[0], size[1])  # Set the minimum size of the window

        # Resizing configurations
        self.master.columnconfigure(0, weight=1)  # Allow column 0 to resize horizontally
        self.master.rowconfigure(1, weight=1)  # Allow row 1 to resize vertically

        # Widgets setup
        # Visualiser setup
        self.visualiser = Visualiser(self.master)  # Create the visualiser widget
        self.visualiser.grid(column=0, row=1, sticky='nsew')  # Place the visualiser in the grid layout

        # Menu bar setup
        self.menu_bar = MenuBar(self.master, self.visualiser)  # Create the menu bar
        self.master.config(menu=self.menu_bar)  # Configure the menu bar

        # Tabs setup
        self.tabs = Tabs(self.master, self.visualiser)  # Create the tabs widget
        self.tabs.grid(column=0, row=0, sticky='nsew')  # Place the tabs in the grid layout

        # Run the Tkinter main loop
        self.master.mainloop()
