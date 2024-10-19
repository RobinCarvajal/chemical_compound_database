# Imports
import io
from tkinter import ttk
from PIL import Image, ImageTk


class DataTable(ttk.Treeview):
    def __init__(self, master):
        """
        Initialize a DataTable instance.

        :param master: The parent widget.
        """
        super().__init__(master)

        # Create a Scrollbar for vertical scrolling
        y_scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.yview)
        y_scrollbar.pack(side="right", fill="y")
        self.configure(yscrollcommand=y_scrollbar.set)

        # Create a Scrollbar for horizontal scrolling
        x_scrollbar = ttk.Scrollbar(self, orient="horizontal", command=self.xview)
        x_scrollbar.pack(side="bottom", fill="x")
        self.configure(xscrollcommand=x_scrollbar.set)

        self.stored_headings = list()
        self.stored_data = list()
        self.stored_photo_images = list()

        # Configure style to set row height
        style = ttk.Style()
        style.configure("Custom.Treeview", rowheight=250)
        self.configure(style="Custom.Treeview")

    def draw_table(self, headings, data):
        """
        Draw the table with the given headings and data.

        :param headings: List containing column headings.
        :param data: List of lists containing data to be displayed.
        """
        # Set title for column 0 as "Image" and its width to 270
        self.heading("#0", text="Image")
        self.column("#0", width=260, minwidth=260, stretch=False)

        self.delete(*self.get_children())
        headings.pop(1)
        columns = headings
        self["columns"] = columns
        ## HEADINGS
        for col in columns:
            self.column(col, anchor="center", stretch=False)
            self.heading(col, text=col, anchor="center")

        ## ROWS
        # Create PhotoImage objects from image bytes

        self.stored_photo_images = []
        for row in data:
            row = list(row)
            image_bytes = row[1]
            # Convert bytes to an image object
            image = Image.open(io.BytesIO(image_bytes))
            # Resize image to 200x200 pixels
            image = image.resize((200, 200))
            # Create Tkinter PhotoImage object
            photo = ImageTk.PhotoImage(image)

            self.stored_photo_images.append(photo)  # Append the PhotoImage object to the list

        # Insert data into rows
        for index, row in enumerate(data):
            key = ""

            modified_values = list(row)
            modified_values.pop(1)  # Remove the item at index 1 (image)

            self.insert("", "end",
                        image=self.stored_photo_images[index],
                        text="", values=modified_values)

        return None
