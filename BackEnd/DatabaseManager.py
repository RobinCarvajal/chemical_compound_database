#####################################################################################################################
#                                          DATABASE OBJECTS AND METHODS                                             #
#####################################################################################################################
import sqlite3


class DatabaseManager:
    def __init__(self, db_path):
        """
        Initialize a DatabaseManager instance.

        :param db_path: The path to the SQLite database file.
        """
        self._db_path = db_path

    # getter
    @property
    def db_path(self):
        """
        Getter for the database path property.

        :return: The path to the SQLite database file.
        """
        return self._db_path

    def create_db(self):
        """
        Create the drug_database table if it does not exist.
        """
        connection = sqlite3.connect(self.db_path)
        cur = connection.cursor()
        sql = """
            CREATE TABLE IF NOT EXISTS chemicals (
                ID INT PRIMARY KEY,
                Image BLOB,
                Smiles TEXT,
                Molecule TEXT,
                MW REAL,
                HBA INT, 
                HBD INT,
                LogP REAL, 
                LogD REAL, 
                PSA REAL,
                Rings INT,
                FusedAromaticRings INT, 
                RotatableBonds INT
            )
            """
        cur.execute(sql)
        cur.close()
        connection.close()

    def insert_db(self, params):
        """
        Insert data into the chemicals table in the database.

        :param params: List of tuples containing data to be inserted.
        """
        sql = '''
            INSERT INTO chemicals (ID, Image, Smiles, Molecule, MW, HBA, HBD, LogP, LogD, PSA, Rings, FusedAromaticRings, RotatableBonds)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        '''
        error = False

        connection = sqlite3.connect(self.db_path)
        cur = connection.cursor()

        try:
            cur.executemany(sql, params)

        except Exception as e:
            print(f'An exception has been encountered: {e}')
            error = True
            connection.rollback()

        if not error:
            connection.commit()

        cur.close()
        connection.close()

    def query_db(self):
        """
        Fetch all data from the chemicals table in the database.

        :return: List of tuples containing query results.
        """
        sql = "SELECT * FROM chemicals"

        connection = sqlite3.connect(self.db_path)
        cur = connection.cursor()
        cur.execute(sql)
        results = cur.fetchall()
        cur.close()
        connection.close()

        return results

    def header_db(self):
        """
        Retrieve the header of the chemicals table.

        :return: List of column names.
        """
        connection = sqlite3.connect(self.db_path)
        cur = connection.cursor()
        cur.execute(f"PRAGMA table_info(chemicals)")
        columns_info = cur.fetchall()
        column_names = [column[1] for column in columns_info]
        cur.close()
        connection.close()

        return column_names


