import sqlite3 # Import library of your database
from numpy.random import uniform #import random generator

# EXAMPLE 1
connection = sqlite3.connect ("survey.db")  # Establish a connection to database
cursor = connection.cursor()                # Create a cursor to keep track: A pointer into a database that keeps track of outstanding operations.
cursor.execute("SELECT Site.lat, Site.long FROM Site;") # Use cursor to ask database to execute query
results = cursor.fetchall()                 # Return the result of our query on list form
for r in results:                           # Loop over list
    print(r)                                # Print results
cursor.close()
connection.close()

# EXAMPLE 2
def get_name(database_file, person_id): # Create a function
    query = "SELECT personal || ' ' || family FROM Person WHERE id=?;" #Use the question mark to provide a safe list of our query

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query, [person_id])
    results = cursor.fetchall()
    cursor.close()
    connection.close()

    return results[0][0]

print("Full name for dyer:", get_name('survey.db', 'dyer')) #Where you indicate the name to extract

#EXAMPLE 3:
def add_name(database_file, new_person):
    query = "INSERT INTO Person (id, personal, family) VALUES (?, ?, ?);"

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query, list(new_person))
    cursor.close()
    connection.commit() # When you add something you need to commit before closing to save
    connection.close()
# Insert a new name
add_name('survey.db', ('barrett', 'Mary', 'Barrett'))
# Check it exists
print("Full name for barrett:", get_name('survey.db', 'barrett'))

# EXAMPLE 4:
random_numbers = uniform(low=10.0, high=25.0, size=100000) #  To get random numbers

connection = sqlite3.connect ("original.db")  # Create database
cursor = connection.cursor()                # Create a cursor to keep track: A pointer into a database that keeps track of outstanding operations.
cursor.execute("CREATE TABLE Pressure (reading float not null);") # Create table
query = "INSERT INTO PRESSURE (reading) VALUES (?);" #To create the query to insert those numbers

for n in random_numbers:                           # Loop over list
    cursor.execute(query,[n])                                # Insert random numbers

cursor.close()
connection.commit() #Save changes
connection.close()

# EXAMPLE 5 way 1:
connection_original = sqlite3.connect ("original.db")  # Search database
cursor_original = connection_original.cursor()                # Create a cursor to keep track
cursor_original.execute("SELECT * FROM Pressure;") # Select all values table
results = cursor_original.fetchall()        # Call valuesa and save as results
cursor_original.close()
connection_original.close()

connection_backup = sqlite3.connect("backup.db")
cursor_backup = connection_backup.cursor()
cursor_backup.execute("CREATE TABLE Pressure (reading float not null)")
query = "INSERT INTO Pressure (reading) VALUES (?);" #To create the query to insert those numbers

for entry in results: #Call those results back
    if entry [0] > 20.0: #Is a list call 1 by 1
        cursor_backup.execute (query, entry)

cursor_backup.close()
connection_backup.commit() #Save changes
connection_backup.close()

#EXAMPLE 5 way 2:
connection_original = sqlite3.connect("original.db")
cursor_original = connection_original.cursor()
cursor_original.execute("SELECT * FROM Pressure WHERE reading > 20.0;") #ADD Conditioning here
results = cursor_original.fetchall()
cursor_original.close()
connection_original.close()

connection_backup = sqlite3.connect("backup.db")
cursor_backup = connection_backup.cursor()
cursor_backup.execute("CREATE TABLE Pressure (reading float not null)")
query = "INSERT INTO Pressure (reading) VALUES (?);"

for entry in results:
    cursor_backup.execute(query, entry)

cursor_backup.close()
connection_backup.commit()
connection_backup.close()
