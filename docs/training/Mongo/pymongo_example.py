import pprint

from pymongo import MongoClient

# Establish a connection
client = MongoClient("localhost", 27017)
# URI format too: client = MongoClient("mongodb://localhost:27017")

# Access the database
db = client.rptutorials
## Note: You can also use the following syntax to access a database:
# db = client["rptutorials"]
# This dictionary-style syntax is useful if the database name contains special characters

# Access a collection
tutorial = db.tutorial

# Dictionaries to create documents
# Insert a single document
tutorial1 = {
    "title": "Working With JSON Data in Python",
    "author": "Lucas",
    "contributors": ["Aldren", "Dan", "Joanna"],
    "url": "https://realpython.com/python-json/",
}
result = tutorial.insert_one(tutorial1)
print(f"One tutorial: {result.inserted_id}")
# If only one document is insrted you can use
# tutoria1 = db.tutorial
# this way you spefcify which collection you want to use
# then you can use the insert_one method to insert a single document into that collection.
# result = tutorial.insert_one(tutorial1)
# print(f"One tutorial: {result.inserted_id}")

# Insert several documents
tutorial2 = {
    "title": "Python’s Requests Library (Guide)",
    "author": "Alex",
    "contributors": ["Aldren", "Brad", "Joanna"],
    "url": "https://realpython.com/python-requests/",
}

tutorial3 = {
    "title": "Object-Oriented Programming (OOP) in Python 3",
    "author": "David",
    "contributors": ["Aldren", "Joanna", "Jacob"],
    "url": "https://realpython.com/python3-object-oriented-programming/",
}

new_result = tutorial.insert_many([tutorial2, tutorial3])
print(f"Multiple tutorials: {new_result.inserted_ids}")

# Retrieve all documents (loop)
for doc in tutorial.find():
    pprint.pprint(doc)

# Retrieve a single document
jon_tutorial = tutorial.find_one({"author": "Jon"})
pprint.pprint(jon_tutorial)

# Close a connection protocol
with MongoClient() as client:
    db = client.rptutorials
    for doc in db.tutorial.find():
        pprint.pprint(doc)
