#!/usr/bin/env python

''' This script defines the restricted vocabulary of the database for Enum values. Instructions
followed and explanations kept to help following through

'''
from enum import Enum, unique
# All Enum values as they are constant should be capitalised for good practice
# https://realpython.com/python-enum/

# ---- Helper function to convert string to Enum ----

def enum_from_str(enum_class: type[Enum], label: str) -> Enum:
    """
    Convert a string label into an Enum member of the given enum_class.
    Case-insensitive and trims spaces.
    
    Example:
        enum_from_str(DatabaseType, "seq") -> DatabaseType.SEQUENCE
    """
    label = label.strip().upper()
    try:
        return enum_class[label]
    except KeyError:
        valid = [e.name for e in enum_class]
        raise ValueError(
            f"Invalid value {label!r} for {enum_class.__name__}. "
            f"Valid options: {valid}"
        )

# ---- EnumS ----
@unique
class DatabaseType(Enum):
    """Enum for classifying types of external biological databases."""
    STRUCTURE = "STRUCTURE"
    SEQUENCE = "SEQUENCE"
    FUNCTION = "FUNCTION"
    TAXONOMY = "TAXONOMY"
    # OTHER = "OTHER" Unsure if needed?



# ---- Testing ----

if __name__ == "__main__":
    print("Testing DatabaseType Enum...\n")

    # List all members for testing this works. Remove after
    for db_type in DatabaseType:
        print(f"{db_type.name} = {db_type.value}")

    print("\nTesting helper function:")
    print("\nDatabaseType tests:")
    for label in ["seq", "STRUCTURE", " tax ", "invalid"]:
        try:
            result = enum_from_str(DatabaseType, label)
            print(f"{label!r} -> {result}")
        except ValueError as e:
            print(f"{label!r} -> Error: {e}")