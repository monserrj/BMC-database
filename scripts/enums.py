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

@unique
class StructProtType(Enum):
    """Enum for classifying types of structural proteins."""
    HEXAMER= "HEXAMER"
    TRIMER = "TRIMER"
    PENTAMER = "PENTAMER"
    OTHER = "OTHER"

@unique
class ModificationType(Enum):
    """Enum for classifying types of protein modifications."""
    TRUNCATED = "TRUNCATED"
    EXTENDED = "EXTENDED"
    FUSION = "FUSION"
    SYNTHETIC = "SYNTHETIC"
    MUTATED = "MUTATED"
    DOMESTICATED = "DOMESTICATED"


@unique
class ComplexSource(Enum):
    """Enum for classifying sources of protein complexes."""
    NATIVE = "NATIVE"
    ENGINEERED = "ENGINEERED"
    PREDICTED = "PREDICTED"
    THEORETICAL = "THEORETICAL" # same than predicted? need to describe differences

# ---- Testing ----

if __name__ == "__main__":
    print("Testing DatabaseType Enum...\n")

    # List all members for testing this works. Remove after
    for db_type in ModificationType:
        print(f"{db_type.name} = {db_type.value}")

    print("\nTesting helper function:")
    print("\nDatabaseType tests:")
    for label in ["BMC-H", "MUTATED", " fuse ", "FUSION", "invalid"]:
        try:
            result = enum_from_str(ModificationType, label)
            print(f"{label!r} -> {result}")
        except ValueError as e:
            print(f"{label!r} -> Error: {e}")