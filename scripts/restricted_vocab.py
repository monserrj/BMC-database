#!/usr/bin/env python

''' This script defines the restricted vocabulary of the database for Enum values. Instructions
followed and explanations kept to help following through

'''

from enum import Enum, auto, unique

# All Enum values as they are constant should be capitalised for good practice
# https://realpython.com/python-enum/
# Another option is to create empty enums and assign values later
# Unsure which is better practice

@unique # Values assigned must be unique otherwise ValueError is raised
class DatabaseType(Enum):
    '''Enum for the different database types'''
    STRUCT = "STRUCTURE"
    SEQ = "SEQUENCE"
    FUNC = "FUNCTION"
    TAX = "TAXONOMY"
    # OTHER = "Other" Unsure if needed?

    '''or'''
    STRUCT, SEQ, FUNC, TAX = auto()

    # OTHER = 5 Unsure if needed?

'''OR'''
Enum(
    value,
    names,
    *,
    module=None,
    qualname=None,
    type=None,
    start=1
)
# Example:
HTTPMethod = Enum(
...     "HTTPMethod", ["GET", "POST", "PUSH", "PATCH", "DELETE"]
... )

# Exampkle 2:
>>> HTTPStatusCode = Enum(
...     value="HTTPStatusCode",
...     names=[
...         ("OK", 200),
...         ("CREATED", 201),
...         ("BAD_REQUEST", 400),
...         ("NOT_FOUND", 404),
...         ("SERVER_ERROR", 500),
...     ],
... )
# I think this would be the better ption for the DatabaseType Enum
# THis is an example for now

>>> names = []
>>> while True:
...     name = input("Member name: ")
...     if name in {"q", "Q"}:
...         break
...     names.append(name.upper())
...
Member name: YES
Member name: NO
Member name: q

>>> DynamicEnum = Enum("DynamicEnum", names)
>>> list(DynamicEnum)
[<DynamicEnum.YES: 1>, <DynamicEnum.NO: 2>]

# Note to self: Use database type to sort protdxref vs genexref issue?
# Why I though of that:
>>> class Semaphore(Enum):
...     RED = 1
...     YELLOW = 2
...     GREEN = 3
...

>>> def handle_semaphore(light):
...     match light:
...         case Semaphore.RED:
...             print("You must stop!")
...         case Semaphore.YELLOW:
...             print("Light will change to red, be careful!")
...         case Semaphore.GREEN:
...             print("You can continue!")
...

# Do I want Enum or IntEnum?

#Brief note on IntFlag in case is handy for something I am not yet aware
>>> from enum import IntFlag

>>> class Role(IntFlag):
...     OWNER = 8
...     POWER_USER = 4
...     USER = 2
...     SUPERVISOR = 1
...     ADMIN = OWNER | POWER_USER | USER | SUPERVISOR
