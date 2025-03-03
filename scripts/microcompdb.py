#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db import create_db
from data import load_data


def main():
    create_db()
    load_data()


# The syntax here is "special"
# Every file in Python has a __name__ attribute
# If the file is being run as the main program, __name__ is set to "__main__"
# This happens if the file is run from the command line as
# python scripts/monolith.py
# If the file is being imported, __name__ is set to the name of the file
# This is a common Python idiom
if __name__ == "__main__":
    main()
