#!/usr/bin/env python

'''This script is for testing connection to bmc.db as I am getting eerros with the microcomp_db.y script
'''
import sys
from sqlalchemy import create_engine

engine = create_engine(
    "sqlite:///C:/Users/monserrj/BMC-database/scripts/test_simple.db"
)

with engine.connect() as conn:
    print("Connected OK")
