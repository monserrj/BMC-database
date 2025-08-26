#!/usr/bin/env python

''' This script extracts and prints all docstrings from Python files in a specified project folder.'''
import ast
from pathlib import Path

def extract_docstrings(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        tree = ast.parse(f.read(), filename=file_path)

    docs = []

    # module-level docstring
    mod_doc = ast.get_docstring(tree)
    if mod_doc:
        docs.append(("module", mod_doc))

    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            name = node.name
            docstring = ast.get_docstring(node)
            if docstring:
                docs.append((name, docstring))
    return docs


def extract_from_project(project_dir):
    project_dir = Path(project_dir)
    all_docs = {}

    for file in project_dir.rglob("*.py"):
        all_docs[file.name] = extract_docstrings(file)

    return all_docs


# Example usage
# Replace "your_project_folder" with your code folder path
docs = extract_from_project(".")

for fname, entries in docs.items():
    print(f"\n📄 {fname}")
    for name, doc in entries:
        print(f"\n--- {name} ---\n{doc}\n")
