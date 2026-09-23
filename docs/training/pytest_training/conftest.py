import pytest
import requests

# conftext.py is a special configuration file for pytest that allows you to 
# define fixtures, hooks, and other configurations that can be shared across
# multiple test files in the same directory or subdirectories.
# THis is useful for setting up common test data, mocking external dependencies,
# or configuring test behavior in a centralized way.
# For example in your database project, you might have a conftest.py file that
# defines fixtures for creating test database connections, populating test data,
# or mocking external API calls. This way, you can reuse these fixtures across
# multiple test files without duplicating code.
# here an exampl of marking could be running all test that use a database
# (that you configure with your conftest.py file) with pytest -m database. 
# This will run all tests that are marked with the database marker, allowing
# you to selectively run tests that require a database connection.

@pytest.fixture(autouse=True)
def disable_network_calls(monkeypatch):
    def stunted_get():
        raise RuntimeError("Network access not allowed during testing!")
    monkeypatch.setattr(requests, "get", lambda *args, **kwargs: stunted_get())