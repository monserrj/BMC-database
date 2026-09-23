import pytest
from format_data import format_data_for_display

@pytest.fixture
def example_people_data():
    return [
        {
            "given_name": "Alfonsa",
            "family_name": "Ruiz",
            "title": "Senior Software Engineer",
        },
        {
            "given_name": "Sayid",
            "family_name": "Khan",
            "title": "Project Manager",
        },
    ]

# ...
def test_format_data_for_display(example_people_data):
    # people = [
    #     {
    #         "given_name": "Alfonsa",
    #         "family_name": "Ruiz",
    #         "title": "Senior Software Engineer",
    #     },
    #     {
    #         "given_name": "Sayid",
    #         "family_name": "Khan",
    #         "title": "Project Manager",
    #     },
    # ]
# Test the expected output of the format_data_for_display function with the 
# given input data. The expected output is a list of formatted strings, each
# representing a person in the input list.
    assert format_data_for_display(example_people_data) == [
        "Alfonsa Ruiz: Senior Software Engineer",
        "Sayid Khan: Project Manager",
    ]

# Useful pytest plugins:

# pytest-randomly: Randomizes the order of test execution to help identify inter-test dependencies.
# You can use the seed value to reproduce test in the same order to identify the issue.

# pytest-cov: Measures code coverage of your tests and generates reports to help you identify untested code paths.

# pytest-django: Provides additional functionality for testing Django applications, including database fixtures and test client support.

# ppytest-bdd: Provides support for behavior-driven development (BDD) testing using the Gherkin syntax.

# pytest-mock: Provides a simple way to create mock objects for testing, allowing you to isolate and test specific parts of your code.
# For example, you can use pytest-mock to mock external API calls or database queries, allowing you to test your code without relying on external dependencies.

# pytest-xdist: Allows you to run tests in parallel across multiple CPUs or machines, which can significantly speed up test execution time for large test suites.

# pytest-timeout: Allows you to set timeouts for individual tests or test suites, which can help prevent tests from hanging indefinitely and improve test reliability.

# pytest-html: Generates HTML reports for your test results, which can be useful for sharing test results with stakeholders or for tracking test progress over time.

# pytest-sugar: Provides a more visually appealing output format for test results, including progress bars and color-coded output.

# pytest-faker: Provides a simple way to generate fake data for testing, which can be useful for testing edge cases or for generating large amounts of test data quickly.

# pytest-metadata: Allows you to add metadata to your test results, such as the test environment or test configuration, which can be useful for tracking test results across different environments or configurations.

# pytest-rerunfailures: Allows you to automatically rerun failed tests a specified number of times, which can be useful for dealing with flaky tests or intermittent failures.

# Summary of key commands for pytest:
# pytest: Run all tests in the current directory and its subdirectories.
# pytest <test_file.py>: Run tests in a specific test file.
# pytest <test_file.py>::<test_function>: Run a specific test function in a test file.
# pytest -k <expression>: Run tests that match a specific expression (e.g., pytest -k "test_uppercase" will run only the test_uppercase function).
# pytest -m <marker>: Run tests that are marked with a specific marker (e.g., pytest -m "slow" will run all tests marked with the "slow" marker).
# pytest --maxfail=<num>: Stop after a specified number of test failures (e.g., pytest --maxfail=3 will stop after 3 test failures).
# pytest --tb=<style>: Set the traceback style (e.g., pytest --tb=short will show a shorter traceback).
# pytest -v: Increase verbosity of test output (e.g., pytest -v will show more detailed output for each test).
# pytest --disable-warnings: Disable warning messages in test output.