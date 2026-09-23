# If I include a function with test prefix this is a pytest
# pytest will run all the tests in this file or folder?
# Results:
# . test passed
# F test failed
# E unexpected exception

# A test must be simple and isolated. It should not depend on other tests or external factors. Each test should be able to run independently and produce the same
# result every time. Each test should be easy to understand and maintain. It should clearly state what it is testing and what the expected outcome is.


def test_always_passes():
    assert True

def test_always_fails():
    assert False

# ASSERT examples
# Assert is a statement that checks if a condition is true. If the condition is false, it raises an AssertionError exception.
# In pytest, assert statements are used to verify that the code behaves as expected.
def test_uppercase():
    assert "loud noises".upper() == "LOUD NOISES"

def test_reversed():
    assert list(reversed([1, 2, 3, 4])) == [4, 3, 2, 1]

def test_some_primes():
    assert 37 in {
        num
        for num in range(2, 50)
        if not any(num % div == 0 for div in range(2, num))
    }
    
# Fixtures
# Fixture are functions that run before each test function to which it is applied. They are used to set up some state or context for the tests.
# pytest fixtures are functions that can create data, test doubles, or initialize system state for the test suite.
# Any test that wants to use a fixture must explicitly use this fixture function as an argument to the test function, so dependencies are always stated up front
# Fixtures can also make use of other fixtures, again by declaring them explicitly as dependencies. That means that, over time, your fixtures can become bulky and modular.
import pytest

@pytest.fixture
def example_fixture():
    return 1

def test_with_fixture(example_fixture):
    assert example_fixture == 1

# Filtering tests by name-based filtering (-k) and marker-based filtering (-m)
# pytest -k "test_uppercase" will run only the test_uppercase function
# Marker-based filtering allows you to group tests by markers and run them selectively. You can mark a test with a custom marker using the @pytest.mark decorator
# and then run tests with that marker using the -m option. For example, you can mark a test with @pytest.mark.slow and then run all slow tests with pytest -m slow.
# Other ways of filtering are by directory scope, file scope, class scope, and module scope and test categories. You can also use the -v option to get more verbose
# output, which will show you which tests are being run and their results.

# Test parameterization
# Test parameterization allows you to run the same test function with different input values.
# You can use the @pytest.mark.parametrize decorator to specify the input values for the test function.
# For example, you can use @pytest.mark.parametrize("input,expected", [(1, 2), (2, 3), (3, 4)]) to run the test function with three different input values
# and their expected results.

def test_is_palindrome_empty_string():
    assert is_palindrome("")

def test_is_palindrome_single_character():
    assert is_palindrome("a")

def test_is_palindrome_mixed_casing():
    assert is_palindrome("Bob")

def test_is_palindrome_with_spaces():
    assert is_palindrome("Never odd or even")

def test_is_palindrome_with_punctuation():
    assert is_palindrome("Do geese see God?")

def test_is_palindrome_not_palindrome():
    assert not is_palindrome("abc")

def test_is_palindrome_not_quite():
    assert not is_palindrome("abab")

# Example above repeats a lot the same test in different ways.
# Instead we can use parameterization to run the same test function with different
# input values and expected results.
@pytest.mark.parametrize("palindrome", [
    "",
    "a",
    "Bob",
    "Never odd or even",
    "Do geese see God?",
])
def test_is_palindrome(palindrome):
    assert is_palindrome(palindrome)

@pytest.mark.parametrize("non_palindrome", [
    "abc",
    "abab",
])
def test_is_palindrome_not_palindrome(non_palindrome):
    assert not is_palindrome(non_palindrome)

# Further parametrization can be achieved by
@pytest.mark.parametrize("maybe_palindrome, expected_result", [
    ("", True),
    ("a", True),
    ("Bob", True),
    ("Never odd or even", True),
    ("Do geese see God?", True),
    ("abc", False),
    ("abab", False),
])
def test_is_palindrome(maybe_palindrome, expected_result):
    assert is_palindrome(maybe_palindrome) == expected_result

# Parametrization can be risky if use to much, as it can make the test suite
# harder to understand and maintain. It is important to strike a balance between
# parametrization and test clarity.