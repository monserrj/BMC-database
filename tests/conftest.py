import pytest

@pytest.fixture 
def cleanup():
    to_delete = []
    yield to_delete
    for item in to_delete:
        if item.is_file():
            item.unlink()
