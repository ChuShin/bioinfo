import pytest
from he_caller.src.call_he import find_seeds

def test_find_seeds():
    results = find_seeds([])
    assert len(results) == 0