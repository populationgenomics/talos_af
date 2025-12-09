import pytest


def test_something():
    with pytest.raises(TypeError):
        raise TypeError('pass')
