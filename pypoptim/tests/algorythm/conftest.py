import numpy as np
import pytest

from pypoptim.algorythm.solution import Solution


@pytest.fixture()
def square_solution():
    class SquareSolution(Solution):
        def update(self):
            self._y = np.sum(self.x ** 2)

        def is_valid(self):
            return self.is_updated()

    return SquareSolution


@pytest.fixture()
def maxabs_solution():
    class MaxAbsSolution(Solution):
        def update(self):
            self._y = np.max(np.abs(self.x))

        def is_valid(self):
            return self.is_updated()

    return MaxAbsSolution
