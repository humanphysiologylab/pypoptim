import itertools

import numpy as np
import pytest

from ...algorythm.solution import Solution


class TestSolution:
    def test_init(self):

        invalid_xs = [42, "foo", [], [[0, 1], [2, 3]]]
        for x in invalid_xs:
            with pytest.raises(ValueError):
                Solution(x)

        sol = Solution([42])
        assert not sol.is_updated()
        with pytest.raises(NotImplementedError):
            sol.is_valid()
        with pytest.raises(NotImplementedError):
            sol.update()

    def test_len(self):
        n = 42
        assert len(Solution(np.random.random(n))) == n

    def test_data(self):

        x = [0]

        sol = Solution(x)
        assert "smth" not in sol
        with pytest.raises(KeyError):
            smth = sol["smth"]

        kw = dict(a=1, b="foo")
        sol = Solution(x, **kw)
        for key, value in kw.items():
            assert sol[key] == value
            assert key in sol

        forty_two = 42
        sol["c"] = forty_two
        assert sol._data["c"] == forty_two

    def test_solution_update(self, square_solution):

        x = [1, 2]
        sol = square_solution(x)
        assert len(sol) == len(x)
        assert not sol.is_updated()
        assert not sol.is_valid()

        sol.update()
        assert sol.is_updated()
        assert sol.is_valid()
        assert sol.y == sum(xi ** 2 for xi in x)

        sol.x = x
        assert not sol.is_updated()
        assert not sol.is_valid()

    def test_comparators(self, maxabs_solution):
        def zip_product(xs, sols):
            return zip(itertools.product(xs, xs), itertools.product(sols, sols))

        xs = [1, 2, 2, 3]
        sols = [maxabs_solution([x]) for x in xs]

        for _, sol_pair in zip_product(xs, sols):
            sol_1, sol_2 = sol_pair
            with pytest.raises(ValueError):
                sol_1 == sol_2
            with pytest.raises(ValueError):
                sol_1 > sol_2

        for sol in sols:
            sol.update()

        for x_pair, sol_pair in zip_product(xs, sols):
            x_1, x_2 = x_pair
            sol_1, sol_2 = sol_pair
            assert (x_1 == x_2) == (sol_1 == sol_2)
            assert (x_1 != x_2) == (sol_1 != sol_2)
            assert (x_1 > x_2) == (sol_1 > sol_2)
            assert (x_1 < x_2) == (sol_1 < sol_2)
            assert (x_1 >= x_2) == (sol_1 >= sol_2)
            assert (x_1 <= x_2) == (sol_1 <= sol_2)
