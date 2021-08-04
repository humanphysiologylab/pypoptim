import numpy as np
import pytest

from ....algorythm.ga import GA
from ....algorythm.solution import Solution


class TestGA:
    def test_init_invalid(self):
        class SquareSolution(Solution):
            def update(self):
                self._y = np.sum(self.x ** 2)

            def is_valid(self):
                return self.is_updated()

        invalids = Solution, "smth"
        for invalid in invalids:
            with pytest.raises(TypeError):
                GA(SolutionSubclass=invalid, bounds=...)

        bounds_invalid = (
            [],
            [[], []],
            [0, 1],
            [[0, 1, 2], [0, 1, 2], [0, 1, 2]],
            [[1, 0], [0, 1]],
        )
        for bounds in bounds_invalid:
            with pytest.raises(ValueError):
                GA(SolutionSubclass=SquareSolution, bounds=bounds)

        bounds = [[-3, 5], [1, 10]]

        gammas_invalid = [1], [1, 1, 1], [-1, 1]
        for gammas in gammas_invalid:
            with pytest.raises(ValueError):
                GA(SolutionSubclass=SquareSolution, bounds=bounds, gammas=gammas)

        mask_log10_scale_invalids = (
            [False],
            [False, False, False],
            [True, False],
            [True, True],
        )
        for mask_log10_scale in mask_log10_scale_invalids:
            with pytest.raises(ValueError):
                GA(
                    SolutionSubclass=SquareSolution,
                    bounds=bounds,
                    mask_log10_scale=mask_log10_scale,
                )

        invalids = -1, 1.1
        for invalid in invalids:
            with pytest.raises(ValueError):
                GA(
                    SolutionSubclass=SquareSolution,
                    bounds=bounds,
                    mutation_rate=invalid,
                )
            with pytest.raises(ValueError):
                GA(
                    SolutionSubclass=SquareSolution,
                    bounds=bounds,
                    crossover_rate=invalid,
                )

        with pytest.raises(TypeError):
            GA(SolutionSubclass=SquareSolution, bounds=bounds, selection_force=2.0)
        with pytest.raises(ValueError):
            GA(SolutionSubclass=SquareSolution, bounds=bounds, selection_force=1)

        invalids = "abs", np.array(["foo", "bar"])
        for invalid in invalids:
            with pytest.raises(TypeError):
                GA(
                    SolutionSubclass=SquareSolution,
                    bounds=bounds,
                    keys_data_transmit=invalid,
                )

        invalids = "rng", 42
        for invalid in invalids:
            with pytest.raises(TypeError):
                GA(SolutionSubclass=SquareSolution, bounds=bounds, rng=invalid)

    def test_is_inside_bounds(self, ga_optim_default):
        ga_optim = ga_optim_default
        bounds = ga_optim.bounds

        sol = ga_optim_default.generate_solution()
        sol.x = [0, 3]
        assert ga_optim.is_solution_inside_bounds(sol, bounds)
        assert ga_optim.is_solution_inside_bounds(sol)

        sol.x = [0, 0]
        assert not ga_optim.is_solution_inside_bounds(sol, bounds)
        assert not ga_optim.is_solution_inside_bounds(sol)

    def test_generate_solution(self, ga_optim_fabric):

        ga_optim_1 = ga_optim_fabric()
        sol_1 = ga_optim_1.generate_solution()

        ga_optim_2 = ga_optim_fabric()
        sol_2 = ga_optim_2.generate_solution()
        sol_2_another = ga_optim_2.generate_solution()

        assert np.all(sol_1.x == sol_2.x)
        assert np.all(sol_2.x != sol_2_another.x)  # super likely
        assert not sol_1.is_updated()
        assert not sol_1.is_valid()

        for sol in sol_1, sol_2, sol_2_another:
            assert ga_optim_1.is_solution_inside_bounds(sol)

    def test_generate_population(self, ga_optim_default):
        n = 42
        population = ga_optim_default.generate_population(n)
        assert len(population) == n
        assert all(
            ga_optim_default.is_solution_inside_bounds(sol) for sol in population
        )

    def test_update_population(self, ga_optim_default):
        n = 42
        population = ga_optim_default.generate_population(n)
        assert all(not sol.is_updated() for sol in population)

        ga_optim_default.update_population(population)
        assert all(sol.is_updated() for sol in population)

    def test_filter_population(self, ga_optim_for_is_valid):
        n = 42
        population = ga_optim_for_is_valid.generate_population(n)

        population_filtered = ga_optim_for_is_valid.filter_population(population)
        assert len(population_filtered) == 0  # no one is updated

        bounds = ga_optim_for_is_valid.bounds
        margin = 0.1
        assert np.all(bounds[:, 1] > margin)
        for sol in population:
            sol.x = bounds[:, 1] - margin
            sol.update()
            assert ga_optim_for_is_valid.is_solution_inside_bounds(sol)

        population_filtered = ga_optim_for_is_valid.filter_population(population)
        assert len(population) == len(population_filtered)
        assert population == population_filtered

        population[0].x = population[0].x  # makes this solution not updated
        assert not population[0].is_updated()
        population[1]._x = np.full_like(
            len(bounds), -0.1
        )  # makes this solution invalid, yet updated
        assert not population[1].is_valid() and population[1].is_updated()
        population[2]._x = (
            bounds[:, 1] + 1
        )  # makes this solution outside bounds, yet updated and valid
        assert (
            not ga_optim_for_is_valid.is_solution_inside_bounds(population[2])
            and population[2].is_updated()
            and population[2].is_valid()
        )

        population_filtered = ga_optim_for_is_valid.filter_population(population)
        assert len(population) == len(population_filtered) + 3
        assert population[3:] == population_filtered

    def test_get_mutants(self, ga_optim_with_data):

        for n in range(2):
            with pytest.raises(ValueError):
                population = ga_optim_with_data.generate_population(n)
                ga_optim_with_data.get_mutants(population)

        n = 3
        population = ga_optim_with_data.generate_population(n)
        ga_optim_with_data.update_population(population)
        for i, sol in enumerate(sorted(population)):
            sol["state"] = str(i)

        for size in "1.0", 1.0:
            with pytest.raises(TypeError):
                ga_optim_with_data.get_mutants(population, size)

        with pytest.raises(ValueError):
            ga_optim_with_data.get_mutants(population, -1)

        ga_optim_with_data._crossover_rate = 0
        ga_optim_with_data._mutation_rate = 0

        mutants = ga_optim_with_data.get_mutants(population, 2)
        ga_optim_with_data.update_population(mutants)
        population_sorted = sorted(population)
        assert min(mutants) == population_sorted[0] and min(mutants)["state"] == "0"
        assert max(mutants) == population_sorted[1] and max(mutants)["state"] == "1"

        ga_optim_with_data._crossover_rate = 1
        mutants = ga_optim_with_data.get_mutants(population, 2)
        ga_optim_with_data.update_population(population)
        for sol in mutants:
            assert sol["state"] == "0"

    def test_get_elites(self, ga_optim_default):

        n = 10
        population = ga_optim_default.generate_population(n)

        for size in 1.0, "1.0":
            with pytest.raises(TypeError):
                ga_optim_default.get_elites(population, size)

        for size in -1, n + 1:
            with pytest.raises(ValueError):
                ga_optim_default.get_elites(population, size)

        ga_optim_default.update_population(population)
        elites = ga_optim_default.get_elites(population, size=n)
        assert len(population) == len(elites)
        assert elites == sorted(elites)

        elites = ga_optim_default.get_elites(population, size=0)
        assert len(elites) == 0

        elites = ga_optim_default.get_elites(population, size=2)
        assert elites[0] <= elites[1]

    def test_transmit_solution_data(self, ga_optim_with_data):
        parent = ga_optim_with_data.generate_solution()
        child = ga_optim_with_data.generate_solution()

        parent["state"] = "parent_state"  # will be transmitted
        parent["spam"] = "parent_spam"

        child["spam"] = "child_spam"
        child["foo"] = "child_foo"

        ga_optim_with_data._transmit_solution_data(sol_parent=parent, sol_child=child)

        assert child["state"] == "parent_state"
        assert child["spam"] == "child_spam"
        assert child["foo"] == "child_foo"

        ga_optim_with_data._keys_data_transmit = ["state", "spam"]
        parent["state"] = "parent_state_another"
        parent["spam"] = ["s", "p", "a", "m"]
        parent["foo"] = 42

        ga_optim_with_data._transmit_solution_data(sol_parent=parent, sol_child=child)

        assert child["state"] == "parent_state_another"
        assert child["spam"] == parent["spam"]
        assert child["foo"] != parent["foo"]

        parent["spam"].append("X")
        assert child["spam"] != parent["spam"]  # child had copy of the `spam`
