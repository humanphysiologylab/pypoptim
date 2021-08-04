import numpy as np
import pytest

from pypoptim.algorythm.ga import GA


@pytest.fixture()
def ga_optim_fabric(square_solution):
    def _ga_optim_fabric(**kw):
        bounds = np.asfarray([[-1, 1], [2, 4]])
        rng = np.random.default_rng(42)
        return GA(SolutionSubclass=square_solution, bounds=bounds, rng=rng, **kw)

    return _ga_optim_fabric


@pytest.fixture()
def ga_optim_default(ga_optim_fabric):
    return ga_optim_fabric()


@pytest.fixture()
def ga_optim_with_data(ga_optim_fabric):
    keys_data_transmit = ["state"]
    return ga_optim_fabric(keys_data_transmit=keys_data_transmit)


@pytest.fixture()
def ga_optim_for_is_valid(ga_optim_fabric):
    ga_optim = ga_optim_fabric()

    class SS(ga_optim._SolutionSubclass):
        def is_valid(self):
            return np.all(self.x > 0)

    ga_optim._SolutionSubclass = SS
    return ga_optim


@pytest.fixture()
def population(square_solution):
    def _population(n, bounds, **kw):
        ga_optim = GA(SolutionSubclass=square_solution, bounds=bounds, **kw)
        return ga_optim.generate_population(n)

    return _population
