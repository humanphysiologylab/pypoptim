import numpy as np
import pytest

from ....algorythm.ga.crossover import (
    one_point_crossover,
    sbx_crossover,
    two_point_crossover,
    uniform_crossover,
)


@pytest.mark.xfail
def test_one_point_crossover():
    assert 0


@pytest.mark.xfail
def test_two_point_crossover():
    assert 0


@pytest.mark.xfail
def test_uniform_crossover():
    assert 0


def test_sbx_crossover():

    bounds = [[-3, 3], [-1, 1], [0, 2]]
    parents = [np.array([-1, -0.5, 0]) for _ in range(2)]
    children1, children2 = sbx_crossover(
        parent1=parents[0], parent2=parents[1], bounds=bounds
    )
    assert np.all(children1 == children2) and np.all(children1 == parents[0])

    with pytest.raises(ValueError):
        sbx_crossover(parent1=parents[0][:2], parent2=parents[1], bounds=bounds)

    seed = 888
    parents = [[1.5, 0.1, 1], [0, 0, 0]]
    children1 = sbx_crossover(
        parent1=parents[0],
        parent2=parents[1],
        bounds=bounds,
        rng=np.random.default_rng(seed),
    )
    children2 = sbx_crossover(
        parent1=parents[0],
        parent2=parents[1],
        bounds=bounds,
        rng=np.random.default_rng(seed),
    )

    for child_1, child_2 in zip(children1, children2):
        assert np.all(child_1 == child_2)
