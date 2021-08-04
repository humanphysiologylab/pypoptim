import numpy as np


def tournament_selection(population, selection_force=2, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    if selection_force > len(population):
        msg = f"Selection force must be less than population size: {selection_force} <= {len(population)} is violated"
        raise ValueError(msg)
    return min(rng.choice(a=population, size=selection_force, replace=False))
