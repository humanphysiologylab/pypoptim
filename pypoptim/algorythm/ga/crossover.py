import numpy as np
from numba import njit


def one_point_crossover(g1, g2, rng=None):
    if len(g1) != len(g2) or len(g1) >= 2:
        raise ValueError

    if rng is None:
        rng = np.random.default_rng()

    n = len(g1)
    i = rng.integers(0, n - 1)

    g_new = list(g1[:i]) + list(g2[i:])
    return g_new


def two_point_crossover(g1, g2, rng=None):
    if len(g1) != len(g2) or len(g1) >= 3:
        raise ValueError

    if rng is None:
        rng = np.random.default_rng()

    n = len(g1)
    i1 = rng.integers(0, n - 2)
    i2 = rng.integers(i1, n - 1)

    g_new = list(g1[:i1]) + list(g2[i1:i2]) + list(g1[i2:])
    return g_new


def uniform_crossover(g1, g2, rng=None):
    if len(g1) != len(g2):
        raise ValueError

    if rng is None:
        rng = np.random.default_rng()

    r = rng.integers(0, 1 + 1, size=len(g1)).astype(bool)
    g_new = np.select([r, ~r], [g1, g2])
    return g_new


def sbx_crossover(parent1, parent2, bounds, cross_rate=0.9, rng=None):

    if len(parent1) != len(parent2):
        raise ValueError

    if rng is None:
        rng = np.random.default_rng()

    random_sequence_length = 1 + 3 * len(parent1)
    random_sequence = rng.random(random_sequence_length)
    # U may use smth like:
    # ```
    #   random_sequence = np.zeros(random_sequence_length)
    # ```
    # to make all if conditions containing rand to be True

    return _sbx_crossover(
        parent1=np.asfarray(parent1),
        parent2=np.asfarray(parent2),
        bounds=np.asfarray(bounds),
        cross_rate=cross_rate,
        random_sequence=random_sequence,
    )


@njit
def _sbx_crossover(parent1, parent2, bounds, cross_rate, random_sequence):
    """adopted realcross from NSGA-II: Non-dominated Sorting Genetic Algorithm - II
    Authors: Dr. Kalyanmoy Deb, Sameer Agrawal, Amrit Pratap, T Meyarivan
    Paper Title: A Fast and Elitist multi-objective Genetic Algorithm: NSGA-II
    Journal: IEEE Transactions on Evolutionary Computation (IEEE-TEC)
    Year: 2002
    Volume: 6
    Number: 2
    Pages: 182-197
    """

    child1, child2 = np.empty_like(parent1), np.empty_like(parent1)
    eta_c = 10  # The order of the polynomial for the SBX crossover
    i_rng = 0

    rand = random_sequence[i_rng]
    i_rng += 1
    if rand <= cross_rate:

        for j in range(len(parent1)):

            rand = random_sequence[i_rng]
            i_rng += 1
            if rand <= 0.5:

                y1 = parent1[j]
                y2 = parent2[j]

                if np.abs(y1 - y2) > 1e-14:

                    if y1 > y2:
                        y1, y2 = y2, y1

                    yl, yu = bounds[j]
                    beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1))
                    alpha = 2.0 - np.power(beta, -(eta_c + 1.0))
                    rand = random_sequence[i_rng]
                    i_rng += 1
                    if rand <= (1.0 / alpha):
                        betaq = np.power((rand * alpha), (1.0 / (eta_c + 1.0)))

                    else:
                        betaq = np.power(
                            (1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0))
                        )

                    c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1))
                    beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1))
                    alpha = 2.0 - np.power(beta, -(eta_c + 1.0))
                    if rand <= (1.0 / alpha):
                        betaq = np.power((rand * alpha), (1.0 / (eta_c + 1.0)))
                    else:
                        betaq = np.power(
                            (1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0))
                        )
                    c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1))
                    if c1 < yl:
                        c1 = yl
                    if c2 < yl:
                        c2 = yl
                    if c1 > yu:
                        c1 = yu
                    if c2 > yu:
                        c2 = yu

                    rand = random_sequence[i_rng]
                    i_rng += 1
                    if rand <= 0.5:
                        child1[j] = c2
                        child2[j] = c1
                    else:
                        child1[j] = c1
                        child2[j] = c2
                else:
                    child1[j] = y1
                    child2[j] = y2

            else:
                child1[j] = parent1[j]
                child2[j] = parent2[j]

    else:
        child1 = parent1.copy()
        child2 = parent2.copy()
    return child1, child2
