import numpy as np
from sklearn.metrics import mean_squared_error as MSE


def RMSE(x, y, *, sample_weight=None, multioutput="uniform_average"):
    return MSE(
        x, y, squared=False, sample_weight=sample_weight, multioutput=multioutput
    )


def calculate_RMSE(x, y) -> float:
    assert len(x) == len(y)  # TODO
    return np.sqrt(np.mean(((x - y) ** 2)))


def calculate_RMSE_balanced(x, y) -> float:
    assert len(x) == len(y)  # TODO
    x = (x - y.min(axis=0)) / y.ptp(axis=0)
    y = (y - y.min(axis=0)) / y.ptp(axis=0)
    return np.sqrt(np.sum(((x - y) ** 2)) / len(x))


def calculate_RMSE_weightened(x, y, weights) -> float:
    return float(np.sum(MSE(x, y, squared=False, multioutput="raw_values") * weights))
