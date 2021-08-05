import numpy as np
import pytest

from ..helpers import (
    argmax,
    argmin,
    calculate_reflection,
    random_value_from_bounds,
    transform_genes_bounds,
    transform_genes_bounds_back,
    uniform_vector,
)


def test_uniform_vector():

    n = 42
    uv = uniform_vector(n)

    assert len(uv) == n
    assert np.isclose(np.linalg.norm(uv), 1)


def test_random_value_from_bounds():
    bounds_invalid = [[0, 0], [3, 1], [-2, -1]]
    log_scale = [None, None, True]
    for bounds, log in zip(bounds_invalid, log_scale):
        with pytest.raises(ValueError):
            random_value_from_bounds(bounds=bounds, log_scale=log)

    bounds_correct = [1, 2]
    for log_scale in True, False:
        value = random_value_from_bounds(bounds=bounds_correct, log_scale=log_scale)
        assert bounds_correct[0] <= value <= bounds_correct[1]


def test_transform_genes():

    genes = np.array([1.0, 1.0, -1.0, 100.0])

    bounds = np.array([[0.0, 0.1, -3.0, 10.0], [2.0, 10.0, 0.0, 1000.0]]).T

    gammas = np.array([1.0, 1.0, 2.0, 0.5])

    mask_multipliers = np.array([False, True, False, True])

    genes_transformed, bounds_transformed = transform_genes_bounds(
        genes, bounds, gammas, mask_multipliers
    )

    assert np.allclose(genes_transformed, np.array([0.25, 0.25, 0.16666667, 0.5]))

    assert np.all(bounds_transformed[:, 0] == 0)

    assert np.allclose(bounds_transformed[:, 1], np.array([0.5, 0.5, 0.25, 1.0]))

    genes_back = transform_genes_bounds_back(
        genes_transformed, bounds_transformed, bounds, mask_multipliers
    )

    assert np.allclose(genes, genes_back)


@pytest.mark.parametrize(
    "shifts,genes_expected,ub,lb,values",
    [
        (0.25, 0.75, 1, 0, 0.5),
        (1.3, 0.2, 1, 0, 0.5),
        (-5.3, 0.8, 1, 0, 0.5),
    ],  # TODO: add multidimensional test cases
)
def test_reflection(shifts, genes_expected, ub, lb, values):
    genes_after = calculate_reflection(ub=ub, lb=lb, values=values, shifts=shifts)
    assert np.allclose(genes_after, genes_expected)


@pytest.mark.xfail
def test_strip_comments():
    assert 0


@pytest.mark.xfail
def test_autoscaling():
    assert 0


@pytest.mark.parametrize(
    "test_input,expected_output", [([0], 0), ([0, 1], 0), ([1, 0], 1), ([0, 1, 0], 0)]
)
def test_argmin(test_input, expected_output):
    assert argmin(test_input) == expected_output


@pytest.mark.parametrize(
    "test_input,expected_output", [([0], 0), ([0, 1], 1), ([1, 0], 0), ([1, 0, 1], 0)]
)
def test_argmax(test_input, expected_output):
    assert argmax(test_input) == expected_output


@pytest.mark.xfail
def test_calculate_gammas_from_bounds():
    assert 0
