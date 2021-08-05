# pypoptim
>Population-based algorithms for global optimization of the cardiac action potential models.

## Getting started

### Install
```sh
pip install .
# or pip install -e .
```

### Test
```sh
pytest --pyargs pypoptim
```

### Usage
See [examples/](./examples) folder.

### Pre-commit
1. Install:
```sh
pip install pre-commit && pre-commit install
```
2. Run:
```sh
git add FILES
pre-commit run
git add FILES_FIXED  # if necessary
git commit -m "MESSAGE"
```
