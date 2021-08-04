import copy
from typing import final

import numpy as np


class Solution:
    def __init__(self, x, **kwargs_data):

        x = np.asfarray(x)
        if x.ndim != 1 or x.shape[0] == 0:
            raise ValueError

        self._x = x
        self._y = None
        self._data = copy.deepcopy(kwargs_data)

    def __repr__(self):
        s = "Solution = {\n"
        s += "    x    = {},\n".format(self.x)
        s += "    y    = {},\n".format(self._y)
        s += "    data = {}\n".format(self._data)
        s += "}"
        return s

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if not (self.is_updated() and other.is_updated()):
            raise ValueError("Can't compare invalid Solutions")
        return self.y == other.y

    def __gt__(self, other):
        if not (self.is_updated() and other.is_updated()):
            raise ValueError("Can't compare invalid Solutions")
        return self.y > other.y

    def __ne__(self, other):
        return not self.__eq__(other)

    def __le__(self, other):
        return not self.__gt__(other)

    def __ge__(self, other):
        return self.__gt__(other) or self.__eq__(other)

    def __lt__(self, other):
        return not self.__ge__(other)

    def __contains__(self, item):
        return item in self.data

    def __getitem__(self, item):
        if not self.__contains__(item):
            raise KeyError
        else:
            return self.data[item]

    def __setitem__(self, key, value):
        self._data[key] = value

    def __len__(self):
        return len(self._x)

    # x
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x_new):
        x_new = np.asfarray(x_new)
        if x_new.ndim != 1 or x_new.shape[0] == 0:
            raise ValueError
        self._x = x_new
        self._y = None

    # y
    @property
    def y(self):
        return self._y

    # data
    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        if not isinstance(data, dict):
            raise ValueError
        self._data = data

    def update(self, *args, **kwargs) -> None:
        raise NotImplementedError("You must implement this method on your side!")

    @final
    def is_updated(self) -> bool:
        return self._y is not None

    def is_valid(self) -> bool:
        raise NotImplementedError("You must implement this method on your side!")
