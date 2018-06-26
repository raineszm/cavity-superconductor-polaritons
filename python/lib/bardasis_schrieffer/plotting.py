import numpy as np


class EigenPlotter:
    def __init__(self, idxs, evals):
        self.idxs = idxs
        self._evals = np.array(evals)
        self.row, self.col = np.indices(self._evals.shape)

    @property
    def evals(self):
        return self.ordered(self._evals)

    def reorder(self,):
        self.row, self.col = np.indices(self._evals.shape)

        for i in range(2, self._evals.shape[0]):
            dx = (
                self._evals[i - 1, self.col[i - 1, :]]
                - self._evals[i - 2, self.col[i - 2, :]]
            ) / (self.idxs[i - 1] - self.idxs[i - 2])
            for k in range(3):
                j = np.argmin(
                    np.abs(
                        (self._evals[i, k] - self._evals[i - 1, self.col[i - 1, :]])
                        / (self.idxs[i] - self.idxs[i - 1])
                        - dx
                    )
                )
                self.col[i, j] = k

    def ordered(self, ary):
        return ary[self.row, self.col]
