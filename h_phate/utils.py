import numpy as np
from contextlib import contextmanager

def setdiff1d(x, y):
    if len(x) > 0 and len(y) > 0:
        return np.setdiff1d(x, y)
    else:
        return x
    
class Lock():
    def __init__(self):
        self._lock = False

    @contextmanager
    def lock(self):
        skip_error = ValueError("Global lock is on")

        def lock_fn():
            if self._lock:
                raise skip_error
            else:
                self._lock = True

        try:
            yield lock_fn
            self._lock = False
        except ValueError as err:
            if err != skip_error:
                self._lock = False
                raise
        except:
            self._lock = False
            raise