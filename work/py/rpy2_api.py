import rpy2
from rpy2 import robjects
from rpy2.robjects import packages
from pathlib import Path


class Rpy2:
    def __init__(self, r_file_path: Path, *args, **kwargs):
        """
        This class is responsible for providing access to rpy2
        """
        self._r_file_path: Path = r_file_path
        self._r_file_read: str = open(self._r_file_path, "r").read()
        self._args = args
        self._kwargs = kwargs

    def call_function(self, r_function_name):
        """
        Call an R function
        Taken in part from: https://gist.github.com/indraniel/da11c4f79c79b5e6bfb8
        """
        func_ = rpy2.robjects.packages.SignatureTranslatedAnonymousPackage(self._r_file_read, "func_")
        call_func_ = getattr(func_, r_function_name)
        results = call_func_(*self._args, **self._kwargs)
        return results