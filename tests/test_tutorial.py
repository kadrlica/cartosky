#!/usr/bin/env python
"""
Test notebooks in tutorial.

Adapted from:
https://blog.thedataincubator.com/2016/06/testing-jupyter-notebooks/
"""
__author__ = "Alex Drlica-Wagner"

import os
import unittest
import glob

import subprocess
import tempfile

import nbformat

def _notebook_run(path):
    """Execute a notebook via nbconvert and collect output.
       :returns (parsed nb object, execution errors)
    """
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
          "--ExecutePreprocessor.timeout=60", "--log-level=WARN",
          "--output", fout.name, path]
        subprocess.check_call(args)

        fout.seek(0)
        nb = nbformat.read(fout, nbformat.current_nbformat)

    errors = [output for cell in nb.cells if "outputs" in cell
                     for output in cell["outputs"]\
                     if output.output_type == "error"]

    return nb, errors


class TestTutorial(unittest.TestCase):
    """ Test the tutorial notebooks. """

    path = './tutorial'
    
    def _test_notebook(self, path):
        """ Test a notebook. 

        Parameters
        ----------
        path : path to notebook

        Returns
        -------
        None
        """
        nb, errors = _notebook_run(path)
        assert errors == []
        
    def _test_tutorial(self, i=1):
        """ Test a specific tutorial of the tutorial. """
        path = glob.glob(self.path+'/tutorial%i_*.ipynb'%i)[0]
        abspath = os.path.abspath(path)
        self._test_notebook(abspath)

    def test_tutorial1(self):
        return self._test_tutorial(1)

    def test_tutorial2(self):
        return self._test_tutorial(2)

    def test_tutorial3(self):
        return self._test_tutorial(3)

    def test_tutorial4(self):
        return self._test_tutorial(4)

    def test_tutorial5(self):
        return self._test_tutorial(5)


if __name__ == "__main__":
    unittest.main()
