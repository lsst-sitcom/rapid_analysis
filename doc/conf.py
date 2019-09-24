"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documentation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.rapid.analysis


_g = globals()
_g.update(build_package_configs(
    project_name='rapid_analysis',
    version=lsst.rapid.analysis.version.__version__))
