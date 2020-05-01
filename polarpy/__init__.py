from polarpy.polar_response import PolarResponse
from polarpy.polarlike import PolarLike

from ._version import get_versions

try:
    import ROOT

    has_root = True

except ImportError:

    has_root = False

if has_root:
    from polarpy.polar2hdf5 import polar_polarization_to_hdf5, polar_spectra_to_hdf5

__version__ = get_versions()['version']
del get_versions
