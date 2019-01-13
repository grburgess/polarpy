from polarpy.polarlike import PolarLike
from polarpy.polar_response import PolarResponse

try:
    import ROOT
    
    has_root = True

except ImportError:

    has_root = False

if has_root:
    from polarpy.polar2hdf5 import polar_polarization_to_hdf5, polar_spectra_to_hdf5
