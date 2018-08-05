import numpy as np
from polarpy.modulation_curve_file import ModulationCurveFile




def test_constructor():

    scattering_bins = np.linspace(0,360,5)
    counts = np.ones(4)

    count_errors = np.ones(4) * 0.5

    m = ModulationCurveFile(counts=counts,
                            scattering_bins=scattering_bins,
                            exposures=1.


    )

    scattering_bins = np.linspace(0,360,5)
    counts = np.ones((2,4))

    count_errors = np.ones((2,4)) * 0.5

    m = ModulationCurveFile(counts=counts,
                            scattering_bins=scattering_bins,
                            exposures=[1.,1.]


    )

def test_mod_curve_file_writing():

    
    scattering_bins = np.linspace(0,360,5)
    counts = np.ones(4)

    count_errors = np.ones(4) * 0.5

    m = ModulationCurveFile(counts=counts,
                            scattering_bins=scattering_bins,
                            exposures=1.


    )

    m.writeto('test_mc.h5')

    
