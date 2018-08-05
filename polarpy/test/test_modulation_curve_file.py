import numpy as np
import os
from polarpy.modulation_curve_file import ModulationCurveFile

from threeML.io.file_utils import sanitize_filename
data_path = sanitize_filename(os.environ.get('POLAR_TEST_DATA_DIR'),abspath=True)

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

    write_path = os.path.join(data_path,'test_write_mc.h5')
    
    m.writeto(write_path)

    
