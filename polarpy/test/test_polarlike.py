import numpy as np
import os

from polarpy.polar_response import PolarResponse
from polarpy.polarlike import PolarLike
from polarpy.polar2hdf5 import polar_polarization_to_hdf5

from threeML.io.file_utils import sanitize_filename
data_path = sanitize_filename(os.environ.get('POLAR_TEST_DATA_DIR'),abspath=True)



root_file = os.path.join(data_path,'fold_spec.root')
    
outfile = 'testrsp.h5'
    
polar_polarization_to_hdf5(polarization_root_file=root_file, hdf5_out_file=outfile)

# now make sure the polar response works

prsp = PolarResponse(outfile)



def test_countructor():

    obs = os.path.join(data_path,'test_mcf.h5')
    bak = os.path.join(data_path,'test_mcf.h5')
    
    polarlike = PolarLike('test',observation=obs, background=bak,response=prsp, interval_number=1)
    
