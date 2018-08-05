from threeML.io.file_utils import sanitize_filename
import os


data_path = sanitize_filename(os.environ.get('POLAR_TEST_DATA_DIR'),abspath=True)



def test_response_conversion():

    root_file = os.path.join(data_path,'fold_spec.root')
    
    outfile = 'testrsp.h5'
    
    polar_polarization_to_hdf5(polarization_root_file=root_file, hdf5_out_file=outfile)
