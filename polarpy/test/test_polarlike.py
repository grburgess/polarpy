import numpy as np
import os

from polarpy.polar_response import PolarResponse
from polarpy.polarlike import PolarLike
from polarpy.polar2hdf5 import polar_polarization_to_hdf5

from threeML.io.file_utils import sanitize_filename
data_path = sanitize_filename(os.environ.get('POLAR_TEST_DATA_DIR'), abspath=True)

from astromodels import LinearPolarization, Powerlaw, Model, PointSource, SpectralComponent

root_file = os.path.join(data_path, 'fold_spec.root')

outfile = 'testrsp.h5'

polar_polarization_to_hdf5(polarization_root_file=root_file, hdf5_out_file=outfile)

# now make sure the polar response works

prsp = PolarResponse(outfile)


def test_countructor():

    obs = os.path.join(data_path, 'test_mcf.h5')
    bak = os.path.join(data_path, 'test_mcf.h5')

    polarlike = PolarLike('test', observation=obs, background=bak, response=prsp, interval_number=1)


def test_setting_model():

    obs = os.path.join(data_path, 'test_mcf.h5')
    bak = os.path.join(data_path, 'test_mcf.h5')

    polarlike = PolarLike('test', observation=obs, background=bak, response=prsp, interval_number=1)

    pl = Powerlaw()
    pz = LinearPolarization(10, 10)

    sc = SpectralComponent(pl, pz)

    ps = PointSource('test', 0, 0, spectral_components=[sc])
    model = Model(ps)

    polarlike.set_model(model)

    polarlike.get_log_like()

    sim = polarlike.get_simulated_dataset()

    assert isinstance(sim, PolarLike)
