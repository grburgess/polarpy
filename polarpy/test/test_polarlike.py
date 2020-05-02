import os

import numpy as np

from astromodels import (LinearPolarization, Model, PointSource, Powerlaw,
                         SpectralComponent)
from polarpy.modulation_curve_file import ModulationCurveFile
from polarpy.polar2hdf5 import polar_polarization_to_hdf5
from polarpy.polar_response import PolarResponse
from polarpy.polarlike import PolarLike
from threeML.io.file_utils import sanitize_filename
from threeML.utils.polarization.binned_polarization import \
    BinnedModulationCurve

data_path = sanitize_filename(os.environ.get(
    'POLAR_TEST_DATA_DIR'), abspath=True)


root_file = os.path.join(data_path, 'fold_spec.root')

outfile = 'testrsp.h5'

polar_polarization_to_hdf5(
    polarization_root_file=root_file, hdf5_out_file=outfile)

# now make sure the polar response works

prsp = PolarResponse(outfile)


def test_countructor():

    obs = os.path.join(data_path, 'test_mcf.h5')
    bak = os.path.join(data_path, 'test_mcf.h5')

    polarlike = PolarLike('test', observation=obs, background=bak,
                          response=prsp, interval_number=1, verbose=True)

    assert isinstance(polarlike._observation, BinnedModulationCurve)
    assert isinstance(polarlike._background, BinnedModulationCurve)

    scattering_bins = np.linspace(0, 360, 31)
    counts = np.ones(30)

    count_errors = np.ones(30) * 0.5

    m = ModulationCurveFile(counts=counts,
                            scattering_bins=scattering_bins,
                            exposures=1.


                            )

    obs = m.to_binned_modulation_curve(interval=0)

    bak = m.to_binned_modulation_curve(interval=0)

    polarlike = PolarLike('test', observation=obs,
                          background=bak, response=prsp)

    assert isinstance(polarlike._observation, BinnedModulationCurve)
    assert isinstance(polarlike._background, BinnedModulationCurve)

    # now test with rsp as file neame

    polarlike = PolarLike('test', observation=obs,
                          background=bak, response=outfile)

    polarlike.use_effective_area_correction()

    polarlike.fix_effective_area_correction(1.)


def test_setting_model():

    obs = os.path.join(data_path, 'test_mcf.h5')
    bak = os.path.join(data_path, 'test_mcf.h5')

    polarlike = PolarLike('test', observation=obs,
                          background=bak, response=prsp, interval_number=1)

    pl = Powerlaw()
    pz = LinearPolarization(10, 10)

    sc = SpectralComponent('pl', pl, pz)

    ps = PointSource('test', 0, 0, components=[sc])
    model = Model(ps)

    polarlike.set_model(model)

    polarlike.get_log_like()

    polarlike.display()

    polarlike.display(show_total=True)


def test_rebin():

    obs = os.path.join(data_path, 'test_mcf.h5')
    bak = os.path.join(data_path, 'test_mcf.h5')

    polarlike = PolarLike('test', observation=obs,
                          background=bak, response=prsp, interval_number=1)

    pl = Powerlaw()
    pz = LinearPolarization(10, 10)

    sc = SpectralComponent('pl', pl, pz)

    ps = PointSource('test', 0, 0, components=[sc])
    model = Model(ps)

    polarlike.set_model(model)

    polarlike.rebin_on_background(5)

    polarlike.get_log_like()

#    polarlike.display()

#    polarlike.display(show_total=True)

    polarlike.remove_rebinning()

    polarlike.rebin_on_source(5)

    polarlike.get_log_like()

#    polarlike.display()

#    polarlike.display(show_total=True)

    polarlike.remove_rebinning()


def test_simulations():

    scattering_bins = np.linspace(0, 360, 31)
    counts = np.ones(30)
    count_errors = np.ones(30) * 0.5

    m = ModulationCurveFile(counts=counts,
                            scattering_bins=scattering_bins,
                            exposures=1.
                            )

    obs = m.to_binned_modulation_curve(interval=0)

    bak = m.to_binned_modulation_curve(interval=0)

    polarlike = PolarLike('test', observation=obs,
                          background=bak, response=prsp)

    pl = Powerlaw()
    pz = LinearPolarization(10, 10)

    sc = SpectralComponent('pl', pl, pz)

    ps = PointSource('test', 0, 0, components=[sc])
    model = Model(ps)

    polarlike.set_model(model)

    polarlike.get_log_like()

    sim = polarlike.get_simulated_dataset()

    assert isinstance(sim, PolarLike)

    m = ModulationCurveFile(counts=counts,
                            scattering_bins=scattering_bins,
                            exposures=1.,
                            count_errors=count_errors


                            )

    bak = m.to_binned_modulation_curve(interval=0)

    polarlike = PolarLike('test', observation=obs,
                          background=bak, response=prsp, interval_number=0)

    pl = Powerlaw()
    pz = LinearPolarization(10, 10)

    sc = SpectralComponent('pl', pl, pz)

    ps = PointSource('test', 0, 0, components=[sc])
    model = Model(ps)

    polarlike.set_model(model)

    polarlike.get_log_like()

    sim = polarlike.get_simulated_dataset()

    assert isinstance(sim, PolarLike)


def test_writing():

    scattering_bins = np.linspace(0, 360, 31)
    counts = np.ones(30)
    count_errors = np.ones(30) * 0.5

    m = ModulationCurveFile(counts=counts,
                            scattering_bins=scattering_bins,
                            exposures=1.
                            )

    obs = m.to_binned_modulation_curve(interval=0)

    bak = m.to_binned_modulation_curve(interval=0)

    polarlike = PolarLike('test', observation=obs,
                          background=bak, response=prsp)

    polarlike.writeto('dummy')
