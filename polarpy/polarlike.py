import numpy as np
import scipy.interpolate as interpolate
import h5py

from threeML import PluginPrototype
from threeML.io.plotting.step_plot import step_plot
from threeML.utils.statistics.likelihood_functions import poisson_observed_poisson_background, poisson_observed_gaussian_background
from threeML.utils.polarization.binned_polarization import BinnedModulationCurve
from astromodels import Parameter, Uniform_prior
import matplotlib.pyplot as plt

from polarpy.polar_response import PolarResponse
from polarpy.modulation_curve_file import ModulationCurveFile
import collections


class PolarLike(PluginPrototype):
    """
    Preliminary POLAR polarization plugin
    """

    def __init__(self, name, observation, background, response, interval_number=None, verbose=False):
        """



        :param name:
        :param observation:
        :param background:
        :param response:
        :param background_exposure:
        :param verbose:

        """

        if isinstance(observation, str):

            assert interval_number is not None, 'must specify an interval number'

            # this is a file
            read_file = ModulationCurveFile.read(observation)

            observation = read_file.to_binned_modulation_curve(interval=interval_number)

        if isinstance(background, str):

            assert interval_number is not None, 'must specify an interval number'

            # this is a file
            read_file = ModulationCurveFile.read(background)

            background = read_file.to_binned_modulation_curve(interval=interval_number)

        assert isinstance(observation, BinnedModulationCurve), 'The observation must be a BinnedModulationCurve'
        assert isinstance(background, BinnedModulationCurve), 'The observation must be a BinnedModulationCurve'

        #assert len(observation) == len(background)

        # attach the required variables

        self._observation = observation
        self._background = background

        self._total_counts = observation.counts
        self._background_counts = background.counts
        self._scale = observation.exposure / background.exposure
        self._exposure = observation.exposure
        self._background_exposure = background.exposure

        self._n_synthetic_datasets = 0

        # set up the effective area correction

        self._nuisance_parameter = Parameter(
            "cons_%s" % name,
            1.0,
            min_value=0.8,
            max_value=1.2,
            delta=0.05,
            free=False,
            desc="Effective area correction for %s" % name)

        nuisance_parameters = collections.OrderedDict()
        nuisance_parameters[self._nuisance_parameter.name] = self._nuisance_parameter

        super(PolarLike, self).__init__(name, nuisance_parameters)

        self._source_name = None

        self._verbose = verbose

        # we can either attach or build a response

        assert isinstance(response, str) or isinstance(
            response, PolarResponse), 'The response must be a file name or a PolarResponse'

        if isinstance(response, PolarResponse):

            self._response = response

        else:

            self._response = PolarResponse(response)

        # attach the interpolators to the

        self._all_interp = self._response.interpolators

        assert self._response.n_scattering_bins == len(self._observation.counts), 'observation counts shape does not agree with response shape'

        assert self._response.n_scattering_bins == len(self._background.counts), 'background counts shape does not agree with response shape'

        
    def use_effective_area_correction(self, lower=0.5, upper=1.5):
        """
        Use an area constant to correct for response issues


        :param lower:
        :param upper:
        :return:
        """

        self._nuisance_parameter.free = True
        self._nuisance_parameter.bounds = (lower, upper)
        self._nuisance_parameter.prior = Uniform_prior(lower_bound=lower, upper_bound=upper)
        if self._verbose:

            print('Using effective area correction')

    def fix_effective_area_correction(self, value=1):
        """

        fix the effective area correction to a particular values

        :param value:
        :return:
        """

        # allow the value to be outside the bounds
        if self._nuisance_parameter.max_value < value:

            self._nuisance_parameter.max_value = value + 0.1

        elif self._nuisance_parameter.min_value > value:

            self._nuisance_parameter.min_value = value = 0.1

        self._nuisance_parameter.fix = True
        self._nuisance_parameter.value = value

        if self._verbose:

            print('Fixing effective area correction')

    @property
    def effective_area_correction(self):

        return self._nuisance_parameter

    def get_simulated_dataset(self, new_name=None, **kwargs):
        """
        Returns another Binned instance where data have been obtained by randomizing the current expectation from the
        model, as well as from the background (depending on the respective noise models)

        :return: an BinnedSpectrum or child instance
        """

        assert self._likelihood_model is not None, "You need to set up a model before randomizing"

        # Keep track of how many syntethic datasets we have generated

        self._n_synthetic_datasets += 1

        # Generate a name for the new dataset if needed
        if new_name is None:
            new_name = "%s_sim_%i" % (self.name, self._n_synthetic_datasets)

        # Generate randomized data depending on the different noise models

        # We remove the mask temporarily because we need the various elements for all channels. We will restore it
        # at the end

        # Get the source model for all channels (that's why we don't use the .folded_model property)

        source_model_counts = self._get_model_counts()

        if self._background.is_poisson:
            _, background_model_counts = poisson_observed_poisson_background(
                self._total_counts, self._background_counts, self._scale, source_model_counts)
        else:

            _, background_model_counts = poisson_observed_gaussian_background(
                self._total_counts, self._background_counts, self._background.count_errors, source_model_counts)

        # Now randomize the expectations

        # Randomize expectations for the source

        randomized_source_counts = np.random.poisson(source_model_counts + background_model_counts)

        randomized_background_counts = np.random.poisson(background_model_counts)

        new_observation = self._observation.clone(new_counts=randomized_source_counts)

        new_background = self._background.clone(new_counts=randomized_background_counts)

        new_plugin = PolarLike(
            name=new_name,
            observation=new_observation,
            background=new_background,
            response=self._response,
            verbose=False,)

        return new_plugin

    def set_model(self, likelihood_model_instance):
        """
        Set the model to be used in the joint minimization. Must be a LikelihoodModel instance.
        :param likelihood_model_instance: instance of Model
        :type likelihood_model_instance: astromodels.Model
        """

        if likelihood_model_instance is None:

            return

        if self._source_name is not None:

            # Make sure that the source is in the model
            assert self._source_name in likelihood_model_instance.sources, \
                                                "This XYLike plugin refers to the source %s, " \
                                                "but that source is not in the likelihood model" % (self._source_name)

        for k, v in likelihood_model_instance.free_parameters.items():

            if 'polarization.degree' in k:

                self._pol_degree = v

            if 'polarization.angle' in k:

                self._pol_angle = v

        # now we need to get the intergal flux

        _, integral = self._get_diff_flux_and_integral(likelihood_model_instance)

        self._integral_flux = integral

        self._likelihood_model = likelihood_model_instance

    def _get_diff_flux_and_integral(self, likelihood_model):

        if self._source_name is None:

            n_point_sources = likelihood_model.get_number_of_point_sources()

            # Make a function which will stack all point sources (OGIP do not support spatial dimension)

            def differential_flux(energies):
                fluxes = likelihood_model.get_point_source_fluxes(0, energies, tag=self._tag)

                # If we have only one point source, this will never be executed
                for i in range(1, n_point_sources):
                    fluxes += likelihood_model.get_point_source_fluxes(i, energies, tag=self._tag)

                return fluxes

        else:

            # This SpectrumLike dataset refers to a specific source

            # Note that we checked that self._source_name is in the model when the model was set

            try:

                def differential_flux(energies):

                    return likelihood_model.sources[self._source_name](energies, tag=self._tag)

            except KeyError:

                raise KeyError("This XYLike plugin has been assigned to source %s, "
                               "which does not exist in the current model" % self._source_name)

        # The following integrates the diffFlux function using Simpson's rule
        # This assume that the intervals e1,e2 are all small, which is guaranteed
        # for any reasonable response matrix, given that e1 and e2 are Monte-Carlo
        # energies. It also assumes that the function is smooth in the interval
        # e1 - e2 and twice-differentiable, again reasonable on small intervals for
        # decent models. It might fail for models with too sharp features, smaller
        # than the size of the monte carlo interval.

        def integral(e1, e2):
            # Simpson's rule

            return (e2 - e1) / 6.0 * (differential_flux(e1) + 4 * differential_flux(
                (e1 + e2) / 2.0) + differential_flux(e2))

        return differential_flux, integral

    def _get_model_rate(self):

        # first we need to get the integrated expectation from the spectrum

        summed_hist = np.zeros(self._response.n_scattering_bins)
        intergal_spectrum = np.array(
            [self._integral_flux(emin, emax) for emin, emax in zip(self._response.ene_lo, self._response.ene_hi)])
        eval_points = np.array(
            [[ene, self._pol_angle.value, self._pol_degree.value] for ene in self._response.energy_mid])
        expectation = []

        for i, interpolator in enumerate(self._all_interp):

            rate = np.dot(interpolator(eval_points), intergal_spectrum)

            expectation.append(rate)

        return np.array(expectation)

    def _get_model_counts(self):

        model_rate = self._get_model_rate()

        return self._nuisance_parameter.value * self._exposure * model_rate

    def get_log_like(self):

        model_counts = self._get_model_counts()

        if self._background.is_poisson:

            loglike, bkg_model = poisson_observed_poisson_background(self._total_counts, self._background_counts,
                                                                     self._scale, model_counts)

        else:

            loglike, bkg_model = poisson_observed_gaussian_background(self._total_counts, self._background_counts,
                                                                      self._background.count_errors, model_counts)

        return np.sum(loglike)

    def inner_fit(self):

        return self.get_log_like()

    def writeto(self, file_name):
        """
        Write the data to HDF5 modulation curve files. Both background and observation
        files are created
        :param file_name: the file name header. The .h5 extension is added automatically
        """
        # first create a file container
        observation_file = ModulationCurveFile.from_binned_modulation_curve(self._observation)

        background_file = ModulationCurveFile.from_binned_modulation_curve(self._background)

        observation_file.writeto("%s.h5" % file_name)

        background_file.writeto("%s_bak.h5" % file_name)

    def display(self, ax=None, show_data=True, show_model=True, model_kwargs={}, data_kwargs={}):

        if ax is None:

            fig, ax = plt.subplots()

        else:

            fig = ax.get_figure()

        if show_data:

            net_rate = (self._total_counts / self._exposure) - self._background_counts / self._background_exposure

            errors = np.sqrt((self._total_counts / self._exposure) + (self._background_counts /
                                                                      self._background_exposure))

            ax.hlines(net_rate, self._response.scattering_bins_lo, self._response.scattering_bins_hi, **data_kwargs)
            ax.vlines(self._response.scattering_bins, net_rate - errors, net_rate + errors, **data_kwargs)

            # step_plot(ax=ax,
            #           xbins=np.vstack([self._scattering_bins_lo, self._scattering_bins_hi]).T,
            #           y=net_rate,
            #           **data_kwargs
            # )

        if show_model:

            step_plot(
                ax=ax,
                xbins=np.vstack([self._response.scattering_bins_lo, self._response.scattering_bins_hi]).T,
                y=self._get_model_counts() / self._exposure,
                **model_kwargs)

        ax.set_xlabel('Scattering Angle')
        ax.set_ylabel('Net Rate (cnt/s/bin)')

        return fig
