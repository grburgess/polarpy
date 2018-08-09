import collections
from contextlib import contextmanager
import matplotlib.pyplot as plt
import numpy as np

from astromodels import Parameter, Uniform_prior
from threeML import PluginPrototype
from threeML.io.plotting.step_plot import step_plot
from threeML.utils.binner import Rebinner
from threeML.utils.polarization.binned_polarization import BinnedModulationCurve
from threeML.utils.statistics.likelihood_functions import poisson_observed_poisson_background, \
    poisson_observed_gaussian_background

from polarpy.modulation_curve_file import ModulationCurveFile
from polarpy.polar_response import PolarResponse


class PolarLike(PluginPrototype):
    """
    Preliminary POLAR polarization plugin
    """

    def __init__(self, name, observation, background, response, interval_number=None, verbose=False):
        """



        :param interval_number:
        :param name:
        :param observation:
        :param background:
        :param response:

        :param verbose:

        """

        # if we pass a string, there may be multiple time intervals
        # saved so we must specify a time interval

        if isinstance(observation, str):
            assert interval_number is not None, 'must specify an interval number'

            # this is a file
            read_file = ModulationCurveFile.read(observation)

            # create the bmc
            observation = read_file.to_binned_modulation_curve(interval=interval_number)

        # the same applies for the background
        if isinstance(background, str):
            assert interval_number is not None, 'must specify an interval number'

            # this is a file
            read_file = ModulationCurveFile.read(background)

            background = read_file.to_binned_modulation_curve(interval=interval_number)

        assert isinstance(observation, BinnedModulationCurve), 'The observation must be a BinnedModulationCurve'
        assert isinstance(background, BinnedModulationCurve), 'The observation must be a BinnedModulationCurve'

        # attach the required variables

        self._observation = observation
        self._background = background

        self._observed_counts = observation.counts
        self._background_counts = background.counts
        self._background_count_errors = background.count_errors
        self._scale = observation.exposure / background.exposure
        self._exposure = observation.exposure
        self._background_exposure = background.exposure

        self._likelihood_model = None
        self._rebinner = None
        
        # now do some double checks

        assert len(self._observed_counts) == len(self._background_counts)

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

        # pass to the plugin proto

        super(PolarLike, self).__init__(name, nuisance_parameters)


        # The following vectors are the ones that will be really used for the computation. At the beginning they just
        # point to the original ones, but if a rebinner is used and/or a mask is created through set_active_measurements,
        # they will contain the rebinned and/or masked versions

        self._current_observed_counts = self._observed_counts
        self._current_background_counts = self._background_counts
        self._current_background_count_errors = self._background_count_errors


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

        # we also make sure the lengths match up here
        assert self._response.n_scattering_bins == len(
            self._observation.counts), 'observation counts shape does not agree with response shape'

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

        
        # We remove the mask temporarily because we need the various elements for all channels. We will restore it
        # at the end

        original_rebinner = self._rebinner

        with self._without_rebinner():

            # Get the source model for all channels (that's why we don't use the .folded_model property)

        
            source_model_counts = self._get_model_counts()

            if self._background.is_poisson:
                _, background_model_counts = poisson_observed_poisson_background(
                            self._current_observed_counts, self._current_background_counts, self._scale, source_model_counts)
            else:

                _, background_model_counts = poisson_observed_gaussian_background(
                            self._current_observed_counts, self._current_background_counts, self._current_background_count_errors, source_model_counts)

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
                verbose=False,
            )

            # Apply the same selections as the current data set
            if original_rebinner is not None:

                # Apply rebinning, which also applies the mask
                new_plugin._apply_rebinner(original_rebinner)

            
            return new_plugin

    def set_model(self, likelihood_model_instance):
        """
        Set the model to be used in the joint minimization. Must be a LikelihoodModel instance.
        :param likelihood_model_instance: instance of Model
        :type likelihood_model_instance: astromodels.Model
        """

        if likelihood_model_instance is None:
            return

        # if self._source_name is not None:

        #     # Make sure that the source is in the model
        #     assert self._source_name in likelihood_model_instance.sources, \
        #                                         "This XYLike plugin refers to the source %s, " \
        #                                         "but that source is not in the likelihood model" % (self._source_name)

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

        n_point_sources = likelihood_model.get_number_of_point_sources()

        # Make a function which will stack all point sources (OGIP do not support spatial dimension)

        def differential_flux(scattering_edges):
            fluxes = likelihood_model.get_point_source_fluxes(0, scattering_edges, tag=self._tag)

            # If we have only one point source, this will never be executed
            for i in range(1, n_point_sources):
                fluxes += likelihood_model.get_point_source_fluxes(i, scattering_edges, tag=self._tag)

            return fluxes

        # The following integrates the diffFlux function using Simpson's rule
        # This assume that the intervals e1,e2 are all small, which is guaranteed
        # for any reasonable response matrix, given that e1 and e2 are Monte-Carlo
        # scattering_edges. It also assumes that the function is smooth in the interval
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

        intergal_spectrum = np.array(
            [self._integral_flux(emin, emax) for emin, emax in zip(self._response.ene_lo, self._response.ene_hi)])

        # we evaluate at the center of the bin. the bin widths are already included
        eval_points = np.array(
            [[ene, self._pol_angle.value, self._pol_degree.value] for ene in self._response.energy_mid])

        expectation = []

        # create the model counts by summing over energy

        for i, interpolator in enumerate(self._all_interp):
            rate = np.dot(interpolator(eval_points), intergal_spectrum)

            expectation.append(rate)

        return np.array(expectation)

    def _get_model_counts(self):


        if self._rebinner is None:
            model_rate = self._get_model_rate()

        else:

            model_rate, = self._rebinner.rebin(self._get_model_rate)

        
        return self._nuisance_parameter.value * self._exposure * model_rate

    def get_log_like(self):

        model_counts = self._get_model_counts()

        if self._background.is_poisson:

            loglike, bkg_model = poisson_observed_poisson_background(        self._current_observed_counts, self._current_background_counts,
                                                                     self._scale, model_counts)

        else:

            loglike, bkg_model = poisson_observed_gaussian_background(        self._current_observed_counts, self._current_background_counts,
                                                                      self._current_background_count_errors, model_counts)

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



    @property
    def scattering_boundaries(self):
        """
        Energy boundaries of channels currently in use (rebinned, if a rebinner is active)

        :return: (sa_min, sa_max)
        """

        scattering_edges = np.array(self._observation.edges)

        sa_min, sa_max = scattering_edges[:-1], scattering_edges[1:]

        if self._rebinner is not None:
            # Get the rebinned chans. NOTE: these are already masked

            sa_min, sa_max = self._rebinner.get_new_start_and_stop(sa_min, sa_max)

 
        return sa_min, sa_max


        
    def display(self, ax=None, show_data=True, show_model=True, show_total=False, model_kwargs={}, data_kwargs={}):
        """

        :param ax:
        :param show_data:
        :param show_model:
        :param show_total:
        :param model_kwargs:
        :param data_kwargs:
        :return:
        """

        sa_min, sa_max = self.scattering_boundaries
        
        if show_total:
            show_model = False
            show_data = False

        if ax is None:

            fig, ax = plt.subplots()

        else:

            fig = ax.get_figure()

        if show_total:

            total_rate =         self._current_observed_counts / self._exposure
            bkg_rate = self._current_background_counts / self._background_exposure

            total_errors = np.sqrt(total_rate)

            if self._background.is_poisson:

                bkg_errors = np.sqrt(bkg_rate)

            else:

                bkg_errors = self._current_background_count_errors

            ax.hlines(
                total_rate,
                sa_min,
                sa_max,
                color='#7D0505',
                **data_kwargs)
            ax.vlines(
                np.mean([self.scattering_boundaries],axis=1),
                total_rate - total_errors,
                total_rate + total_errors,
                color='#7D0505',
                **data_kwargs)

            ax.hlines(
                bkg_rate,
                sa_min,
                sa_max,
                color='#0D5BAE',
                **data_kwargs)
            ax.vlines(
                np.mean([self.scattering_boundaries],axis=1),
                bkg_rate - bkg_errors,
                bkg_rate + bkg_errors,
                color='#0D5BAE',
                **data_kwargs)

        if show_data:

            net_rate = (        self._observed_counts / self._exposure) - self._background_counts / self._background_exposure

            if self._background.is_poisson:

                errors = np.sqrt((        self._observed_counts / self._exposure) +
                                 (self._background_counts / self._background_exposure))

            else:

                errors = np.sqrt((        self._observed_counts / self._exposure) +
                                 (self._background.count_errors / self._background_exposure)**2)

            ax.hlines(net_rate, self._response.scattering_bins_lo, self._response.scattering_bins_hi, **data_kwargs)
            ax.vlines(self._response.scattering_bins, net_rate - errors, net_rate + errors, **data_kwargs)

        if show_model:
            step_plot(
                ax=ax,
                xbins=np.vstack([self._response.scattering_bins_lo, self._response.scattering_bins_hi]).T,
                y=self._get_model_counts() / self._exposure,
                **model_kwargs)

        ax.set_xlabel('Scattering Angle')
        ax.set_ylabel('Net Rate (cnt/s/bin)')

        return fig

    @property
    def observation(self):
        return self._observation

    @property
    def background(self):
        return self._background

    @contextmanager
    def _without__rebinner(self):

        # Store rebinner for later use

        rebinner = self._rebinner

        # Clean mask and rebinning

        self.remove_rebinning()


        # Execute whathever

        yield

        # Restore mask and rebinner (if any)



        if rebinner is not None:

            # There was a rebinner, use it. Note that the rebinner applies the mask by itself

            self._apply_rebinner(rebinner)




            
    def rebin_on_background(self, min_number_of_counts):
        """
        Rebin the spectrum guaranteeing the provided minimum number of counts in each background bin. This is usually
        required for spectra with very few background counts to make the Poisson profile likelihood meaningful.
        Of course this is not relevant if you treat the background as ideal, nor if the background spectrum has
        Gaussian errors.

        The observed spectrum will be rebinned in the same fashion as the background spectrum.

        To neutralize this completely, use "remove_rebinning"

        :param min_number_of_counts: the minimum number of counts in each bin
        :return: none
        """

        # NOTE: the rebinner takes care of the mask already

        assert self._background is not None, "This data has no background, cannot rebin on background!"

        rebinner = Rebinner(self._background_counts, min_number_of_counts, mask = None)

        self._apply_rebinner(rebinner)

    def rebin_on_source(self, min_number_of_counts):
        """
        Rebin the spectrum guaranteeing the provided minimum number of counts in each source bin.

        To neutralize this completely, use "remove_rebinning"

        :param min_number_of_counts: the minimum number of counts in each bin
        :return: none
        """

        # NOTE: the rebinner takes care of the mask already



        rebinner = Rebinner(self._observed_counts, min_number_of_counts, self._mask)

        self._apply_rebinner(rebinner)

    def _apply_rebinner(self, rebinner):

        self._rebinner = rebinner

        # Apply the rebinning to everything.
        # NOTE: the output of the .rebin method are the vectors with the mask *already applied*

        self._current_observed_counts, = self._rebinner.rebin(self._observed_counts)

        if self._background is not None:

            self._current_background_counts, = self._rebinner.rebin(self._background_counts)

            if self._background_count_errors is not None:
                # NOTE: the output of the .rebin method are the vectors with the mask *already applied*

                self._current_background_count_errors, = self._rebinner.rebin_errors(self._background_count_errors)

        if self._verbose:
            print("Now using %s bins" % self._rebinner.n_bins)

    def remove_rebinning(self):
        """
        Remove the rebinning scheme set with rebin_on_background.

        :return:
        """

        self._rebinner = None
