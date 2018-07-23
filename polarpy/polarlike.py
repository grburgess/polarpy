import numpy as np
import scipy.interpolate as interpolate
import 

from threeML import XYLike, PluginPrototype
from threeML.utils.statistics.likelihood_functions import poisson_observed_poisson_background

class PolarLike(PluginPrototype):
    """
    Preliminary POLAR polarization plugin
    """
    
    
    def __init__(self,name, observation, background, response, exposure=1., background_exposure=1.):
        """
        
        """
    

        
        self._total_counts = observation
        self._background_counts = background
        self._scale = exposure/background_exposure
        
        
        self._interpolate_rsp(rsp_file)
        
        
        self._nuisance_parameter = Parameter("cons_%s" % name, 1.0, min_value=0.8, max_value=1.2, delta=0.05,
                                             free=False, desc="Effective area correction for %s" % name)

        nuisance_parameters = collections.OrderedDict()
        nuisance_parameters[self._nuisance_parameter.name] = self._nuisance_parameter
         
         
        super(PolarLike, self).__init__(name, nuisance_parameters)
        
    
    def _interpolate_rsp(self, rsp_file):
        """
        Builds the interpolator for the response. This is currently incredibly slow
        and should be improved

        """

        
        # create a functions to get the histograms
        def get_hist(ene,degree,angle):
    
    
            with h5py.File(rsp_file,'r') as f:

                tmp =  np.array(f['ene_%d'% int(ene)]['deg_%d'% int(degree)]['ang_%d'% int(angle)].value)

                bins = np.array(f['bins'].value)

            return bins, tmp
        
        
        # now go through the response and extract things
        
        with h5py.File(rsp_file,'r') as f:
        
            ene = [int(str(x[4:])) for x in f.keys() if 'ene'  in x]
            
            
            energy = np.sort(np.array(ene))
            
            ene_lo, ene_hi = [],[]

            for ene in energy:

                ene_lo.append(ene-2.5)
                ene_hi.append(ene+2.5)
            
            
            pol_ang = np.array(f['pol_ang'].value)
        
            pol_deg = np.array(f['pol_deg'].value)
        
        
            bins = np.array(f['bins'].value)
            
            bins +=12

            
            bin_center = 0.5 *(bins[:-1] + bins[1:])

            all_interp = []
            
            for i, bm in enumerate(bin_center):
    
                data = []
                #energy = get_energy()

                for ene in energy:
                    for ang in pol_ang:
                        for deg in pol_deg:

                            _, hist = get_hist(ene,deg,ang)


                            data.append(hist[i])
                data = np.array(data).reshape((energy.shape[0],
                                               pol_ang.shape[0],
                                               pol_deg.shape[0]))




                this_interpolator = interpolate.RegularGridInterpolator((energy,pol_ang,pol_deg), data)


                all_interp.append(this_interpolator)
                
            self._all_interp = all_interp
            
            self._ene_lo = ene_lo
            self._ene_hi = ene_hi
            self._energy_mid = energy
            
            self._n_scatter_bins = len(bin_center)
            self._scattering_bins = bin_center


            
            

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

        
        for k,v in likelihood_model_instance.free_parameters.items():
    
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

            return (e2 - e1) / 6.0 * (differential_flux(e1)
                                      + 4 * differential_flux((e1 + e2) / 2.0)
                                      + differential_flux(e2))

        return differential_flux, integral

   def _get_total_expectation(self):

        # first we need to get the integrated expectation from the spectrum
        
        summed_hist = np.zeros(self._n_scatter_bins)
        intergal_spectrum = np.array([self._integral_flux(emin, emax) for emin, emax in zip(self._ene_lo, self._ene_hi)])
        eval_points = np.array([ [ene, self._pol_angle.value, self._pol_degree.value] for ene in self._energy_mid])
        expectation = []
        
        for i, interpolator in enumerate(self._all_interp):
            
            counts = np.dot(interpolator(eval_points), intergal_spectrum)
        
            expectation.append(counts)
        
    
        return np.array(expectation)
    
    def get_log_like(self):
        
        model_rate = self._get_total_expectation()

        model_counts = self._nuisance_parameter.value * self._exposure


        
        
        
        loglike, bkg_model = poisson_observed_poisson_background(self._total_counts,
                                                                 self._background_counts,
                                                                 self._scale,
                                                                 model_counts)

        return np.sum(loglike)

    def inner_fit(self):

        return self.get_log_like()
