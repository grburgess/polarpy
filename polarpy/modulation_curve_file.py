import h5py
import numpy as np

from threeML.utils.polarization.binned_polarization import BinnedModulationCurve
#from polarpy.polarlike import PolarLike

class ModulationCurveFile(object):

    def __init__(self, counts, scattering_bins, exposures, count_errors=None, sys_errors=None, scale_factor=1.,
                 mission=None, instrument=None, tstart=None, tstop=None):
        """

        A int for modulation curve to facilitate  the reading and writing of files

        :param counts: a matrix of counts where the rows are time intervals and the columns are scattering bins
        :param scattering_bins: and array of scattering bin edges
        :param exposures: and array of exposures for each time bin
        :param count_errors: a matrix of counts where the rows are time intervals and the columns are scattering bins
        :param sys_errors: a matrix of counts where the rows are time intervals and the columns are scattering bins
        :param scale_factor: the scale factor of the background
        :param mission: the space mission
        :param instrument: the instrument on the mission
        :param tstart: an array of start times for the intervals
        :param tstop: an array of stop times for the intervals
        """


        # make sure that all the arrays are the correct shape

        counts = np.atleast_2d(counts)
        exposures = np.atleast_1d(exposures)

        assert len(exposures.shape) == 1

        assert len(counts.shape) == 2

        # extract the shape so we know the interval size
        # as well as the number of bins

        n_intervals, n_bins = counts.shape

        assert len(exposures) == n_intervals

        assert len(scattering_bins) == n_bins + 1, 'The shape of the counts is incorrect'

        assert np.all(exposures >= 0.), 'exposures are not positive'

        if count_errors is not None:

            self._is_poisson = False
            count_errors = np.atleast_2d(count_errors)

            assert np.all(count_errors.shape == counts.shape)

        else:

            self._is_poisson = True

        self._count_errors = count_errors

        # correct start and stops
        tmp = [tstart,tstop]

        for i,_ in enumerate(tmp):

            if tmp[i] is not None:
                tmp[i] = np.atleast_1d(tmp[i])
                assert len(tmp[i].shape) == 1
                assert len(tmp[i]) == n_intervals

        ## fix this later
        self._sys_errors = sys_errors

        if instrument is None:
            instrument = 'Unknown'

        if mission is None:
            mission = 'Unknown'

        self._counts = counts
        self._scattering_bins = scattering_bins
        self._exposures = exposures
        self._n_intervals = n_intervals
        self._scale_factor = scale_factor
        self._tstart = tstart
        self._tstop = tstop
        self._instrument = instrument
        self._mission = mission

    @classmethod
    def read(cls, file_name):

        with h5py.File(file_name, 'r') as f:

            n_intervals = int(f.attrs['n_intervals'])
            is_poisson = bool(f.attrs['is_poisson'])
            scattering_bins = f['scattering_bins'].value
            mission = f.attrs['mission']
            instrument = f.attrs['instrument']
            scale_factor = f.attrs['scale_factor']

            counts = []
            exposures = []

            if not is_poisson:

                count_errors = []

            else:

                count_errors = None

            tstart_flag = True
            tstop_flag = True
            tstart = []
            tstop = []

            for interval in range(n_intervals):
                int_grp = f['interval_%d' % interval]

                counts.append(int_grp['counts'].value)
                exposures.append(int_grp.attrs['exposure'])

                if not is_poisson:
                    count_errors.append(int_grp['count_errors'].value)

                if tstart_flag:
                    try:

                        tstart.append(int_grp.attrs['tstart'])

                    except:

                        tstart = None
                        tstart_flag = False

                if tstop_flag:
                    try:

                        tstop.append(int_grp.attrs['tstop'])

                    except:

                        tstop = None
                        tstop_flag = False

        return cls(counts=counts,
                   exposures=exposures,
                   scattering_bins=scattering_bins,
                   count_errors=count_errors,
                   sys_errors=None,
                   scale_factor=scale_factor,
                   mission=mission,
                   instrument=instrument,
                   tstart=tstart,
                   tstop=tstop)

    def writeto(self, file_name):

        with h5py.File(file_name, 'w') as f:

            for interval in range(self._n_intervals):

                int_grp = f.create_group('interval_%d' % interval)

                int_grp.create_dataset('counts', data=self._counts[interval], compression='lzf')
                int_grp.attrs['exposure'] = self._exposures[interval]

                if self._count_errors is not None:
                    int_grp.create_dataset('count_errors', data=self._count_errors[interval], compression='lzf')

                if self._tstart is not None:
                    int_grp.attrs['tstart'] = self._tstart[interval]

                if self._tstop is not None:
                    int_grp.attrs['tstop'] = self._tstop[interval]

            f.create_dataset('scattering_bins', data=self._scattering_bins, compression='lzf')
            f.attrs['n_intervals'] = self._n_intervals
            f.attrs['is_poisson'] = self._is_poisson
            f.attrs['instrument'] = self._instrument
            f.attrs['mission'] = self._mission
            f.attrs['scale_factor'] = self._scale_factor

    def to_binned_modulation_curve(self, interval=0):
        """

        Create a 3ML BinnedModulationCurve from the contents of the interface.
        An interval with c-style indexing must be specified.

        :param interval: The interval to extract. Starting at zero
        :return: A 3ML BinnedModulationCurve
        """

        # this keeps us from trying to call create things
        # not stored

        count_errors = None
        sys_errors = None
        if self._count_errors is not None:
            count_errors = self._count_errors[interval]

        if self._sys_errors is not None:
            sys_errors = self._sys_errors[interval]

        return BinnedModulationCurve(counts=self._counts[interval],
                                     exposure=self._exposures[interval],
                                     abounds=self._scattering_bins,
                                     count_errors=count_errors,
                                     sys_errors=sys_errors,
                                     scale_factor=self._scale_factor,
                                     is_poisson=self._is_poisson,
                                     instrument=self._instrument,
                                     tstart=self._tstart,
                                     tstop=self._tstop)

    @classmethod
    def from_binned_modulation_curve(cls, binned_mod_curve):
        """
        Create a modulation curve file from a binned modulation curve
        instance. This is really a simple pass thru for writing to
        disk
        :param binned_mod_curve: a 3ML BinnedModulationCurve
        :return: a modulation curve fit file memory instance

        """
        return cls(counts=binned_mod_curve.counts,
                   exposures=binned_mod_curve.exposure,
                   scattering_bins=binned_mod_curve.edges,
                   count_errors=binned_mod_curve.count_errors,
                   sys_errors=binned_mod_curve.sys_errors,
                   scale_factor=binned_mod_curve.scale_factor,
                   mission=binned_mod_curve.mission,
                   instrument=binned_mod_curve.instrument,
                   tstart=binned_mod_curve.tstart,
                   tstop=binned_mod_curve.tstop)

    @classmethod
    def from_list_of_plugins(cls, plugins, kind='total'):

        assert (kind=='total') or (kind=='background'), 'invalid kind. Must be totla or background'

        out= dict(tstart = [],
                  tstop  = [],
                  counts = [],
                  count_errors =[],
                  sys_errors = [],
                  exposures = [],
        )
        for plugin in plugins:

#            assert isinstance(plugin, PolarLike), 'must be a PolarLike plugin'

            if kind=='total':

                bmc = plugin.observation

            elif kind =='background':

                bmc = plugin.background

            
            
            out['tstart'].append(bmc.tstart)
            out['tstop'].append(bmc.tstop)
            out['counts'].append(bmc.counts)
            out['count_errors'].append(bmc.count_errors)
            out['sys_errors'].append(bmc.sys_errors)
            out['exposures'].append(bmc.exposure)

        
            
        for k, v in out.items():

            if np.all(np.array(v) == None):
                out[k] = None
                

            
        return cls(counts=out['counts'],
                   exposures=out['exposures'],
                   scattering_bins=bmc.edges,
                   count_errors=out['count_errors'],
                   sys_errors=out['sys_errors'],
                   scale_factor=bmc.scale_factor,
                   mission=bmc.mission,
                   instrument=bmc.instrument,
                   tstart=out['tstart'],
                   tstop=out['tstop']
        )

        
    
    @property
    def counts(self):
        return self._counts

    @property
    def scattering_bins(self):
        return self._scattering_bins

    @property
    def exposures(self):
        return self._exposures


    
