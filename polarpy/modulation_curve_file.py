import h5py
import numpy as np

from threeML.utils.polarization.binned_polarization import BinnedModulationCurve

class ModulationCurveFile(object):

    def __init__(self, counts, scattering_bins, exposures, count_errors=None, sys_errors=None,  scale_factor=1.,mission=None, instrument=None, tstart=None, tstop=None):

        counts = np.atleast_2d(counts)
        exposures = np.atleast_1d(exposures)

        assert len(exposures.shape) == 1 

        assert len(counts.shape) == 2

        n_intervals, n_bins = counts.shape

        assert len(exposures) == n_intervals
        
        assert len(scattering_bins) == n_bins + 1, 'The shape of the counts is incorrect' 

        assert np.all(exposures >= 0.), 'exposures are not positive'


        if count_errors is not None:

            self._is_poisson = False
            count_errors = np.atleast_2d(count_errors)

            assert np.all(count_errors.shape == counts.shape)

        if tstart is not None:

            tstart = np.atleast_1d(tstart)
            assert len(tstart.shape) == 1
            assert len(tstart) == n_intervals

        if tstop is not None:

            tstop = np.atleast_1d(tstop)
            assert len(tstop.shape) == 1
            assert len(tstop) == n_intervals
 
            

        ## fix this later
        self._sys_errors = sys_errors
        
            
        
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
            tstop =[]
                
            for interval in range(n_intervals):
                int_grp = f['interval_%d' % interval]

                counts.append(int_grp['counts'].value)
                exposures.append(int_grp['exposure'].value)

                if not is_poisson:

                    count_errors.append(int_grp['count_errors'].value)

                if tstart_flag:
                    try:

                        tstart.append(int_grp['tstart'].value)

                    except:

                        tstart = None
                        tstart_flag = False

                if tstop_flag:
                    try:

                        tstop.append(int_grp['tstop'].value)

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

            
                int_grp.create_dataset('counts',data = self._counts[interval], compression='lzf' )
                int_grp.create_dataset('exposure',data = self._exposure[interval], compression='lzf' )

                if self._count_errors is not None:
                    int_grp.create_dataset('count_errors',data = self._count_errors[interval], compression='lzf' )

                if self._tstart is not None:
                    int_grp.create_dataset('tstart',data = self._tstart[interval], compression='lzf' )

                if self._tstop is not None:
                    int_grp.create_dataset('tstop',data = self._tstop[interval], compression='lzf' )
                

                
            f.create_dataset('scattering_bins',data = self._scattering_bins, compression='lzf')
            f.attrs['n_intervals'] = self._n_intervals
            f.attrs['is_poisson'] = self._is_poisson
            f.attrs['instrument'] = self._instrument
            f.attrs['mission'] = self._mission

    def to_binned_modulation_curve(self, interval=0):

        return BinnedModulationCurve(counts=self._counts[interval],
                                     exposure=self._exposures[interval],
                                     abounds=self._scattering_bins,
                                     count_errors = self._count_errors[interval],
                                     sys_errors = self._sys_errors[interval],
                                     scale_factors=self._scale_factor,
                                     is_poisson=self._is_poisson,
                                     instrument = self._instrument,
                                     tstart = self._tstart,
                                     tstop = self._tstop)
    @classmethod
    def from_binned_modulation_curve(cls, binned_mod_curve):
        """
        Create a modulation curve file from a binned modulation curve
        instance. This is really a simple pass thru for writing to
        disk

        """
        return cls(counts=binned_mod_curve.counts,
                   exposures=binned_mod_curve.exposures,
                   scattering_bins=binned_mod_curve.abounds,
                   count_errors=binned_mod_curve.count_errors,
                   sys_errors=binned_mod_curve.sys_errors,
                   scale_factor=binned_mod_curve.scale_factor,
                   mission=binned_mod_curve.mission,
                   instrument=binned_mod_curve.instrument,
                   tstart=binned_mod_curve.tstart,
                   tstop=binned_mod_curve.tstop)


        
        return cls(counts=binned_mod_curve.counts,
                   exposure=binned_mod_curve.exposure


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
