import h5py
import numpy as np


class ModulationCurveFile(object):

    def __init__(self, counts, scattering_bins, exposures):

        counts = np.atleast_2d(counts)
        exposures = np.atleast_1d(exposures)

        assert len(exposures.shape) == 1 

        assert len(counts.shape) == 2

        n_intervals, n_bins = counts.shape

        assert len(exposures) == n_intervals
        
        assert len(scattering_bins) == n_bins + 1, 'The shape of the counts is incorrect' 

        assert exposure >= 0.
        
        self._counts = counts
        self._scattering_bins = scattering_bins
        self._exposures = exposures
        self._n_intervals = n_intervals
        
    @classmethod
    def read(cls, file_name):

        with h5py.File(file_name, 'r') as f:

            n_intervals = int(f.attrs['n_intervals'])

            scattering_bins = f['scattering_bins'].value
            counts = []
            exposures = []
            
            for interval in range(n_intervals):
                int_grp = f['interval_%d' % interval]

                counts.append(int_grp['counts'].value)
                exposures.append(int_grp['exposure'])

        return cls(counts=counts,
                   exposures=exposures,
                   scattering_bins=scattering_bins)
            

    def writeto(self, file_name):

        with h5py.File(file_name, 'w') as f:

            for interval in range(self._n_intervals):

                int_grp = f.create_group('interval_%d' % interval)

            
                int_grp.create_dataset('counts',data = self._counts[interval], compression='lzf' )
                int_grp.create_dataset('exposure',data = self._exposure[interval], compression='lzf' )
                
            f.create_dataset['scattering_bins',data = self._scattering_bins, compression='lzf']
            f.attrs['n_intervals'] = self._n_intervals

    @property
    def counts(self):
        return self._counts

    @property
    def scattering_bins(self):
        return self._scattering_bins

    @property
    def exposures(self):
        return self._exposures
