import numpy as np

from threeML.utils.OGIP.response import InstrumentResponse

import h5py


class POLARData(object):
    def __init__(self, polar_hdf5_file, polar_hdf5_response, reference_time=0.):
        """
        container class that converts raw POLAR HDF5 data into useful python
        variables


        :param polar_root_file: path to polar event file
        :param reference_time: reference time of the events (tunix?)
        :param rsp_file: path to rsp file
        """


        

        
        with h5py.File(polar_hdf5_file,'r') as f:

            rsp_grp = f['rsp']

            matrix = rsp_grp['matrix'].value
            ebounds = rsp_grp['ebounds'].value
            mc_low = rsp_grp['mc_low'].value 
            mc_high = rsp_grp['mc_high'].value


            # open the event file
        
            

            # extract the pedestal corrected ADC channels
            # which are non-integer and possibly
            # less than zero
            pha = f['energy'].value

            # non-zero ADC channels are invalid
            idx = pha >= 0
            #pha = pha[idx]

            idx2 = (pha <= ebounds.max()) & (pha >= ebounds.min())

            pha = pha[idx2 & idx]

            # get the dead time fraction
            self._dead_time_fraction = (f['dead_ratio'].value)[idx & idx2]

            # get the arrival time, in tunix of the events
            self._time = (f['time'].value)[idx & idx2] - reference_time

            # digitize the ADC channels into bins
            # these bins are preliminary


            # now do the scattering angles

            scattering_angles = f['scatter_angle'].value

            # clear the bad scattering angles
            idx = scattering_angles != -1

            self._scattering_angle_time = (f['time'].value)[idx] - reference_time
            self._scattering_angle_dead_time_fraction = (f['dead_ratio'].value)[idx]
            self._scattering_angles = scattering_angles[idx]
            

        # build the POLAR response

        mc_energies = np.append(mc_low, mc_high[-1])

        self._rsp = InstrumentResponse(matrix=matrix,
                                       ebounds=ebounds,
                                       monte_carlo_energies=mc_energies)

        # bin the ADC channels

        self._binned_pha = np.digitize(pha, ebounds)

        # bin the scattering_angles

        with h5py.File(polar_hdf5_response, 'r') as f:

            scatter_bounds = f['bins'].value


        self._scattering_bins = scatter_bounds
        self._binned_scattering_angles = np.digitize(self._scattering_angles, scatter_bounds)
        

        
    @property
    def pha(self):
        return self._binned_pha

    @property
    def time(self):
        return self._time

    @property
    def dead_time_fraction(self):
        return self._dead_time_fraction

    @property
    def rsp(self):
        return self._rsp

    @property
    def n_channels(self):

        return len(self._rsp.ebounds) - 1

    @property
    def scattering_angles(self):

        return self._binned_scattering_angles

    @property
    def scattering_angle_time(self):

        return self._scattering_angle_time

    @property
    def scattering_angle_dead_time_fraction(self):
        return self._scattering_angle_dead_time_fraction

    @property
    def n_scattering_bins(self):

        return len(self._scattering_bins) - 1
