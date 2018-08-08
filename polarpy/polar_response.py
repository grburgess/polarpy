import numpy as np
import h5py
import scipy.interpolate as interpolate



class PolarResponse(object):

    def __init__(self, response_file):

        self._rsp_file = response_file

        self._interpolate_rsp()

    # def _get_hist(self, ene,degree,angle):


    #     with h5py.File(self._rsp_file,'r') as f:



    #        # bins = np.array(f['bins'].value)

    #     return  tmp


    def _interpolate_rsp(self):
        """
        Builds the interpolator for the response. This is currently incredibly slow
        and should be improved

        """

        
        # create a functions to get the histograms
        
        
        # now go through the response and extract things
        
        with h5py.File(self._rsp_file,'r') as f:
        
            ene = [int(str(x[4:])) for x in f.keys() if 'ene'  in x]
            
            
            energy = np.sort(np.array(ene))
            
            ene_lo, ene_hi = [],[]

            for ene in energy:

                ene_lo.append(ene-2.5)
                ene_hi.append(ene+2.5)
            
            
            pol_ang = np.array(f['pol_ang'].value)
        
            pol_deg = np.array(f['pol_deg'].value)
        
        
            bins = np.array(f['bins'].value)

            
            bin_center = 0.5 *(bins[:-1] + bins[1:])

            all_interp = []



            for i, bm in enumerate(bin_center):

                data = []
                #energy = get_energy()

                for ene in energy:

                    for ang in pol_ang:

                        for deg in pol_deg:

                            tmp = np.array(f['ene_%d' % int(ene)]['deg_%d' % int(deg)]['ang_%d' % int(ang)].value)

                            hist = self._get_hist(ene,deg,ang)


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
            self._scattering_bins_lo = bins[:-1]
            self._scattering_bins_hi = bins[1:]


    @property
    def ene_lo(self):
        return self._ene_lo

    @property
    def ene_hi(self):
        return self._ene_hi

    @property
    def energy_mid(self):
        return self._energy_mid

    @property
    def n_scattering_bins(self):
        return self._n_scatter_bins

    @property
    def scattering_bins(self):
        return self._scattering_bins

    @property
    def scattering_bins_lo(self):
        return self._scattering_bins_lo

    @property
    def scattering_bins_hi(self):
        return self._scattering_bins_hi

    @property
    def interpolators(self):
        return self._all_interp
        
