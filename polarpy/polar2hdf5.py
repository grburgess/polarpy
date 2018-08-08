import h5py
import numpy as np


try:
    import ROOT
    from threeML.io.cern_root_utils.io_utils import open_ROOT_file
    from threeML.io.cern_root_utils.tobject_to_numpy import tree_to_ndarray, th2_to_arrays

    has_root = True

except(ImportError):

    has_root = False

    

if has_root:

    def polar_polarization_to_hdf5(polarization_root_file, hdf5_out_file):
        """

        :param polarization_root_file: The ROOT file from which to build the response
        :param hdf5_out_file: The output HDF5 file name
        """
        energy = []
        degree = []
        angle = []

        energy_str = []
        degree_str = []
        angle_str = []

        
        with open_ROOT_file(polarization_root_file) as f:

            tmp = [key.GetName() for key in f.GetListOfKeys()]
            tmp = filter(lambda x: 'sim' in x, tmp)
            for tmp2 in tmp:
                _, x, y, z = tmp2.split('_')

                energy.append(float(x))
                degree.append(float(y))
                angle.append(float(z))

                energy_str.append(x)
                degree_str.append(y)
                angle_str.append(z)


                

            energy = np.array(np.unique(energy))
            degree = np.array(np.unique(degree))
            angle  = np.array(np.unique(angle))

            energy_str = np.array(np.unique(energy_str))
            degree_str = np.array(np.unique(degree_str))
            angle_str  = np.array(np.unique(angle_str))

            # just to get the bins
            # must change this from ints later

            file_string = 'sim_%s_%s_%s' % (energy_str[1], degree_str[1], angle_str[1])

            bins, _, hist = th2_to_arrays(f.Get(file_string))

            out_matrix = np.zeros((len(energy), len(degree), len(angle), len(hist)))

            with h5py.File(hdf5_out_file, 'w', libver='latest') as database:

                for i,x in enumerate(energy_str):



                    for j, y in enumerate(degree_str):



                        for k, z in enumerate(angle_str):

                            file_string = 'sim_%s_%s_%s' % (x, y, z)

                            _ , _, hist = th2_to_arrays(f.Get(file_string))

                            out_matrix[i,j,k,:] = hist



                database.create_dataset('matrix',data=out_matrix,compression='lzf')


                if np.min(bins) < 0:
                    # we will try to automatically correct for the badly specified bins
                    bins = np.array(bins)

                    bins += -np.min(bins)

                    assert np.min(bins) >=0, 'The scattering bins have egdes less than zero'
                    assert np.max(bins) <=360, 'The scattering bins have egdes greater than 360'

                    
                
                database.create_dataset('bins', data=bins, compression='lzf')
                database.create_dataset('pol_ang', data=angle, compression='lzf')
                database.create_dataset('pol_deg', data=degree, compression='lzf')
                database.create_dataset('energy',data=energy,compression='lzf')

    def polar_spectra_to_hdf5(polar_root_file, polar_rsp_root_file, hdf5_out_file):
        """

        :param polar_root_file: The spectral ROOT file
        :param polar_rsp_root_file: The response ROOT file
        :param hdf5_out_file: the name of the output HDF5 file
        """

        # extract the info from the crappy root file

        with h5py.File(hdf5_out_file, 'w') as outfile:

            # first we do the RSP

            rsp_grp = outfile.create_group('rsp')

            with open_ROOT_file(polar_rsp_root_file) as f:

                matrix = th2_to_arrays(f.Get('rsp'))[-1]

                rsp_grp.create_dataset('matrix', data=matrix, compression='lzf')

                ebounds = th2_to_arrays(f.Get('EM_bounds'))[-1]

                rsp_grp.create_dataset('ebounds', data=ebounds, compression='lzf')

                mc_low = th2_to_arrays(f.Get('ER_low'))[-1]

                rsp_grp.create_dataset('mc_low', data=mc_low, compression='lzf')

                mc_high = th2_to_arrays(f.Get('ER_high'))[-1]

                rsp_grp.create_dataset('mc_high', data=mc_high, compression='lzf')

            # now we get the spectral informations
            keys_to_use = ['polar_out']

            f = ROOT.TFile(polar_root_file)

            extra_grp = outfile.create_group('extras')

            for key in f.GetListOfKeys():

                name = key.GetName()

                if name not in keys_to_use:
                    try:

                        # first we see if it is a TTree and then
                        # add a new group and attach its data

                        tree = tree_to_ndarray(f.Get(name))

                        new_grp = extra_grp.create_group(name)

                        for new_name in tree.dtype.names:
                            new_grp.create_dataset(new_name, data=tree[new_name], compression='lzf')

                    except:

                        # in this case we just want the actual data

                        data = th2_to_arrays(f.Get(name))[-1]

                        extra_grp.create_dataset(name, data=data, compression='lzf')

            # now we will deal with the data that is important

            tmp = tree_to_ndarray(f.Get('polar_out'))

            outfile.create_dataset('energy', data=tmp['Energy'], compression='lzf')

            outfile.create_dataset('scatter_angle', data=tmp['scatter_angle'], compression='lzf')
            
            outfile.create_dataset('dead_ratio', data=tmp['dead_ratio'], compression='lzf')

            outfile.create_dataset('time', data=tmp['tunix'], compression='lzf')

            
            
            f.Close()
