from threeML.io.cern_root_utils.io_utils import open_ROOT_file
from threeML.io.cern_root_utils.tobject_to_numpy import tree_to_ndarray, th2_to_arrays
import ROOT
import h5py

def polar_root2hdf5(polar_root_file, polar_rsp_root_file ,hdf5_out_file):
    """
    
    :param polar_root_file: The spectral ROOT file
    :param polar_rsp_root_file: The response ROOT file
    :param hdf5_out_file: the name of the output HDF5 file
    """


    # extract the info from the crappy root file

    with h5py.File(hdf5_out_file,'w') as outfile:


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


        f=ROOT.TFile(polar_root_file)


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
        
        outfile.create_dataset('energy',data=tmp['Energy'],compression='lzf')

        outfile.create_dataset('dead_ratio',data=tmp['dead_ratio'],compression='lzf')

        outfile.create_dataset('time',data=tmp['tunix'],compression='lzf')

        f.Close()
