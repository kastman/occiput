
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA



import nibabel
from occiput.Core.Conversion import nipy_to_occiput, occiput_from_array
from occiput.DataSources.FileSources.LookupTable import load_freesurfer_lut_file
import os
import numpy
import dicom


import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import nipy




def load_image_file(filename): 
    nip = nipy.load_image(filename)
    return nipy_to_occiput(nip)
    
    
def load_mask_file(filename, lookup_table_filename=None): 
    # Load file 
    nip = nipy.load_image(filename)
    occ = nipy_to_occiput(nip)
    occ.set_mask_flag(1) 
    
    # Load the lookup table. If not specified, try to load from file with the same name as 
    # the mask image file. 
    if lookup_table_filename == None: 
        f = []
        f.append(os.path.splitext(filename)[0]+'.lut')
        f.append(os.path.splitext(os.path.splitext(filename)[0])[0]+'.lut') #This includes .nii.gz files 
        for lookup_table_filename in f: 
            try: 
                lut = load_freesurfer_lut_file(lookup_table_filename) 
            except: 
                lut = None 
                
    else: 
        lut = load_freesurfer_lut_file(lookup_table_filename) 
    if lut != None: 
        occ.set_lookup_table(lut) 
    return occ


def load_dicom_series(path, files_start_with=None, files_end_with=None, exclude_files_end_with=['.dat','.txt','.py','.pyc','.nii','.gz'] ):
        """Rudimentary file to load dicom serie from a directory. """ 
        N=0 
        paths  = []
        slices = []
        files = os.listdir(path)
 
        for file_name in files: 
            file_valid = True
            if files_start_with!=None: 
                if not file_name.startswith(files_start_with): 
                    file_valid = False
            if files_end_with!=None: 
                if not file_name.endswith(files_end_with): 
                    file_valid = False 
            for s in exclude_files_end_with: 
                if file_name.endswith(s): 
                    file_valid = False           
            if file_valid: 
                full_path = path+os.sep+file_name
                # read moco information from files 
                paths.append(full_path)
                f = dicom.read_file(full_path) 
                slice = f.pixel_array
                slices.append(slice)
                N+=1
                instance_number    = f.get(0x00200013).value 
                creation_time      = f.get(0x00080013).value
                #print "Instance number:    ",instance_number
                #print "Creation time:      ",creation_time 
        array = numpy.zeros((slices[0].shape[0],slices[0].shape[1],N),dtype=numpy.float32)
        for i in range(N):  
            slice = numpy.float32(slices[i])  #FIXME: handle other data types 
            array[:,:,i] = slice 
        #return occiput_from_array(array)
        return array
    
