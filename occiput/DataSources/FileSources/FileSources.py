

# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA

__all__ = ['FileSource']


import nibabel
import Image
from DisplayNode import DisplayNode
import numpy 

class FileSource(object):
    def __init__(self,filename=None):
        self.nifti_image=None
        if filename != None: 
            self.load(filename)

    def load(self,file): 
        self.nifti_image = nibabel.load(file) 

    def get_data(self): 
        return self.nifti_image.get_data() 
 
    def get_affine(self):
        return self.nifti_image.get_affine()
        
    def get_info(self): 
        return self.nifti_image.get_info() 

    def save(self,file):
        return self.nifti_image.save(file)  
 
    def get_data_dtype(self):
        return self.nifti_image.get_data_type() 
         
    def get_filename(self):
        return self.nifti_image.get_filename() 
         
    def get_header(self):
        return self.nifti_image.get_header()  
        
    def get_qform(self):
        return self.nifti_image.get_qform()  
        
    def get_sform(self):
        return self.nifti_image.get_sform()  
        
    def get_shape(self):
        return self.nifti_image.get_shape()  

    def set_filename(self,filename):
        return self.nifti_image.set_filename(filename) 
        
    def set_qform(self,q):
        return self.nifti_image.set_qform(q)
        
    def set_sform(self,s):
        return self.nifti_image.set_sform(s)
        
    def is_instance_to_filename(self):
        return self.nifti_image.is_instance_to_filename()

    def to_image(self,index,axis=0,normalise=True,scale_factor=None,shrink=None,rotate=0): 
        #FIXME: handle 4D, 5D, ..
        if axis == 0:
            a = numpy.float32(self.get_data()[index,:,:].reshape((self.get_shape()[1],self.get_shape()[2]))) 
        elif axis == 1: 
            a = numpy.float32(self.get_data()[:,index,:].reshape((self.get_shape()[0],self.get_shape()[2])))
        else: 
            a = numpy.float32(self.get_data()[:,:,index].reshape((self.get_shape()[0],self.get_shape()[1])))
        if normalise: 
            if scale_factor==None:
                a = a*255/(a.max()+1e-9)
            else: 
                a = a * scale_factor 
        im =  Image.fromarray(a).convert("RGB") 
        if shrink != None: 
            # scale in order to make the largest dimension equal to 'shrink' 
            shrink = int(shrink) 
            (h,w) = im.size
#FIXME: downsample the volume (with the GPU, if available) before converting to images, it will save a lot of time, conversion to Image and to RGB is slow
            if h>shrink or w>shrink: 
                if h>w: 
                    im = im.resize( (shrink, int(shrink*w/(1.0*h))) )
                else: 
                    im = im.resize( (int(shrink*h/(1.0*w)),shrink) )
            if rotate: 
                im = im.rotate(rotate)
        return im
        

    # display
    def display_in_browser(self,axis=0,shrink=256,rotate=90,scale_factor=None): 
        if scale_factor == None: 
            scale_factor = 255/(self.get_data().max()+1e-9)
        D = DisplayNode()
        images = []
        for i in range(self.get_shape()[axis]): 
            im = self.to_image(i,axis,normalise=True,scale_factor=scale_factor,shrink=shrink,rotate=rotate)
            images.append(im) 
        D.display_in_browser('tipix', images)  
        
    def display(self,axis=0,shrink=256,rotate=90,scale_factor=None): 
        if scale_factor == None: 
            scale_factor = 255/(self.get_data().max()+1e-9)
        D = DisplayNode()
        images = []
        for i in range(self.get_shape()[axis]): 
            im = self.to_image(i,axis,normalise=True,scale_factor=scale_factor,shrink=shrink,rotate=rotate)
            images.append(im) 
        return D.display('tipix', images)   

    def _repr_html_(self): 
        return self.display()._repr_html_()



