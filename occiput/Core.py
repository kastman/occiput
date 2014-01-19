
# occiput  
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# Dec. 2013, Boston, MA


__all__ = ['Volume','RigidVolume']


import DisplayNode 
import numpy
import Image
from occiput.Visualization import ProgressBar


XYZ_ROTATION = 0



class Volume(object):
    pass 


class RigidVolume(Volume): 
    def __init__(self,data=None, scale=None, translation=None, rotation=None, rotation_order=None):
        if data==None:  
            self.data = numpy.asarray([]) 
        else: 
            # FIXME: check type, raise error
            self.data = data 
        if scale==None: 
            self.set_scale(numpy.asarray([1.0, 1.0, 1.0]))
        else: 
            self.set_scale(scale)
        if translation==None:
            self.set_translation(numpy.asarray([0.0, 0.0, 0.0]))
        else: 
            self.set_translation(translation)
        if rotation==None:
            self.set_rotation(numpy.asarray([0.0, 0.0, 0.0]))
        else: 
            self.set_rotation(rotation)
        if rotation_order==None:
            self.set_rotation_order(XYZ_ROTATION)
        else: 
            self.set_rotation_order(rotation_order)
        

    def get_shape(self): 
        return self.data.shape 
    
    def reshape(self,shape): 
        self.data.respahe(shape)

    def get_scale(self): 
        return self._scale
        
    def set_scale(self,scale): 
        self._scale = scale 
    
    def get_translation(self): 
        return self._translation 
    
    def set_translation(self,tr):
        self._translation = tr 
    
    def get_rotation(self):
        return self._rotation 
        
    def set_rotation(self,rot):
        self._rotation = rot 
    
    def get_rotation_order(self):
        return self._rotation_order
    
    def set_rotation_order(self,order):
        self._rotation_order = order 
    
    def get_sform_matrix(self): 
        sform = 0            # FIXME 
        return sform 
        
    def set_sform_matrix(self,sform): 
        self._translation = 0 #FIXME: allow only rigid transformation matrix, compute translation, rotation
        self._rotation    = 0 

    def get_data(self): 
        return self.data 
    
    def set_data(self,data):
        self.data = data 


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

    def display_in_browser(self,axis=0,shrink=256,rotate=90,scale_factor=None): 
        self.display(axis,shrink,rotate,scale_factor,open_browser=True) 
        
    def display(self,axis=0,shrink=256,rotate=90,scale_factor=None,open_browser=False): 
        if scale_factor == None: 
            scale_factor = 255/(self.get_data().max()+1e-9)
        D = DisplayNode.DisplayNode()
        images = []
        #bar = ProgressBar()
        
        for i in range(self.get_shape()[axis]): 
            im = self.to_image(i,axis,normalise=True,scale_factor=scale_factor,shrink=shrink,rotate=rotate)
            images.append(im) 
            #bar.set_percentage(i*100.0/(self.get_shape()[axis]))
        return D.display('tipix', images,open_browser)   

    def _repr_html_(self): 
        return self.display()._repr_html_()

    shape          = property(get_shape, reshape, None,"Shape of the volume")
    scale          = property(get_scale, set_scale, None,"Scale of the volume")
    translation    = property(get_translation, set_translation, None,"Translation of the volume")
    rotation       = property(get_rotation, set_rotation, None,"Rotation of the volume")
    rotation_order = property(get_rotation_order, set_rotation_order, None,"Order of rotation") 
    sform          = property(get_sform_matrix, set_sform_matrix, None, "Affine transformation matrix")
