
from DisplayNode import DisplayNode
import Image
import numpy



__all__ = ['MultipleVolumes','VolumeRenderer','ProgressBar']


import uuid
from IPython.display import HTML, Javascript, display

    
class ProgressBar(): 
    def __init__(self): 
        self.divid = str(uuid.uuid4())
        self.pb = HTML(
        """
        <div style="border: 1px solid white; width:800px; height:6px; background-color:rgb(200,228,246)">
            <div id="%s" style="background-color:rgb(47,128,246); width:0%%; height:6px">&nbsp;</div>
        </div> 
        """ % self.divid)
        self.visible = False
    def show(self):
        display(self.pb)
        self.visible = True 
        
    def set_percentage(self,percentage):
        if not self.visible: 
            self.show()
        if percentage < 1:
            percentage = 1
        if percentage > 100:
            percentage = 100 
        display(Javascript("$('div#%s').width('%i%%')" % (self.divid, percentage)))
        if percentage >= 100:
            self.set_done()

    def set_done(self):
        #display(Javascript("$('div#%s').css({'background-color':'rgb(208,92,92)'});"%self.divid))
        pass 







class MultipleVolumes(): 
    def __init__(self,volumes,axis=0,shrink=256,rotate=90,scale_factors=None): 
        self.volumes = volumes 
        self._axis = axis
        self._shrink = shrink
        self._rotate = rotate
        self._scale_factors = scale_factors
        
    def get_data(self,volume_index): 
        volume = self.volumes[volume_index] 
        if isinstance(volume,numpy.ndarray): 
            return volume
        else:
            # Volume type
            return volume.data
        
    def get_shape(self,volume_index): 
        return self.volumes[volume_index].shape
        
    def to_image(self,volume_index,slice_index,axis=0,normalise=True,scale_factor=None,shrink=None,rotate=0): 
        #FIXME: handle 4D, 5D, ..
        if axis == 0:
            a = numpy.float32(self.get_data(volume_index)[slice_index,:,:].reshape((self.get_shape(volume_index)[1],self.get_shape(volume_index)[2]))) 
        elif axis == 1: 
            a = numpy.float32(self.get_data(volume_index)[:,slice_index,:].reshape((self.get_shape(volume_index)[0],self.get_shape(volume_index)[2])))
        else: 
            a = numpy.float32(self.get_data(volume_index)[:,:,slice_index].reshape((self.get_shape(volume_index)[0],self.get_shape(volume_index)[1])))
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
    def display_in_browser(self,axis=None,shrink=None,rotate=None,scale_factors=None): 
        self.display(axis,shrink,rotate,scale_factors,open_browser=True)

    def display(self,axis=None,shrink=None,rotate=None,scale_factors=None,open_browser=False):
        if axis==None:
            axis = self._axis
        if shrink == None:
            shrink = self._shrink
        if rotate == None: 
            rotate = self._rotate
        if scale_factors == None: 
            scale_factors = self._scale_factors  
        D = DisplayNode()
        images = []
        n=0
        
        #N_slices = 0
        #for j in range(len(self.volumes)): 
        #    N_slices += self.get_shape(j)[axis] 
        #bar = ProgressBar() 
        
        for j in range(len(self.volumes)): 
            if scale_factors == None: 
                scale_factor = 255/(self.get_data(j).max()+1e-9)
            else: 
                scale_factor = scale_factors[j]
            images_inner = []
            for i in range(self.get_shape(j)[axis]): 
                im = self.to_image(j,i,axis,normalise=True,scale_factor=scale_factor,shrink=shrink,rotate=rotate)
                images_inner.append(im) 
                #n+=1
                #progress_bar.set_percentage(n*100.0/N_slices) 
            images.append(images_inner) 
        return D.display('tipix', images, open_browser)   

    def _repr_html_(self): 
        return self.display()._repr_html_()
        




class VolumeRenderer(): 
    def __init__(self,volume): 
        self.volume = volume 

    def render(self): 
        pass 
    
    
    




