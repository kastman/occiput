
from DisplayNode import DisplayNode
import Image
import numpy



__all__ = ['MultipleVolumes','VolumeRenderer','ProgressBar','LIGHT_BLUE','LIGHT_RED','BLUE','RED','GRAY','LIGHT_GRAY']

LIGHT_BLUE   = 'rgb(200,228,246)'
BLUE         = 'rgb(47,128,246)'
LIGHT_RED    = 'rgb(246,228,200)'
RED          = 'rgb(246,128,47)'
LIGHT_GRAY   = 'rgb(246,246,246)'
GRAY         = 'rgb(200,200,200)'


import uuid
from IPython.display import HTML, Javascript, display

    
class ProgressBar(): 
    def __init__(self, height='6px', width='100%%', background_color=LIGHT_BLUE, foreground_color=BLUE): 
        self.divid = str(uuid.uuid4())
        self.pb = HTML(
        """
        <div style="border: 1px solid white; width:%s; height:%s; background-color:%s">
            <div id="%s" style="background-color:%s; width:0%%; height:%s"> </div>
        </div> 
        """ % ( width, height, background_color, self.divid, foreground_color, height))
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







class MultipleVolumes(): 
    def __init__(self,volumes,axis=0,shrink=256,rotate=90,subsample_slices=None,scales=None,open_browser=None): 
        self.volumes = volumes 
        self._axis = axis
        self._shrink = shrink
        self._rotate = rotate
        self._subsample_slices = subsample_slices
        self._scales = scales
        self._open_browser = open_browser
        self._progress_bar = ProgressBar(height='6px', width='100%%', background_color=LIGHT_GRAY, foreground_color=GRAY)

    def get_data(self,volume_index): 
        volume = self.volumes[volume_index] 
        if isinstance(volume,numpy.ndarray): 
            return volume
        else:
            # Volume type
            return volume.data
        
    def get_shape(self,volume_index): 
        return self.volumes[volume_index].shape
        
    def to_image(self,volume_index,slice_index,axis=0,normalise=True,scale=None,shrink=None,rotate=0,global_scale=True): 
        #FIXME: handle 4D, 5D, ..
        
        M = self.get_data(volume_index).max()
        m = self.get_data(volume_index).min()
        
        if axis == 0:
            a = numpy.float32(self.get_data(volume_index)[slice_index,:,:].reshape((self.get_shape(volume_index)[1],self.get_shape(volume_index)[2]))) 
        elif axis == 1: 
            a = numpy.float32(self.get_data(volume_index)[:,slice_index,:].reshape((self.get_shape(volume_index)[0],self.get_shape(volume_index)[2])))
        else: 
            a = numpy.float32(self.get_data(volume_index)[:,:,slice_index].reshape((self.get_shape(volume_index)[0],self.get_shape(volume_index)[1])))

        if m>=0: 
            if normalise: 
                if not global_scale: 
                    if scale==None:
                        a = a * 255/(a.max()+1e-9) 
                    else: 
                        a = a * scale *255/(a.max()+1e-9) 
                else: 
                    if scale==None:
                        a = a * 255/(M+1e-9) 
                    else: 
                        a = a * scale *255/(M+1e-9) 
            im = Image.fromarray(a).convert("RGB") 
        else:
            if normalise: 
                if not global_scale:
                    if scale==None:
                        a = a * 512/(a.max()-a.min()+1e-9) 
                    else: 
                        a = a * scale * 512/(a.max()-a.min()+1e-9) 
                else: 
                    if scale==None:
                        a = a * 512/(M-m+1e-9) 
                    else: 
                        a = a * scale * 512/(M-m+1e-9)                     
            blue = a
            red  = a.copy()
            red[numpy.where(red>=0)]  = 0 
            red  = - red
            blue[numpy.where(blue<0)] = 0 
            green = numpy.zeros(red.shape)
            rgb = numpy.zeros((red.shape[0],red.shape[1],3),dtype=numpy.uint8)
            rgb[:,:,0]=red
            rgb[:,:,1]=green
            rgb[:,:,2]=blue
            im = Image.fromarray(rgb,mode='RGB')
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
    def display_in_browser(self,axis=None,shrink=None,rotate=None,subsample_slices=None,scales=None): 
        self.display(axis,shrink,rotate,subsample_slices,scales,open_browser=True)

    def display(self,axis=None,shrink=None,rotate=None,subsample_slices=None,scales=None,open_browser=None):
        if axis==None:
            axis = self._axis
        if shrink == None:
            shrink = self._shrink
        if rotate == None: 
            rotate = self._rotate
        if  subsample_slices == None: 
            subsample_slices = self._subsample_slices
        if subsample_slices==None: 
            subsample_slices=1
        if scales == None: 
            scales = self._scales 
        if open_browser == None: 
            open_browser = self._open_browser
        if open_browser == None: 
            open_browser = False 
        D = DisplayNode()
        images = []
        n=0
        
        N_slices = 0
        for j in range(len(self.volumes)): 
            N_slices += self.get_shape(j)[axis] 
        self._progress_bar = ProgressBar(height='6px', width='100%%', background_color=LIGHT_GRAY, foreground_color=GRAY)
        
        for j in range(len(self.volumes)): 
            if scales == None: 
                scale = 255/(self.get_data(j).max()+1e-9)
            else: 
                scale = scales[j]*255/(self.get_data(j).max()+1e-9)
            images_inner = []
            for i in range(0,self.get_shape(j)[axis],subsample_slices): 
                im = self.to_image(j,i,axis,normalise=True,scale=scale,shrink=shrink,rotate=rotate)
                images_inner.append(im) 
                n+=1
                self._progress_bar.set_percentage(n*100.0/N_slices*subsample_slices) 
            images.append(images_inner) 
        return D.display('tipix', images, open_browser)   

    def _repr_html_(self): 
        return self.display()._repr_html_()
        




class VolumeRenderer(): 
    def __init__(self,volume): 
        self.volume = volume 

    def render(self): 
        pass 
    
    
    




