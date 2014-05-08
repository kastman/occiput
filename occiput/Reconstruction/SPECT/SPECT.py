
import Scintillators 
import Collimators 
from mMR import UncompressedProjection 

from numpy import *
from numpy.random import randint 

# Import NiftyCore ray-tracers
from occiput.Core.NiftyCore_wrap import SPECT_project_parallelholes, SPECT_backproject_parallelholes, has_NiftyCore
from occiput.Core import Image3D
from occiput.Visualization import ProgressBar, svgwrite, has_svgwrite, ipy_table, has_ipy_table


DEFAULT_ITERATIONS  = 20
DEFAULT_SUBSET_SIZE = 32
EPS                 = 1e-9




class FileNotFound(Exception): 
    def __init__(self,msg,filename): 
        self.msg = str(msg) 
        self.filename = str(filename)
    def __str__(self): 
        return "Cannot find file '%s' (%s)."%(self.filename, self.msg)

class UnknownParameter(Exception): 
    def __init__(self,msg): 
        self.msg = str(msg) 
    def __str__(self): 
        return "Unkwnown parameter: %s"%(self.msg)

class UnexpectedParameter(Exception): 
    def __init__(self,msg): 
        self.msg = str(msg) 
    def __str__(self): 
        return "Unexpected parameter: %s"%(self.msg)





class SubsetGenerator():  
    def __init__(self,N_positions):
        self._N_positions = N_positions

    def new_subset(self,mode,subset_size): 
        if mode=='random': 
            return self._random_no_replacement(subset_size) 
        elif mode=='ordered': 
            raise UnexpectedParameter("'%s' subset selection mode not yet supported."%str(mode))
        else: 
            raise UnexpectedParameter("'mode' parameter %s not recognised."%str(mode))
        
    def all_active(self): 
        return ones((self._N_positions),dtype=uint32) 

    def _random_no_replacement(self,subset_size): 
        if subset_size>=self._N_positions: 
            return self.all_active() 
        M = zeros((self._N_positions),dtype=int32) 
        n = 0
        while n<subset_size:
            active = randint(self._N_positions)
            if M[active] == 0: 
                M[active] = 1
                n+=1
        return M




def deg_to_rad(deg):
    return deg*pi/180.0

def rad_to_deg(rad):
    return rad*180.0/pi
    






class SPECT_Static_Scan(object):
    def __init__(self): 
        self._name         = "Generic SPECT Scanner"     
        self._scanner_type = "SPECT" 
        self._manufacturer = "No manufacturer"
        self._version = "0.0"
        # scanner parameters are named with 'self._p_xxx' 
        self._p_gantry_angular_positions      = 180   #[adim,integer]
        self._p_gantry_angular_position_first = 0.0   #[degrees]
        self._p_gantry_angular_position_last  = 358.0 #[degrees]     
        self._subset_generator = SubsetGenerator(self._p_gantry_angular_positions)                     
        self._p_scan_time_sec = 600.0                 #[seconds]
        self._p_radius_mm     = 300.0                 #[mm]
        self._p_n_pix_x       = 128                   #[adim]
        self._p_n_pix_y       = 128                   #[adim]
        self._p_pix_size_x_mm = 2.5                   #[mm]
        self._p_pix_size_y_mm = 2.5                   #[mm]
        self.set_background_activity(0.0) 
        self.set_background_attenuation(0.0) 
        self.set_use_gpu(True) 
        self.set_truncate_negative(False)
        self.set_scintillator( Scintillators.Ideal() )
        self.set_collimator( Collimators.LEHR() ) 
        self._measurement = None 
        self._need_update_norm = True

    def get_name(self):
        return self._name

    def get_type(self):
        return self._scanner_type

    def get_manufacturer(self):
        return self._manufacturer

    def get_version(self): 
        return self._version

    def _get_parameters(self):
        parameters = {}
        dic = self.__dict__
        for k in dic.keys(): 
            if k.startswith('_p_'):
                parameters[k[3:]]=dic[k]        
        return parameters 

    def get_gantry_angular_positions(self): 
        return (self._p_gantry_angular_position_first, self._p_gantry_angular_position_last, self._p_gantry_angular_positions)

    def set_gantry_angular_positions(self, first_position_deg, last_position_deg, N_positions): 
        if not ( isscalar(first_position_deg) and isscalar(last_position_deg) and isscalar(N_positions) ): 
            raise UnexpectedParameter('Expected scalar values.')
        if not isinstance(N_positions,type(1)): 
            raise UnexpectedParameter('Expected an integer value.') 
        self._p_gantry_angular_position_first = first_position_deg
        self._p_gantry_angular_position_last =  last_position_deg
        self._p_gantry_angular_positions = N_positions
        self._subset_generator = SubsetGenerator(self._p_gantry_angular_positions)

    def get_scan_time(self): 
        return self._p_scan_time_sec

    def set_scan_time(self,scan_time_sec): 
        if not isscalar(scan_time_sec): 
            raise UnexpectedParameter('Expected a scalar value.')
        self._p_scan_time_sec = scan_time_sec

    def get_radius(self): 
        return self._p_radius_mm 

    def set_radius(self,radius_mm): 
        if not isscalar(radius_mm): 
            raise UnexpectedParameter('Expected a scalar value.') 
        self._p_radius_mm = radius_mm

    def get_n_pixels(self): 
        return (self._p_n_pix_x, self._p_n_pix_y)

    def set_n_pixels(self,n_pixels_x,n_pixels_y): 
        if (not isscalar(n_pixels_x)) or (not isscalar(n_pixels_y)): 
            raise UnexpectedParameter('Expected scalar values.') 
        self._p_n_pix_x = n_pixels_x
        self._p_n_pix_y = n_pixels_y
        self._need_update_norm = True

    def get_pixel_size(self): 
        return (self._p_pix_size_x, self._p_pix_size_y)

    def set_pixel_size(self,pixel_size_x,pixel_size_y): 
        if (not isscalar(pixel_size_x)) or (not isscalar(pixel_size_y)): 
            raise UnexpectedParameter('Expected scalar values.') 
        self._p_pix_size_x_mm = pixel_size_x
        self._p_pix_size_y_mm = pixel_size_y

    def get_scintillator(self):
        return self._scintillator 

    def set_scintillator(self,scintillator): 
        if not isinstance(scintillator,Scintillators.BaseScintillatorSPECT): 
            raise UnexpectedParameter('Expected an instance of BaseScintillatorSPECT')
        self._scintillator = scintillator 
        self.__make_psf() 
        self._need_update_norm = True

    def get_collimator(self): 
        return self._collimator 

    def set_collimator(self,collimator): 
        if not isinstance(collimator,Collimators.BaseCollimatorSPECT): 
            raise UnexpectedParameter('Expected an instance of BaseCollimatorSPECT')
        self._collimator = collimator 
        self.__make_psf() 
        self._need_update_norm = True
        
    def set_background_activity(self,value): 
        self._background_activity    = value 
        
    def get_background_activity(self,value): 
        return self._background_activity
    
    def set_background_attenuation(self,value): 
        self._background_attenuation = value 
        
    def get_background_attenuation(self,value): 
        return self._background_attenuation
        
    def set_use_gpu(self, value): 
        self._use_gpu = value 
    
    def set_truncate_negative(self,value): 
        self._truncate_negative = value 

    def project(self, activity, attenuation=None, cameras=None, psf=None, subsets_array=None): 
        if isinstance(activity,ndarray): 
            activity = float32(activity)
        else: 
            activity = float32(activity.data)
        if attenuation!=None:
            if isinstance(attenuation,ndarray):
                attenuation = float32(attenuation)
            else: 
                attenuation = float32(attenuation.data)
        if cameras==None: 
            cameras = float32(linspace(deg_to_rad(self._p_gantry_angular_position_first),deg_to_rad(self._p_gantry_angular_position_last),self._p_gantry_angular_positions).reshape((self._p_gantry_angular_positions,1))) 
        # subsets: 
        if subsets_array!=None: 
            cameras=cameras[where(subsets_array)]
        if psf==None: 
            psf=self._psf
        proj = SPECT_project_parallelholes(activity, cameras, attenuation, psf, self._background_activity, self._background_attenuation, self._use_gpu, self._truncate_negative)
        return UncompressedProjection(proj) 


    def backproject(self, projection, attenuation=None,  cameras=None, psf=None, subsets_array=None):
        if isinstance(projection,ndarray): 
            projection = float32(projection)
        else: 
            projection = float32(projection.data)
        if attenuation!=None:
            if isinstance(attenuation,ndarray):
                attenuation = float32(attenuation)
            else: 
                attenuation = float32(attenuation.data)
        if cameras==None: 
            cameras = float32(linspace(deg_to_rad(self._p_gantry_angular_position_first),deg_to_rad(self._p_gantry_angular_position_last),self._p_gantry_angular_positions).reshape((self._p_gantry_angular_positions,1)))
        # subsets: 
        if subsets_array!=None: 
            cameras=cameras[where(subsets_array)]
        if psf==None: 
            psf=self._psf   
        backproj = SPECT_backproject_parallelholes(projection, cameras, attenuation, psf, self._background_activity, self._background_attenuation, self._use_gpu, self._truncate_negative)
        return Image3D(backproj)

    def scan(self,activity_Bq,scan_time_sec=None): 
        if scan_time_sec==None: 
            scan_time_sec = self.get_scan_time() 
        sinogram = 0
        return sinogram
        
    def __make_probabilistic_graphical_model(self): 
        pass 

    def __make_psf(self): 
        self._psf = None

    def get_normalization(self): 
        if self._need_update_norm: 
            self._compute_normalisation() 
        return self._norm 
        
    def _compute_normalisation(self): 
        subsets_array = self._subset_generator.all_active()
        self._norm = self.backproject(ones(( self._p_n_pix_x,self._p_n_pix_y,self._p_gantry_angular_positions ),dtype=float32, order="F") ).data 
        self._need_update_norm = False 

    def estimate_activity(self, iterations=DEFAULT_ITERATIONS, subset_size=DEFAULT_SUBSET_SIZE, subset_mode='random'): 
        progress_bar = ProgressBar() 
        progress_bar.set_percentage(0.1)
        activity = ones((self._p_n_pix_x,self._p_n_pix_y,self._p_n_pix_x),dtype=float32, order="F")
        for i in range(iterations): 
            # Subsets: 
            if subset_size==None:
                subsets_array=None
                subset_size=self._p_gantry_angular_positions
            elif subset_size>=self._p_gantry_angular_positions: 
                subsets_array=None
                subset_size=self._p_gantry_angular_positions
            else:
                subsets_array = self._subset_generator.new_subset(subset_mode,subset_size)
            if subsets_array!=None: 
                proj = self.project(activity,subsets_array=subsets_array).data
                P = (self._measurement[:,:,where(subsets_array)].reshape((self._p_n_pix_x,self._p_n_pix_y,subset_size))+EPS)/(proj+EPS)
                norm = self.backproject(ones(( self._p_n_pix_x,self._p_n_pix_y,subset_size ),dtype=float32, order="F"), subsets_array=subsets_array).data 
                update = (self.backproject( P ,subsets_array=subsets_array).data+EPS) / (norm +EPS) 
            else: 
                proj = self.project(activity).data
                P = (self._measurement+EPS)/(proj+EPS)   
                norm = self.get_normalization()  
                update = (self.backproject( P ).data+EPS) / (norm +EPS) 
            activity = activity * update #* self.get_mask().data

            progress_bar.set_percentage((i+1)*100.0/iterations) 
            #print "Iteration: %d    max act: %f    min act: %f    max proj: %f    min proj: %f    max norm: %f    min norm: %f"%(i, activity.max(), activity.min(), proj.max(), proj.min(), norm.data.max(), norm.data.min() )
        progress_bar.set_percentage(100.0)
        return Image3D(activity)
            
    def volume_render(self,volume,scale=1.0): 
        # FIXME: use the VolumeRenderer object in occiput.Visualization (improve it), the following is a quick fix: 
        if isinstance(volume,ndarray): 
            volume = float32(volume)
        else: 
            volume = float32(volume.data)
        proj = self.project(volume).data
        proj[where(proj>proj.max()/scale )]=proj.max()/scale
        return UncompressedProjection(proj)

    def load_measurement_file(self,filename): 
        pass 

    def set_measurement(self, measurement): 
        if not ( self._p_n_pix_x == measurement.shape[0] and self._p_n_pix_y == measurement.shape[1] and self._p_gantry_angular_positions == measurement.shape[2] ):
            raise UnexpectedParameter('Measurement size is not compatible with n_pix_x, n_pix_y, gantry_angular_positions. ')
        self._measurement = measurement 

    def get_measurement(self): 
        return Volume(self._measurement)

    def display_measurement(self): 
        return UncompressedProjection(self._measurement)

    def _make_svg(self): 
        if not has_svgwrite: 
            self._svg_string = None 
            return self._svg_string 

        w = '100%'
        h = '100%'
        
        dwg = svgwrite.Drawing('SPECT.svg',size=(w,h), profile='full', debug=True)
        dwg.viewbox(width=100, height=100)

        # DETECTOR 
        # collimator 
        rect = dwg.add(dwg.rect(insert=(12, 30), size=(8, 40), rx=0.5, ry=0.5))
        rect.fill('grey',opacity=0.5).stroke('black',width=0.3,opacity=0.001)

        # scintillator
        rect = dwg.add(dwg.rect(insert=(9, 30), size=(3, 40), rx=0.5, ry=0.5))
        rect.fill('green',opacity=0.1).stroke('none',width=0.3,opacity=0.001)

        # photomultipliers 
        for i in range(8): 
            rect = dwg.add(dwg.rect(insert=(1, 31.2+i*4.8), size=(8, 4), rx=0.3, ry=0.3))
            rect.fill('grey',opacity=0.25).stroke('none',width=0.3,opacity=0.001)
        
        # IMAGING VOLUME
        rect = dwg.add(dwg.rect(insert=(30, 30), size=(40, 40), rx=0.5, ry=0.5))
        rect.fill('grey',opacity=0.02).stroke('grey',width=0.3,opacity=0.02)
        
        # GEOMETRIC NOTATIONS 
        # circle, gantry rotation 
        circle = dwg.add(dwg.circle(center=(50, 50), r=30))
        circle.fill('none').stroke('grey', width=0.1).dasharray([0.5, 0.5]) 
        # center
        circle = dwg.add(dwg.circle(center=(50, 50), r=0.5))
        circle.fill('grey',opacity=0.1).stroke('grey', width=0.1)    
        line = dwg.add(dwg.line(start=(50-1,50), end=(50+1,50)))
        line.stroke('grey', width=0.1) 
        line = dwg.add(dwg.line(start=(50,50-1), end=(50,50+1)))
        line.stroke('grey', width=0.1) 
        
        #line = dwg.add(dwg.polyline([(10, 10), (10, 100), (100, 100), (100, 10), (10, 10)],stroke='black', fill='none'))
        self._svg_string = dwg.tostring() 
        return self._svg_string

    def _repr_svg_(self): 
        self._make_svg()
        return self._svg_string    





class Gantry(): 
    def __init__(self): 
        self.svg_string = self.make_svg()

    def make_svg(self):
        if not has_svgwrite: 
            self._svg_string = None 
            return self._svg_string 

        w = '100%'
        h = '100%'
        
        dwg = svgwrite.Drawing('test.svg',size=(w,h), profile='full', debug=True)
        dwg.viewbox(width=100, height=100)

        # DETECTOR 
        # collimator 
        rect = dwg.add(dwg.rect(insert=(12, 30), size=(8, 40), rx=0.5, ry=0.5))
        rect.fill('grey',opacity=0.5).stroke('black',width=0.3,opacity=0.001)

        # scintillator
        rect = dwg.add(dwg.rect(insert=(9, 30), size=(3, 40), rx=0.5, ry=0.5))
        rect.fill('green',opacity=0.1).stroke('none',width=0.3,opacity=0.001)

        # photomultipliers 
        for i in range(8): 
            rect = dwg.add(dwg.rect(insert=(1, 31.2+i*4.8), size=(8, 4), rx=0.3, ry=0.3))
            rect.fill('grey',opacity=0.25).stroke('none',width=0.3,opacity=0.001)
        
        # IMAGING VOLUME
        rect = dwg.add(dwg.rect(insert=(30, 30), size=(40, 40), rx=0.5, ry=0.5))
        rect.fill('grey',opacity=0.02).stroke('grey',width=0.3,opacity=0.02)
        
        # GEOMETRIC NOTATIONS 
        # circle, gantry rotation 
        circle = dwg.add(dwg.circle(center=(50, 50), r=30))
        circle.fill('none').stroke('grey', width=0.1).dasharray([0.5, 0.5]) 
        # center
        circle = dwg.add(dwg.circle(center=(50, 50), r=0.5))
        circle.fill('grey',opacity=0.1).stroke('grey', width=0.1)    
        line = dwg.add(dwg.line(start=(50-1,50), end=(50+1,50)))
        line.stroke('grey', width=0.1) 
        line = dwg.add(dwg.line(start=(50,50-1), end=(50,50+1)))
        line.stroke('grey', width=0.1) 
        
        #line = dwg.add(dwg.polyline([(10, 10), (10, 100), (100, 100), (100, 10), (10, 10)],stroke='black', fill='none'))
        return dwg.tostring() 

    def _repr_svg_(self): 
        return self.svg_string 




class GE_Infinia(SPECT_Static_Scan):
    def __init__(self): 
        SPECT.__init__(self)
        self._name = "GE Infinia SPECT Scanner with LEHR collimator"  






