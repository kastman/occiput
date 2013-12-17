
import Scintillators
import Collimators





class BaseScanner(): 
    def __init__(self): 
        self._name = "Generic Scanner" 
        self._scanner_type = "Generic scanner"
        self._manufacturer = "No manufacturer"
        self._version = "0.0"

    def get_name(self):
        return self._name

    def get_type(self):
        return self._scanner_type

    def get_manufacturer(self):
        return self._manufacturer

    def get_version(self): 
        return self._version

    def get_parameters(self):
        parameters = {}
        dic = self.__dict__
        for k in dic.keys(): 
            if k.startswith('_p_'):
                parameters[k[3:]]=dic[k]        
        return parameters 




class SPECT(BaseScanner):
    def __init__(self): 
        BaseScanner.__init__(self)
        self._name = "Generic SPECT Scanner"     
        self._scanner_type = "SPECT" 
        #name all parameters with 'self._p_xxx' 
        self._p_gantry_angular_positions      = 180   #[adim,integer]
        self._p_gantry_angular_position_first = 0.0   #[degrees]
        self._p_gantry_angular_position_last  = 358.0 #[degrees]
        self._p_acquisition_pattern = "equal_steps"                              
        self._p_scan_time_sec = 600.0                 #[seconds]
        self._p_radius_mm = 300.0                     #[mm]
        self._p_n_pix_x = 128                         #[adim]
        self._p_n_pix_y = 128                         #[adim]
        self._p_pix_size_x_mm = 2.5                   #[mm]
        self._p_pix_size_y_mm = 2.5                   #[mm]

        self._scintillator = Scintillators.Ideal() 
        self._collimator = Collimators.LEHR() 

    def get_gantry_angular_positions(self): 
        return self._p_gantry_angular_positions

    def set_gantry_angular_positions(self,n): 
        if not np.isscalar(n): 
            raise BadParameter('Expected a scalar value.')
        if not isinstance(n,type(1)): 
            raise BadParameter('Expected an integer value.')
        self._p_gantry_angular_positions = n 

    def get_gantry_angular_position_first(self): 
        return self._p_gantry_angular_position_first

    def set_gantry_angular_position_first(self,position_deg): 
        if not np.isscalar(position_deg): 
            raise BadParameter('Expected a scalar value.')
        self._p_gantry_angular_position_first = position_deg 

    def get_gantry_angular_position_last(self): 
        return self._p_gantry_angular_position_last

    def set_gantry_angular_position_last(self,position_deg): 
        if not np.isscalar(position_deg): 
            raise BadParameter('Expected a scalar value.')
        self._p_gantry_angular_position_last = position_deg 

    def get_acquisition_pattern(self): 
        return self._p_acquisition_pattern

    def set_acquisition_pattern(self,pattern_name):
        if str(pattern_name) == "equal_steps": 
            self._p_acquisition_pattern = str(pattern_name)
        #elif :  ..
        else: 
            raise BadParameter('Unknown acquisition pattern '+str(pattern_name))

    def get_scan_time(self): 
        return self._p_scan_time_sec

    def set_scan_time(self,scan_time_sec): 
        if not np.isscalar(scan_time_sec): 
            raise BadParameter('Expected a scalar value.')
        self._p_scan_time_sec = scan_time_sec

    def get_radius(self): 
        return self._p_radius_mm 

    def set_radius(self,radius_mm): 
        if not np.isscalar(radius_mm): 
            raise BadParameter('Expected a scalar value.') 
        self._p_radius_mm = radius_mm

    def get_n_pixels(self): 
        return (self._p_n_pix_x, self._p_n_pix_y)

    def set_n_pixels(self,n_pixels_x,n_pixels_y): 
        if (not np.isscalar(n_pixels_x)) or (not np.isscalar(n_pixels_y)): 
            raise BadParameter('Expected scalar values.') 
        self._p_n_pix_x = n_pixels_x
        self._p_n_pix_y = n_pixels_y

    def get_pixel_size(self): 
        return (self._p_pix_size_x, self._p_pix_size_y)

    def set_pixel_size(self,pixel_size_x,pixel_size_y): 
        if (not np.isscalar(pixel_size_x)) or (not np.isscalar(pixel_size_y)): 
            raise BadParameter('Expected scalar values.') 
        self._p_pix_size_x_mm = pixel_size_x
        self._p_pix_size_y_mm = pixel_size_y

    def get_scintillator(self):
        return self._scintillator 

    def set_scintillator(self,scintillator): 
        if not isinstance(scintillator,scintillators.BaseScintillatorSPECT): 
            raise BadParameter('Expected an instance of BaseScintillatorSPECT')
        self._scintillator = scintillator 

    def get_collimator(self): 
        return self._collimator 

    def set_collimator(self,collimator): 
        if not isinstance(collimator,collimators.BaseCollimatorSPECT): 
            raise BadParameter('Expected an instance of BaseCollimatorSPECT')
        self._collimator = collimator 
        
    def project(self): 
        pass

    def backproject(self):
        pass

    def scan(self,activity_Bq,scan_time_sec=None): 
        if scan_time_sec==None: 
            scan_time_sec = self.get_scan_time() 
        sinogram = 0
        return sinogram
        



class GE_Infinia(SPECT):
    def __init__(self): 
        SPECT.__init__(self)
        self._name = "GE Infinia SPECT Scanner with LEHR collimator"  






