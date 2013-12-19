
# tomi - Tomographic Inference 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston


__all__ = ['PET_Static_Scan','PET_Dynamic_Scan','PET_Interface_Petlink32','PET_Interface_mMR']


# Import interfile data handling module 
from interfile import Interfile

#Import scanner-specific interfaces: 
from petlink import PET_Interface_Petlink32

try: 
    from mMR import PET_Interface_mMR 
    HAVE_mMR = True
except:
    HAVE_mMR = False

# Import NiftyCore ray-tracers
try:
    from NiftyCore.NiftyRec import PET_project_compressed, PET_backproject_compressed 
    HAVE_NiftyRec = True
except: 
    print "NiftyCore could not be loaded: it will not be possible to reconstruct the PET data. "
    HAVE_NiftyRec = False

# Import ilang (inference language; optimisation) 
from PET_ilang import PET_Static_Poisson, PET_Dynamic_Poisson 

# Import DisplayNode to produce ipython notebook visualisations
from DisplayNode import DisplayNode

# Set verbose level
from tomi.verbose import *
set_verbose_no_printing()
#set_verbose_high()

import Image as PIL 
from numpy import isscalar, linspace, int32
import os



DEFAULT_BINNING = {      "n_axial":                252,
                         "n_azimuthal":            1,
                         "angular_step_axial":     1.0, 
                         "angular_step_azimuthal": 1.0,
                         "size_u":                 344.0,
                         "size_v":                 64.0,
                         "n_u":                    344,
                         "n_v":                    64,   }


DEFAULT_ROI = {          "x":                      0.0, 
                         "y":                      0.0, 
                         "z":                      0.0, 
                         "theta_x":                0.0,
                         "theta_y":                0.0, 
                         "theta_z":                0.0,   }

DEFAULT_N_TIME_BINS = 30

EPS = 0.000000000001



class FileNotFound(Exception): 
    def __init__(self,msg,filename): 
        self.msg = str(msg) 
        self.filename = str(filename)
    def __str__(self): 
        return "Cannot find file '%s' (%s)."%(self.filename, self.msg)




class ROI(): 
    """Region of Interest. """
    def __init__(self,dictionary=None):
        if dictionary: 
            self.load_from_dictionary(dictionary)
        else: 
            self.load_from_dictionary(DEFAULT_ROI)    
            
    def load_from_dictionary(self,dictionary):
        self.x       = dictionary['x']                           # Translation along x
        self.y       = dictionary['y']                           # Translation along y
        self.z       = dictionary['z']                           # Translation along z
        self.theta_x = dictionary['theta_x']                     # Rotation around x
        self.theta_y = dictionary['theta_y']                     # Rotation around y
        self.theta_z = dictionary['theta_z']                     # Rotation around z

    def __repr__(self):
        s = "PET volume location (ROI): \n"
        s = s+" - x:             %f \n"%self.x
        s = s+" - y:             %f \n"%self.y
        s = s+" - z:             %f \n"%self.z
        s = s+" - theta_x:       %f \n"%self.theta_x
        s = s+" - theta_y:       %f \n"%self.theta_y
        s = s+" - theta_z:       %f \n"%self.theta_z
        return s




class Binning(): 
    """PET detectors binning. """
    def __init__(self,dictionary=None): 
        self._initialised = 0
        self.name = "Unknown binning name"
        if dictionary: 
            self.load_from_dictionary(dictionary)        
        else: 
            self.load_from_dictionary(DEFAULT_BINNING)      
            
    def load_from_dictionary(self,dictionary):
        self.N_axial                = dictionary['n_axial']                    # Number of axial bins 
        self.N_azimuthal            = dictionary['n_azimuthal']                # Number of azimuthal bins
        self.angular_step_axial     = dictionary['angular_step_axial']         # Axial angular step, degrees 
        self.angular_step_azimuthal = dictionary['angular_step_azimuthal']     # Azimuthal angular step, degrees 
        self.size_u                 = dictionary['size_u']                     # Size of the detector plane, axis u, [adimensional]
        self.size_v                 = dictionary['size_v']                     # Size of the detector plane, axis v,  [adimensional]
        self.N_u                    = dictionary['n_u']                        # Number of pixels of the detector plan, axis u 
        self.N_v                    = dictionary['n_v']                        # Number of pixels of the detector plan, axis v 
        
    def __repr__(self): 
        s = "PET Binning: \n"        
        s = s+" - N_axial bins:             %d \n"%self.N_axial 
        s = s+" - N_azimuthal bins:         %d \n"%self.N_azimuthal 
        s = s+" - Axial angular step:       %f \n"%self.angular_step_axial 
        s = s+" - Azimuthal angular step:   %f \n"%self.angular_step_azimuthal 
        s = s+" - Size_u:                   %f \n"%self.size_u
        s = s+" - Size_v:                   %f \n"%self.size_v        
        s = s+" - N_u:                      %d \n"%self.N_u
        s = s+" - N_v:                      %d \n"%self.N_v
        return s


        
        
        
class PET_Static_Scan(): 
    """PET Static Scan. """
    def __init__(self): 
        self.interface     = None                               # PET scanner interface 
        self.binning       = Binning(DEFAULT_BINNING)           # PET detectors binning 

        self.time_start    = 0                                  # Time at start of scan 
        self.time_end      = 0                                  # Time at end of scan
        self.N_counts      = 0                                  # Number of photon counts 

        self._offsets          = None                           # 'offsets' and 'locations' define the structure of the sparse measurement (and projection) data
        self._locations        = None                           # 'offsets' and 'locations' define the structure of the sparse measurement (and projection) data
        self._measurement_data = None                           # measurement data, photon counts, locations are defined by 'offsets' and 'locations'

        self.scanner_detected = False                           # True if scanner model has been detected, False if unknown scanner model 
               
        self.enable_gpu_acceleration()                          # change to self.disable_gpu_acceleration() to disable by default 
        
        self._construct_ilang_model() 
        self._display_node = DisplayNode() 
        
    def _construct_ilang_model(self): 
        self.ilang_model = PET_Static_Poisson(self) 
    
    def set_binning(self, binning): 
        if isinstance(binning,Binning): 
            self.binning = binning
        else:
            self.binning = Binning(binning)
        
    def set_interface(self,interface): 
        self.interface = interface 

    def load_listmode_file(self, filename): 
        #FIXME: implement 
        self._construct_ilang_model()        
       
    def project(self,activity=None,roi=None,attenuation=None): 
        # FIXME: resample here according to ROI 
        projection_data = PET_project_compressed(activity,attenuation,self._offsets,self._locations,self.binning,self.gpu_acceleration) 
        return (projection_data, self._locations, self._offsets)

    def backproject(self, projection_data,roi=None,attenuation=None): 
        # FIXME: resample here according to ROI 
        R = PET_backproject_compressed(projection_data, attenuation, self._offsets, self._locations, self.binning, self.gpu_acceleration)
        return R['backprojection']
      
    def get_measurement(self): 
        return (self._measurement_data,self._locations,self._offsets)

    def uncompressed_measurement(self): 
        uncompressed_measurement = self.uncompress(self._measurement_data) 
        return uncompressed_measurement 
               
    def set_measurement_data(self,measurement_data): 
        self._measurement_data = measurement_data 

    def uncompress(self, projection_data): 
        return self.interface.uncompress(self._offsets, projection_data, self._locations, self.binning.N_u, self.binning.N_v)

    def enable_gpu_acceleration(self): 
        self.gpu_acceleration = 1
        
    def disable_gpu_acceleration(self): 
        self.gpu_acceleration = 0 
        
    def load_static_measurement(self, time_bin=None): 
        if time_bin == None: 
            R = self.interface.get_measurement_static() 
        else: 
            R = self.interface.get_measurement(time_bin)
        self.time_start        = R['time_start']
        self.time_end          = R['time_end']
        self.N_counts          = R['N_counts']        
        self._offsets          = R['offsets']
        self._locations        = R['locations'] 
        self._measurement_data = R['counts'] 
        self._construct_ilang_model()
        
    def __repr__(self): 
        s = "Static PET acquisition:  \n" 
        s = s+" - Time_start:                   %d [ms] \n"%self.time_start 
        s = s+" - Time_end:                     %d [ms] \n"%self.time_end 
        s = s+" - N_counts:                     %d      \n"%self.N_counts 
        s = s+" = Interface: \n"
        s = s+"     - Name:                     %s \n"%self.interface.name 
        s = s+"     - Manufacturer:             %s \n"%self.interface.manufacturer 
        s = s+"     - Version:                  %s \n"%self.interface.version 
        s = s+" * Binning: \n"        
        s = s+"     - N_axial bins:             %d \n"%self.binning.N_axial 
        s = s+"     - N_azimuthal bins:         %d \n"%self.binning.N_azimuthal 
        s = s+"     - Axial angular step:       %f \n"%self.binning.angular_step_axial 
        s = s+"     - Azimuthal angular step:   %f \n"%self.binning.angular_step_azimuthal 
        s = s+"     - Size_u:                   %f \n"%self.binning.size_u
        s = s+"     - Size_v:                   %f \n"%self.binning.size_v        
        s = s+"     - N_u:                      %s \n"%self.binning.N_u
        s = s+"     - N_v:                      %s \n"%self.binning.N_v
        return s

    def __del__(self):
        del(self.interface) 


        
        
        
class PET_Dynamic_Scan(): 
    """PET Dynamic Scan. """
    def __init__(self): 
        self.interface     = None                        # PET Scanner interface 
        self.binning       = Binning(DEFAULT_BINNING)    # PET detectors binning 
        self.static        = None                        # Static scan (no time binning)
        self.dynamic       = []                          # Sequence of static scans, one per time bin 
        self.time_bins     = []                          # Time binning (array, [ms])

        self.N_time_bins   = 0                           # Number of time bins         
        self.N_counts      = 0                           # Total number of counts 
        self.time_start    = 0                           # Time at start of scan [ms]
        self.time_end      = 0                           # Time at end of scan [ms] 
        self.scanner_detected = False                    # Scanner model detected (False if unknown scanner model) 

        self._static_measurement_data = None
        self._offsets                 = None
        self._locations               = None 

        self._construct_ilang_model()

    def set_binning(self, binning): 
        if isinstance(binning,Binning): 
            self.binning = binning
        else:    
            self.binning = Binning(binning)

    def set_interface(self,interface): 
        self.interface = interface 

    def _construct_ilang_model(self):
        self.ilang_model = PET_Dynamic_Poisson(self)             
        
    def load_listmode_file(self, hdr_filename, time_bins=None, data_filename=None): 
        """Load measurement data from a listmode file. """
        print_debug("- Loading dynamic PET data from listmode file "+str(hdr_filename) )
        hdr = Interfile.load(hdr_filename) 
        #Extract information from the listmode header
        # 1) Determine the model of the scanner: 
        self.scanner_detected = False
        if hdr.has_key('originating system'): 
            if hdr['originating system']['value'] == 2008: 
                print_debug("- Detected Siemens Biograph mMR scanner. ")
                self.set_interface(PET_Interface_mMR())
                self.scanner_detected = True 
        if not self.scanner_detected: 
            print_debug("- No scanner detected, assuming petlink32 listmode. ")
            self.set_interface(PET_Interface_Petlink32()) 
        
        # 2) Guess the path of the listmode data file, if not specified or mis-specified; 
        if data_filename==None: 
            data_filename      = hdr['name of data file']['value'] 
        if not os.path.exists(data_filename): 
            # cross platform compatibility: 
            data_filename = data_filename.replace("/",os.path.sep).replace("\\",os.path.sep)            
            # search in the same location as the header file: 
            data_filename = os.path.split(hdr_filename)[0]+os.path.sep+os.path.split(data_filename)[-1] 
        if not os.path.exists(data_filename): 
            raise FileNotFound("listmode data",data_filename) 

        # 3) Determine acquisition settings
        n_packets              = hdr['total listmode word counts']['value'] 
        scan_duration          = hdr['image duration']['value']*1000            # now milliseconds
        
        # 4) determine scanner parameters
        n_radial_bins          = hdr['number of projections']['value'] 
        n_angles               = hdr['number of views']['value'] 
        n_rings                = hdr['number of rings']['value'] 
        max_ring_diff          = hdr['maximum ring difference']['value']
        n_sinograms            = n_rings+2*n_rings*max_ring_diff-max_ring_diff**2-max_ring_diff

        # Determine the time binning pattern 
        if time_bins == None: 
            time_bins = int32(linspace(0,scan_duration,DEFAULT_N_TIME_BINS+1))
        elif isscalar(time_bins):    #time_bins in this case indicates the number of time bins
            time_bins = int32(linspace(0,scan_duration,time_bins+1)) 

        # Display information 
        print_debug(" - Number of packets:   %d      "%n_packets )
        print_debug(" - Scan duration:       %d [sec]"%scan_duration )
        print_debug(" - Listmode data file:  %s      "%data_filename )
        print_debug(" - Number of time bins: %d      "%(len(time_bins)-1) )
        print_debug(" - Time start:          %f [sec]"%(time_bins[ 0]/1000.0) )
        print_debug(" - Time end:            %f [sec]"%(time_bins[-1]/1000.0) )

        # Load the listmode data 
        R = self.interface.load_listmode(data_filename,n_packets,time_bins,self.binning,n_radial_bins, n_angles, n_sinograms) 

        self.N_counts     = R['N_counts'] 
        self.N_time_bins  = R['N_time_bins'] 
        self.time_start   = R['time_start'] 
        self.time_end     = R['time_end'] 
        self.time_bins    = time_bins[0:self.N_time_bins+1]  #the actual time bins are less than the requested time bins, truncate time_bins 
        
        # Make list of PET_Static_Scan objects, one per bin
        self.dynamic = [] 
        for t in range(self.N_time_bins): 
            PET_t = PET_Static_Scan() 
            PET_t.set_interface(self.interface) 
            PET_t.set_binning(self.binning) 
            PET_t.load_static_measurement(t) 
            PET_t.scanner_detected = self.scanner_detected
            self.dynamic.append(PET_t) 
        
        # Make a global PET_Static_Scan object 
        self.static = PET_Static_Scan()
        self.static.set_interface(self.interface) 
        self.static.set_binning(self.binning) 
        self.static.load_static_measurement() 
        self.static.scanner_detected = self.scanner_detected 
        
        # Construct ilang model 
        self._construct_ilang_model()

        # Load static measurement data 
        self.load_static_measurement() 

    def load_static_measurement(self): 
        R = self.interface.get_measurement_static() 
        self.time_start               = R['time_start']
        self.time_end                 = R['time_end']
        self.N_counts                 = R['N_counts']        
        self._offsets                 = R['offsets']
        self._locations               = R['locations'] 
        self._static_measurement_data = R['counts'] 

    def get_static_measurement(self): 
        return (self._static_measurement_data,self._locations,self._offsets)

    def uncompressed_measurement(self): 
        uncompressed_measurement = self.uncompress(self._static_measurement_data) 
        return uncompressed_measurement 
               
    def uncompress(self, projection_data): 
        return self.interface.uncompress(self._offsets, projection_data, self._locations, self.binning.N_u, self.binning.N_v) 
                      
    def __repr__(self): 
        """Display information about Dynamic_PET_Scan"""
        s = "Dynamic PET acquisition:  \n" 
        s = s+" - N_time_bins:                  %d      \n"%self.N_time_bins 
        s = s+" - Time_start:                   %d [ms] \n"%self.time_start 
        s = s+" - Time_end:                     %d [ms] \n"%self.time_end 
        s = s+" - N_counts:                     %d      \n"%self.N_counts 
        s = s+" - Mean time bin duration:       %d [ms] \n"%0 #FIXME 
        s = s+" * Interface: \n" 
        s = s+"     - Name:                     %s \n"%self.interface.name 
        s = s+"     - Manufacturer:             %s \n"%self.interface.manufacturer 
        s = s+"     - Version:                  %s \n"%self.interface.version 
        s = s+" * Binning: \n" 
        s = s+"     - N_axial bins:             %d \n"%self.binning.N_axial 
        s = s+"     - N_azimuthal bins:         %d \n"%self.binning.N_azimuthal 
        s = s+"     - Axial angular step:       %f \n"%self.binning.angular_step_axial 
        s = s+"     - Azimuthal angular step:   %f \n"%self.binning.angular_step_azimuthal 
        s = s+"     - Size_u:                   %f \n"%self.binning.size_u
        s = s+"     - Size_v:                   %f \n"%self.binning.size_v        
        s = s+"     - N_u:                      %s \n"%self.binning.N_u
        s = s+"     - N_v:                      %s \n"%self.binning.N_v
        return s

    def __iter__(self): 
        """This method makes the object iterable. """
        return iter(self.dynamic)
    
    def __getitem__(self,i): 
        """This method makes the object addressable like a list. """
        return self.dynamic[i] 

    def __del__(self):
        """Delete interface when the object is deleted: interface needs explicit 
        deletion in order to manage C library memory deallocation  """
        self.interface.free_memory()

        
        
