
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston


__all__ = ['PET_Static_Scan','PET_Dynamic_Scan','PET_Interface_Petlink32','PET_Interface_mMR','Binning']



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
from PET_ilang import PET_Static_Poisson, PET_Dynamic_Poisson, ProbabilisticGraphicalModel
from ilang.Samplers import Sampler 

# Import DisplayNode to produce ipython notebook visualisations
from DisplayNode import DisplayNode

# Import occiput: 
from occiput.Core import Volume, RigidVolume
from occiput.Visualization import ProgressBar
from occiput.DataSources.Synthetic import uniform_cylinder

# Import other modules
import Image as PIL 
import ImageDraw
from numpy import isscalar, linspace, int32, uint32, ones, zeros, pi, float32, where, ndarray
from numpy.random import randint 
import os
import svgwrite
import ipy_table 

# Set verbose level
from occiput.verbose import *
set_verbose_no_printing()
#set_verbose_high()




DEFAULT_BINNING = {      "n_axial":                252,
                         "n_azimuthal":            1,
                         "angular_step_axial":     pi/252, 
                         "angular_step_azimuthal": 0.0,
                         "size_u":                 437.98,
                         "size_v":                 256.0,
                         "n_u":                    344,
                         "n_v":                    64,   }


DEFAULT_ROI = {          "x":                      0.0, 
                         "y":                      0.0, 
                         "z":                      0.0, 
                         "theta_x":                0.0,
                         "theta_y":                0.0, 
                         "theta_z":                0.0,   }

                         
DEFAULT_PROJECTION_PARAMETERS  =  {
                         "N_samples":                128, 
                         "sample_step":              5, 
                         "background_activity":      0.0, 
                         "background_attenuation":   0.0, 
                         "truncate_negative_values": 0,  
                         "gpu_acceleration":         1,      
                         "direction":                2,      # affects performance only, between 1 and 6; 2 and 5 are normally the best value
                         "block_size":               512 }   # affects performance only

DEFAULT_BACKPROJECTION_PARAMETERS  =  {
                         "N_samples":                128, 
                         "sample_step":              5, 
                         "background_activity":      0.0, 
                         "background_attenuation":   0.0, 
                         "truncate_negative_values": 0,  
                         "gpu_acceleration":         1,      
                         "direction":                1,      # affects performance only, between 1 and 6; 2 and 5 are normally the best value
                         "block_size":               512 }   # affects performance only

DEFAULT_N_TIME_BINS = 30

EPS = 1e-6




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
        
        



def millisec_to_min_sec(ms):
    totsec = int(ms)/1000
    sec = totsec%60
    min = int(ms)/1000/60 
    msec = ms%(1000)
    return str(min).zfill(2)+" [min]  "+ str(sec).zfill(2)+" [sec]  " + str(msec).zfill(3)+" [msec]"

def pretty_print_large_number(number):
    s = str(number).ljust(12)
    if number > 0 and number < 1e3: 
        pass 
    elif number >= 1e3 and number < 1e6: 
        s = s + "  (%3.1f Thousand)"%(number*1.0/1e3)
    elif number >= 1e6 and number < 1e9: 
        s = s + "  (%3.1f Million)"%(number*1.0/1e6)
    elif number >= 1e9 and number < 1e12: 
        s = s + "  (%3.1f Billion)"%(number*1.0/1e9)
    elif number >= 1e12 and number < 1e15: 
        s = s + "  (%3.1f Trillion)"%(number*1.0/1e12)
    return s

def print_percentage(number):
    return "%2.2f %%"%((1.0*number)*100)
    






class ROI(): 
    """Region of Interest. """
    def __init__(self,parameters=None): 
        if parameters==None: 
            self.load_from_dictionary(DEFAULT_ROI)  
        elif type(parameters) == dict:  
            self.load_from_dictionary(parameters)        
        elif type(parameters) in [list, tuple]: 
            if len(parameters) == len(DEFAULT_ROI.keys()): 
                self.x       = parameters[0]
                self.y       = parameters[1]
                self.z       = parameters[2]
                self.theta_x = parameters[3]
                self.theta_y = parameters[4]
                self.theta_z = parameters[5]
            else: 
                raise UnknownParameter('Parameter %s specified for the construction of Binning is not compatible. '%str(parameters)) 
        else: 
            raise UnknownParameter('Parameter %s specified for the construction of Binning is not compatible. '%str(parameters)) 
    
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

    def _repr_html_(self):
        table_data = [['x',self.x],['y',self.y],['z',self.z],['theta_x',self.theta_x],['theta_y',self.theta_y],['theta_z',self.theta_z]] 
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic_left')
        #table = ipy_table.set_column_style(0, color='lightBlue')
        table = ipy_table.set_global_style(float_format="%3.3f")        
        return table._repr_html_()






class Binning(): 
    """PET detectors binning. """
    def __init__(self,parameters=None): 
        self._initialised = 0
        self.name = "Unknown binning name"
        if parameters==None: 
            self.load_from_dictionary(DEFAULT_BINNING)  
        elif type(parameters) == dict:  
            self.load_from_dictionary(parameters)        
        elif type(parameters) in [list, tuple]: 
            if len(parameters) == len(DEFAULT_BINNING.keys()): 
                self.N_axial                = parameters[0]
                self.N_azimuthal            = parameters[1]
                self.angular_step_axial     = parameters[2] 
                self.angular_step_azimuthal = parameters[3]
                self.size_u                 = parameters[4]
                self.size_v                 = parameters[5] 
                self.N_u                    = parameters[6]
                self.N_v                    = parameters[7] 
            else: 
                raise UnknownParameter('Parameter %s specified for the construction of Binning is not compatible. '%str(parameters)) 
        else: 
            raise UnknownParameter('Parameter %s specified for the construction of Binning is not compatible. '%str(parameters)) 
            
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
        s = s+" - N_axial_bins:             %d \n"%self.N_axial 
        s = s+" - N_azimuthal_bins:         %d \n"%self.N_azimuthal 
        s = s+" - Angular_step_axial:       %f \n"%self.angular_step_axial 
        s = s+" - Angular_step_azimuthal:   %f \n"%self.angular_step_azimuthal 
        s = s+" - Size_u:                   %f \n"%self.size_u
        s = s+" - Size_v:                   %f \n"%self.size_v        
        s = s+" - N_u:                      %d \n"%self.N_u
        s = s+" - N_v:                      %d \n"%self.N_v
        return s

    def _repr_html_(self):
        table_data = [['N_axial',self.N_axial],['N_azimuthal',self.N_azimuthal],['Angular_step_axial',self.angular_step_axial],
        ['Angular_step_azimuthal',self.angular_step_azimuthal],['Size_u',self.size_u],['Size_v',self.size_v],['N_u',self.N_u],['N_v',self.N_v]] 
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic_left')
        #table = ipy_table.set_column_style(0, color='lightBlue')
        table = ipy_table.set_global_style(float_format="%3.3f")        
        return table._repr_html_()
        
    def _make_svg(self):
        w = '100%'
        h = '100%'
        
        dwg = svgwrite.Drawing('test.svg',size=(w,h), profile='full', debug=True)
        dwg.viewbox(width=100, height=100)
        
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





class ProjectionParameters(): 
    default_parameters = DEFAULT_PROJECTION_PARAMETERS
    def __init__(self,parameters=None): 
        self._initialised = 0
        self.name = "Unknown binning name"
        if parameters==None: 
            self.load_from_dictionary(self.default_parameters)  
        elif type(parameters) == dict:  
            self.load_from_dictionary(parameters)        
        elif type(parameters) in [list, tuple]: 
            if len(parameters) == len(self.default_parameters.keys()): 
                self.N_samples                = parameters[0] 
                self.sample_step              = parameters[1]
                self.background_activity      = parameters[2]
                self.background_attenuation   = parameters[3] 
                self.truncate_negative_values = parameters[4]
                self.gpu_acceleration         = parameters[5]
                self.direction                = parameters[6]
                self.block_size               = parameters[7]
            else: 
                raise UnknownParameter('Parameter %s specified for ProjectionParameters is not compatible. '%str(parameters)) 
        else: 
            raise UnknownParameter('Parameter %s specified for ProjectionParameters is not compatible. '%str(parameters)) 
            
    def load_from_dictionary(self,dictionary):
        self.N_samples                = dictionary['N_samples']                 # Number of samples along a line when computing line integrals
        self.sample_step              = dictionary['sample_step']               # distance between consecutive points along a line when computing line integrals (this is in the same unit measure as the size of the imaging volume (activity_size and attenuation_size))
        self.background_activity      = dictionary['background_activity']       # Activity in voxels outside of the imaging volume
        self.background_attenuation   = dictionary['background_attenuation']    # Attenuation in voxels outside of the imaging volume 
        self.truncate_negative_values = dictionary['truncate_negative_values']  # If set to 1, eventual negative values obtained when projecting are set to 0
                                                                                # (This is meant to remove eventual unwanted small negative values due to FFT-based smoothing within the projection algorithm 
                                                                                # - note that numerical errors may produce small negative numbers when doing FFT-IFFT even if the function is all positive )
        self.gpu_acceleration         = dictionary['gpu_acceleration']          # Whether to use GPU acceleration (if available) or not 
        self.direction                = dictionary['direction']                 # Direction parameter for the projector and b-proj., it affects performance only; 2 and 5 are normally the best values. 
        self.block_size               = dictionary['block_size']                # Number of blocks per thread for the GPU accelerated projector and b-proj. It affects performance only. Typical values: 256, 512, 768
        




class BackprojectionParameters(ProjectionParameters): 
    default_parameters = DEFAULT_BACKPROJECTION_PARAMETERS

      
      
    
    
    
class SubsetGenerator():  
    def __init__(self,N_axial,N_azimuthal):
        self._N_axial     = N_axial
        self._N_azimu     = N_azimuthal

    def new_subset(self,mode,subset_size): 
        if mode=='random': 
            return self._random_no_replacement(subset_size) 
        elif mode=='ordered': 
            raise UnexpectedParameter("'%s' subset selection mode not yet supported."%str(mode))
        else: 
            raise UnexpectedParameter("'mode' parameter %s not recognised."%str(mode))
        
    def all_active(self): 
        return ones((self._N_axial,self._N_azimu),dtype=uint32) 

    def _random_no_replacement(self,subset_size): 
        if subset_size>=self._N_axial*self._N_azimu: 
            return self.all_active() 
        M = zeros((self._N_axial,self._N_azimu),dtype=int32) 
        n = 0
        while n<subset_size:
            active_axial = randint(self._N_axial)
            active_azimu = randint(self._N_azimu) 
            if M[active_axial,active_azimu] == 0: 
                M[active_axial,active_azimu] = 1
                n+=1
        return M



        
class PET_Static_Scan(): 
    """PET Static Scan. """
    def __init__(self): 
        self.interface     = None                                        # PET scanner interface 
        self.binning       = self.set_binning(Binning(DEFAULT_BINNING))  # PET detectors binning 
        
        self.time_start    = 0                                  # Time at start of scan 
        self.time_end      = 0                                  # Time at end of scan
        self.N_counts      = 0                                  # Number of photon counts 
        self.N_locations   = 0                                  # Number of active locations of interaction
        self.compression_ratio = 0.0                            # Compression ratio of projection data 
        self.listmode_loss     = 0.0                            # Loss of data when events are binned from listmode 

        self._offsets          = None                           # 'offsets' and 'locations' define the structure of the sparse measurement (and projection) data
        self._locations        = None                           # 'offsets' and 'locations' define the structure of the sparse measurement (and projection) data
        self._measurement_data = None                           # measurement data, photon counts, locations are defined by 'offsets' and 'locations'

        self.activity_shape    = [128,128,128]  #FIXME: have a default value (from dictionary)
        self.activity_size     = [256,256,256]  #FIXME: have a default value (from dictionary), but adapt to the detector size, and also have a set method 
        self.attenuation_shape = [128,128,128]  #FIXME: have a default value (from dictionary)
        self.attenuation_size  = [256,256,256]  #FIXME: see previous line  
        self.scanner_detected  = False                           # True if scanner model has been detected, False if unknown scanner model 

        self.projection_parameters     = ProjectionParameters()       
        self.backprojection_parameters = BackprojectionParameters()   
        self.enable_gpu_acceleration()                          # change to self.disable_gpu_acceleration() to disable by default      

        self._construct_ilang_model() 
        self._display_node = DisplayNode() 
        
        self._normalization             = None                  # Normalization volume 
        self._need_normalization_update = True                  # If True, the normalization volume needs to be recomputed 



    def _construct_ilang_model(self): 
        # define the ilang probabilistic model 
        self.ilang_model = PET_Static_Poisson(self) 
        # construct a basic Directed Acyclical Graph
        self.dag         = ProbabilisticGraphicalModel(['lambda','alpha','z']) 
        self.dag.set_nodes_given(['z','alpha'],True) 
        self.dag.add_dependence(self.ilang_model,{'lambda':'lambda','alpha':'alpha','z':'z'}) 
        # construct a basic sampler object
        self.sampler     = Sampler(self.dag)
 
    def set_binning(self, binning): 
        if isinstance(binning,Binning): 
            self.binning = binning
        else:
            self.binning = Binning(binning)
        self._subset_generator = SubsetGenerator(self.binning.N_axial,self.binning.N_azimuthal) 
    
    def set_interface(self,interface): 
        self.interface = interface 

    def load_listmode_file(self, hdr_filename, data_filename=None): 
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
        #  1 - see if the specified listmode data file exists 
        if data_filename!=None: 
            data_filename = data_filename.replace("/",os.path.sep).replace("\\",os.path.sep)          # cross platform compatibility 
            if not os.path.exists(data_filename): 
                raise FileNotFound("listmode data",data_filename)  
        #  2 - if the listmode data file is not specified, try with the name (and full path) contained in the listmode header
        data_filename      = hdr['name of data file']['value']
        data_filename = data_filename.replace("/",os.path.sep).replace("\\",os.path.sep)              # cross platform compatibility
        if not os.path.exists(data_filename): 
        #  3 - if it doesn't exist, look in the same path as the header file for the listmode data file with name specified in the listmode header file 
            data_filename = os.path.split(hdr_filename)[0]+os.path.sep+os.path.split(data_filename)[-1]  
            if not os.path.exists(data_filename): 
        #  4 - if it doesn't exist, look in the same path as the header file for the listmode data file with same name as the listmode header file, replacing the extension: ".l.hdr -> .l" 
                if hdr_filename.endswith(".l.hdr"): 
                    data_filename = hdr_filename.replace(".l.hdr",".l") 
                    if not os.path.exists(data_filename):     
                        raise FileNotFound("listmode data",data_filename)  
        #  5 - if it doesn't exist, look in the same path as the header file for the listmode data file with same name as the listmode header file, replacing the extension: ".hdr -> .l" 
                elif hdr_filename.endswith(".hdr"): 
                    data_filename = hdr_filename.replace(".hdr",".l") 
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

        time_bins = int32(linspace(0,scan_duration,2))

        # Display information 
        print_debug(" - Number of packets:   %d      "%n_packets )
        print_debug(" - Scan duration:       %d [sec]"%scan_duration )
        print_debug(" - Listmode data file:  %s      "%data_filename )
        print_debug(" - Number of time bins: %d      "%(len(time_bins)-1) )
        print_debug(" - Time start:          %f [sec]"%(time_bins[ 0]/1000.0) )
        print_debug(" - Time end:            %f [sec]"%(time_bins[-1]/1000.0) )

        # Load the listmode data 
        R = self.interface.load_listmode(data_filename,n_packets,time_bins,self.binning,n_radial_bins, n_angles, n_sinograms) 

        self.N_counts          = R['N_counts'] 
        self.N_locations       = R['N_locations']
        self.compression_ratio = R['compression_ratio']
        self.listmode_loss     = R['listmode_loss']
        self.N_time_bins  = 1
        self.time_start   = R['time_start'] 
        self.time_end     = R['time_end'] 

        self.load_static_measurement(0)
        
        # Construct ilang model 
        self._construct_ilang_model()
        return self 
        
    def project(self,activity,attenuation=None,roi_activity=None,roi_attenuation=None,offsets=None,locations=None,subsets_matrix=None): 
        if isinstance(activity,ndarray): 
            activity = float32(activity)
        else: 
            activity = float32(activity.data)
        if attenuation!=None:
            if isinstance(attenuation,ndarray):
                attenuation = float32(attenuation)
            else: 
                attenuation = float32(attenuation.data)
        if not list(activity.shape) == list(self.activity_shape): 
            raise UnexpectedParameter("Activity must have the same shape as self.activity_shape")
        if attenuation!=None: 
            if not list(attenuation.shape) == list(self.attenuation_shape): 
                raise UnexpectedParameter("Attenuation must have the same shape as self.attenuation_shape")
        if offsets==None: 
            offsets=self._offsets
        if locations==None:
            locations=self._locations
        # By default, the center of the imaging volume is at the center of the scanner 
        if roi_activity == None: 
            roi_activity = ROI((0.5*self.activity_size[0],0.5*self.activity_size[1],0.5*self.activity_size[2],0,0,0))
        # By default, the center of the imaging volume is at the center of the scanner 
        if roi_attenuation == None: 
            roi_attenuation = ROI((0.5*self.attenuation_size[0],0.5*self.attenuation_size[1],0.5*self.attenuation_size[2],0,0,0))    
        # Handle no subsets
        if subsets_matrix==None: 
            subsets_matrix=self._subsets_generator.all_active()    
        # Pass on to the C library all the parameters required for projection 
        projection_data = PET_project_compressed(activity,attenuation,offsets,locations, subsets_matrix, 
            self.binning.N_axial, self.binning.N_azimuthal, 
            self.binning.angular_step_axial, self.binning.angular_step_azimuthal, 
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.activity_size[0], self.activity_size[1], self.activity_size[2], 
            self.attenuation_size[0], self.attenuation_size[1], self.attenuation_size[2], 
            roi_activity.x, roi_activity.y, roi_activity.z, roi_activity.theta_x, roi_activity.theta_y, roi_activity.theta_z, 
            roi_attenuation.x, roi_attenuation.y, roi_attenuation.z, roi_attenuation.theta_x, roi_attenuation.theta_y, roi_attenuation.theta_z, 
            self.projection_parameters.gpu_acceleration, self.projection_parameters.N_samples, self.projection_parameters.sample_step, 
            self.projection_parameters.background_activity, self.projection_parameters.background_attenuation, self.projection_parameters.truncate_negative_values,
            self.projection_parameters.direction, self.projection_parameters.block_size)
#        return (projection_data, self._locations, self._offsets) 
        return projection_data 

    def backproject(self, projection_data, attenuation=None, roi_activity=None, roi_attenuation=None, offsets=None, locations=None, subsets_matrix=None): 
        projection_data = float32(projection_data)
        if attenuation!=None:
            attenuation = float32(attenuation)
        if attenuation!=None: 
            if not list(attenuation.shape) == list(self.attenuation_shape): 
                raise UnexpectedParameter("Activity must have the same shape as self.attenuation_shape")
        # By default, the center of the imaging volume is at the center of the scanner 
        if roi_activity == None: 
            roi_activity = ROI((0.5*self.activity_size[0],0.5*self.activity_size[1],0.5*self.activity_size[2],0,0,0))
        # By default, the center of the imaging volume is at the center of the scanner 
        if roi_attenuation == None: 
            roi_attenuation = ROI((0.5*self.attenuation_size[0],0.5*self.attenuation_size[1],0.5*self.attenuation_size[2],0,0,0))  
        if offsets==None: 
            offsets=self._offsets
        if locations==None:
            locations=self._locations   
        # Handle no subsets
        if subsets_matrix==None: 
            subsets_matrix=self._subsets_generator.all_active()
        # Pass on to the C library all the parameters required for back-projection 
        backprojection = PET_backproject_compressed(projection_data,attenuation,offsets,locations, subsets_matrix, 
            self.binning.N_axial, self.binning.N_azimuthal, self.binning.angular_step_axial, self.binning.angular_step_azimuthal, 
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.activity_shape[0], self.activity_shape[1], self.activity_shape[2], 
            self.activity_size[0], self.activity_size[1], self.activity_size[2], 
            self.attenuation_size[0], self.attenuation_size[1], self.attenuation_size[2], 
            roi_activity.x, roi_activity.y, roi_activity.z, roi_activity.theta_x, roi_activity.theta_y, roi_activity.theta_z, 
            roi_attenuation.x, roi_attenuation.y, roi_attenuation.z, roi_attenuation.theta_x, roi_attenuation.theta_y, roi_attenuation.theta_z, 
            self.backprojection_parameters.gpu_acceleration, self.backprojection_parameters.N_samples, self.backprojection_parameters.sample_step, 
            self.backprojection_parameters.background_activity, self.backprojection_parameters.background_attenuation, 
            self.backprojection_parameters.direction, self.backprojection_parameters.block_size)
        return RigidVolume(backprojection)
      
    def get_measurement(self): 
        return (self._measurement_data,self._locations,self._offsets)

    def uncompressed_measurement(self): 
        uncompressed_measurement = self.uncompress(self._measurement_data) 
        return uncompressed_measurement 
               
    def set_measurement_data(self,measurement_data): 
        self._measurement_data = measurement_data 

    def uncompress(self, projection_data, offsets=None, locations=None, N_u=None, N_v=None):
        if offsets==None: 
            offsets=self._offsets
        if locations==None:
            locations=self._locations 
        if N_u==None:
            N_u=self.binning.N_u
        if N_v==None:
            N_v=self.binning.N_v
        return self.interface.uncompress(offsets, projection_data, locations, N_u, N_v)

    def enable_gpu_acceleration(self): 
        self.projection_parameters.gpu_acceleration = 1 
        self.backprojection_parameters.gpu_acceleration = 1 

    def disable_gpu_acceleration(self): 
        self.projection_parameters.gpu_acceleration = 0 
        self.backprojection_parameters.gpu_acceleration = 0 

    def load_static_measurement(self, time_bin=None): 
        if time_bin == None: 
            R = self.interface.get_measurement_static() 
        else: 
            R = self.interface.get_measurement(time_bin) 
        self.time_start        = R['time_start'] 
        self.time_end          = R['time_end'] 
        self.N_counts          = R['N_counts']         
        self.N_locations       = R['N_locations']
        self.compression_ratio = R['compression_ratio']
        self.listmode_loss     = R['listmode_loss']
        self._offsets          = R['offsets'] 
        self._locations        = R['locations'] 
        self._measurement_data = R['counts'] 
        self._construct_ilang_model() 

    def set_full_sampling(self): 
        R = self.interface.full_sampling(self.binning.N_axial,self.binning.N_azimuthal,self.binning.N_u,self.binning.N_v) 
        self._offsets          = R['offsets'] 
        self._locations        = R['locations'] 
        self.N_locations       = self._locations.shape[0]

#    def reconstruct_activity(self,Nx,Ny,Nz,method='MLEM',iterations=100,roi=None,attenuation=None,roi_attenuation=None): 
#        # Use ilang -inference language - to do the reconstruction
#        # (if emission data has been loaded, then the instance of the object of this class has a 'dag' attribute)
#        if method == 'MLEM': 
#           # initialise the ilang variables 
#            if self.dag.get_node_value('lambda') == None: 
#                self.dag.set_node_value('lambda',ones((Nx,Ny,Nz),dtype=float32)) 
#            if self.dag.get_node_value('alpha') == None: 
#                self.dag.set_node_value('alpha',zeros((Nx,Ny,Nz),dtype=float32)) 
#            # reconstruct
#            self.sampler.set_node_sampling_method_manual('lambda','ExpectationMaximization')
#            self.sampler.sample_node('lambda', nsamples=iterations, trace=False)
#        else: 
#            raise UnknownParameter('method %s unknown'%str(method)) 
#        return self.dag.get_node_value('lambda')
        
    def set_activity(self,activity): 
        self.dag.set_node_value('lambda',activity) 
        # FIXME: how about the roi ? 

    def set_attenuation(self,attenuation): 
        self.dag.set_node_value('alpha',attenuation) 
        # FIXME: how about the roi ? 

    def get_normalization(self): 
        if self._normalization == None: 
           self._update_normalization() 
        elif self._need_normalization_update: 
           self._update_normalization()
        return self._normalization
    
    def _update_normalization(self): 
        self._normalization = self.backproject(ones((1,self.N_locations),dtype=float32)) 
        self._normalization.data = self._normalization.data + EPS
        self._need_normalization_update = False 

    def get_gradient(self,activity): 
        proj = self.project(activity)
        norm = self.get_normalization() 
        return RigidVolume(-norm.data + self.backproject( (self._measurement_data+EPS)/(proj+EPS) ).data)

    def get_mask(self): 
        if not hasattr(self,"_mask"): 
            self._mask = uniform_cylinder(self.activity_shape,self.activity_size,[0.5*self.activity_size[0],0.5*self.activity_size[1],0.5*self.activity_size[2]], 0.5*min(self.activity_size[0],self.activity_size[1]),self.activity_size[2],2,1,0)
        return self._mask
        
    def estimate_activity(self,iterations = 6, subset_size=None, subset_mode='random'): 
        progress_bar = ProgressBar() 
        progress_bar.set_percentage(0.1)
        # Subsets: 
        if subset_size==None:
            subsets_matrix=None
        else: 
            subsets_matrix = self._subset_generator.new_subset(subset_mode,subset_size)
        activity = ones(self.activity_shape,dtype=float32)
        for i in range(iterations):
            proj = self.project(activity,subsets_matrix=subsets_matrix)
            if subset_size==None:
                norm = self.get_normalization()  
            else: 
                norm = self.backproject(ones((1,self.N_locations),dtype=float32), subsets_matrix=subsets_matrix) 
            update = (self.backproject( (self._measurement_data+EPS)/(proj+EPS) ,subsets_matrix=subsets_matrix).data+EPS) / (norm.data +EPS) 
            
            #FIXME: find out why it crashes without this: 
            update[where(update>100000)] = 0.1
            activity = activity * update * self.get_mask().data

            progress_bar.set_percentage((i+1)*100.0/iterations) 
            #print "Iteration: %d    max act: %f    min act: %f    max proj: %f    min proj: %f    max norm: %f    min norm: %f"%(i, activity.max(), activity.min(), proj.max(), proj.min(), norm.data.max(), norm.data.min() )
        progress_bar.set_percentage(100.0)
        return RigidVolume(activity)
            
    def volume_render(self,volume,scale=None): 
        # FIXME: make a VolumeRenderer class and use it, the following is a quick fix: 
        R = self.interface.full_sampling(180,1,256,256) 
        offsets   = R['offsets']
        locations = R['locations']
        #proj = self.project(volume,offsets=offsets,locations=locations) 
        if isinstance(volume,ndarray): 
            volume = float32(volume)
        else: 
            volume = float32(volume.data)
        proj = PET_project_compressed(volume,None,offsets,locations, 
            180, 1, pi/180, 0.0, 256,256, 256.0, 256.0, 
            256.0, 256.0, 256.0,  
            256.0, 256.0, 256.0, 
            128.0, 128.0, 128.0, 0.0, 0.0, 0.0,  
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            1, 256, 1.5, 
            0.0, 0.0, 0,
            2, 512)
            
        proj[where(proj>proj.max()/scale )]=proj.max()/scale
        return self.uncompress(proj,offsets=offsets,locations=locations,N_u=256,N_v=256) 

    def __repr__(self): 
        s = "Static PET acquisition:  \n" 
        s = s+" - Time_start:                   %s \n"%millisec_to_min_sec(self.time_start)
        s = s+" - Time_end:                     %s \n"%millisec_to_min_sec(self.time_end)
        s = s+" - Duration:                     %s \n"%millisec_to_min_sec(self.time_end-self.time_start)
        s = s+" - N_counts:                     %d \n"%self.N_counts 
        s = s+" - N_locations:                  %d \n"%self.N_locations
        s = s+" - compression_ratio:            %d \n"%self.compression_ratio
        s = s+" - listmode_loss:                %d \n"%self.listmode_loss
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

    def _repr_html_(self):
        table_data = [['Time_start',millisec_to_min_sec(self.time_start)],
        ['Time_end',millisec_to_min_sec(self.time_end)],
        ['Duration',millisec_to_min_sec(self.time_end-self.time_start)],
        ['N_counts',pretty_print_large_number(self.N_counts)],
        ['N_locations',pretty_print_large_number(self.N_locations)],
        ['compression_ratio',print_percentage(self.compression_ratio)],
        ['listmode_loss',self.listmode_loss],
        ['Name',self.interface.name],['Manufacturer',self.interface.manufacturer],['Version',self.interface.version], ] 
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic_left')
        #table = ipy_table.set_column_style(0, color='lightBlue')
        table = ipy_table.set_global_style(float_format="%3.3f")        
        return table._repr_html_()

    def __del__(self):
        del(self.interface) 


        
        
        
class PET_Dynamic_Scan(): 
    """PET Dynamic Scan. """
    def __init__(self): 
        self.interface     = None                        # PET Scanner interface 
        self.binning       = Binning(DEFAULT_BINNING)    # PET detectors binning 
        self.static        = None                        # Static scan (no time binning)
        self._dynamic       = []                          # Sequence of static scans, one per time bin 
        self.time_bins     = []                          # Time binning (array, [ms])

        self.N_time_bins   = 0                           # Number of time bins         
        self.N_counts      = 0                           # Total number of counts 
        self.N_locations   = 0                           # Number of active locations of interaction
        self.compression_ratio = 0.0                     # Compression ratio of projection data 
        self.dynamic_inflation = 0.0                     # Inflation of projection data from static compressed to dynamic compressed 
        self.listmode_loss     = 0.0                     # Loss of data when events are binned from listmode 
        
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
        # define the ilang probabilistic model 
        self.ilang_model = PET_Dynamic_Poisson(self)  
        # construct a basic Directed Acyclical Graph
        #self.dag         = ProbabilisticGraphicalModel(['lambda','alpha','z']) 
        #self.dag.set_nodes_given(['counts','alpha'],True) 
        #self.dag.add_dependence(self.ilang_model,{'lambda':'lambda','alpha':'alpha','z':'z'}) 

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
        #  1 - see if the specified listmode data file exists 
        if data_filename!=None: 
            data_filename = data_filename.replace("/",os.path.sep).replace("\\",os.path.sep)          # cross platform compatibility 
            if not os.path.exists(data_filename): 
                raise FileNotFound("listmode data",data_filename)  
        #  2 - if the listmode data file is not specified, try with the name (and full path) contained in the listmode header
        data_filename      = hdr['name of data file']['value']
        data_filename = data_filename.replace("/",os.path.sep).replace("\\",os.path.sep)              # cross platform compatibility
        if not os.path.exists(data_filename): 
        #  3 - if it doesn't exist, look in the same path as the header file for the listmode data file with name specified in the listmode header file 
            data_filename = os.path.split(hdr_filename)[0]+os.path.sep+os.path.split(data_filename)[-1]  
            if not os.path.exists(data_filename): 
        #  4 - if it doesn't exist, look in the same path as the header file for the listmode data file with same name as the listmode header file, replacing the extension: ".l.hdr -> .l" 
                if hdr_filename.endswith(".l.hdr"): 
                    data_filename = hdr_filename.replace(".l.hdr",".l") 
                    if not os.path.exists(data_filename):     
                        raise FileNotFound("listmode data",data_filename)  
        #  5 - if it doesn't exist, look in the same path as the header file for the listmode data file with same name as the listmode header file, replacing the extension: ".hdr -> .l" 
                elif hdr_filename.endswith(".hdr"): 
                    data_filename = hdr_filename.replace(".hdr",".l") 
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

        progress_bar = ProgressBar()
        
        # Load the listmode data 
        R = self.interface.load_listmode(data_filename,n_packets,time_bins,self.binning,n_radial_bins, n_angles, n_sinograms, progress_bar.set_percentage) 
        progress_bar.set_percentage(100) 

        self.N_counts          = R['N_counts'] 
        self.N_locations       = R['N_locations']
        self.compression_ratio = R['compression_ratio']
        self.dynamic_inflation = R['dynamic_inflation']
        self.listmode_loss     = R['listmode_loss']
        self.N_time_bins  = R['N_time_bins'] 
        self.time_start   = R['time_start'] 
        self.time_end     = R['time_end'] 
        self.time_bins    = time_bins[0:self.N_time_bins+1]  #the actual time bins are less than the requested time bins, truncate time_bins 

        # Make list of PET_Static_Scan objects, one per bin
        self._dynamic = [] 
        for t in range(self.N_time_bins): 
            PET_t = PET_Static_Scan() 
            PET_t.set_interface(self.interface) 
            PET_t.set_binning(self.binning) 
            PET_t.load_static_measurement(t) 
            PET_t.scanner_detected = self.scanner_detected 
            # make list of static scans
            self._dynamic.append(PET_t) 
            # also make one attribut for each static scan
            setattr(self,"frame%d"%t,self._dynamic[t]) 

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
        return self 
        
    def load_static_measurement(self): 
        R = self.interface.get_measurement_static() 
        self.time_start               = R['time_start']
        self.time_end                 = R['time_end']
        self.N_counts                 = R['N_counts']    
        self.N_locations              = R['N_locations']
        self.compression_ratio        = R['compression_ratio']
        self.listmode_loss            = R['listmode_loss']    
        self._offsets                 = R['offsets']
        self._locations               = R['locations'] 
        self._static_measurement_data = R['counts'] 

    def get_static_measurement(self): 
        return (self._static_measurement_data,self._locations,self._offsets)

    def uncompressed_measurement(self): 
        uncompressed_measurement = self.uncompress(self._static_measurement_data) 
        return uncompressed_measurement 
               
    def uncompress(self, projection_data, offsets=None, locations=None, N_u=None, N_v=None): 
        if offsets==None:
            offsets=self._offsets
        if locations==None:
            locations=self._locations
        if N_u==None:
            N_u=self.binning.N_u
        if N_v==None:
            N_v=self.binning.N_v
        return self.interface.uncompress(offsets, projection_data, locations, N_u, N_v) 
               
    def display_measurements_in_browser(self): 
        return self.display_sequence(open_browser=True) 
        
    def display_measurements(self,open_browser=False): 
        #im_size = self.uncompressed_measurement().to_image().size
        #IM = PIL.new("RGB",(im_size[0],im_size[1]*self.N_time_bins))
        images = []
        for i in range(self.N_time_bins):
            im = self[i].uncompressed_measurement().to_image() 
            #IM.paste(im,(0,im_size[1]*i)) 
            images.append(im)
        d = DisplayNode() 
        return d.display('tipix',images,open_browser)
        #return d.display('image',IM.rotate(90),open_browser)
    
    def __repr__(self): 
        """Display information about Dynamic_PET_Scan"""
        s = "Dynamic PET acquisition:  \n" 
        s = s+" - N_time_bins:                  %d \n"%self.N_time_bins 
        s = s+" - Time_start:                   %s \n"%millisec_to_min_sec(self.time_start)
        s = s+" - Time_end:                     %s \n"%millisec_to_min_sec(self.time_end)
        s = s+" - N_counts:                     %d \n"%self.N_counts 
        s = s+" - N_locations:                  %d \n"%self.N_locations
        s = s+" - compression_ratio:            %d \n"%self.compression_ratio
        s = s+" - dynamic_inflation:            %d \n"%self.dynamic_inflation
        s = s+" - listmode_loss:                %d \n"%self.listmode_loss
        s = s+" - Mean time bin duration:       %d [sec] \n"%0 #FIXME 
        if self.interface != None: 
            s = s+" * Interface: \n" 
            s = s+"     - Name:                     %s \n"%self.interface.name 
            s = s+"     - Manufacturer:             %s \n"%self.interface.manufacturer 
            s = s+"     - Version:                  %s \n"%self.interface.version 
        if self.binning != None:
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

#    def _repr_html_(self):
#        return self.display_sequence()._repr_html_()

    def _repr_html_(self):
        table_data = [['N_time_bins',self.N_time_bins],
        ['Time_start',millisec_to_min_sec(self.time_start)],
        ['Time_end',millisec_to_min_sec(self.time_end)],
        ['Duration',millisec_to_min_sec(self.time_end-self.time_start)],
        ['N_counts',pretty_print_large_number(self.N_counts)],
        ['N_locations',pretty_print_large_number(self.N_locations)],
        ['compression_ratio',print_percentage(self.compression_ratio)],
        ['dynamic_inflation',self.dynamic_inflation],
        ['listmode_loss',self.listmode_loss], ]
        if self.interface: 
            table_data += [['Name',self.interface.name],['Manufacturer',self.interface.manufacturer],['Version',self.interface.version], ]
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic_left')
        #table = ipy_table.set_column_style(0, color='lightBlue')
        table = ipy_table.set_global_style(float_format="%3.3f")        
        return table._repr_html_()

    def __iter__(self): 
        """This method makes the object iterable. """
        return iter(self._dynamic)
    
    def __getitem__(self,i): 
        """This method makes the object addressable like a list. """
        return self._dynamic[i] 

    def __del__(self):
        """Delete interface when the object is deleted: interface needs explicit 
        deletion in order to manage C library memory deallocation  """
        self.interface.free_memory()

        
        
