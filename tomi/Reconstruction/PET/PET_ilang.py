
# tomi - Tomographic Inference 
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# Dec. 2013, Boston, MA



import ilang
import ilang.Models
from ilang.Models import Model

__all__ = ['PET_Static_Poisson','PET_Dynamic_Poisson']



class PET_Static_Poisson(Model): 
    variables = {'lambda':'continuous','alpha':'continuous','z':'discrete'} 
    dependencies = [['lambda','z','directed'],['alpha','z','directed']]

    def __init__(self, PET_scan, name=None):
        if name == None:  
            name = self.__class__.__name__
        Model.__init__(self, name) 
        # PET scan
        self.PET_scan = PET_scan    
        # variables         
        self._lambda = None
        self._alpha  = None 
        self._roi    = None 
        
    def set_PET_scan(self, PET_scan): 
        self.PET_scan = PET_scan 

    def init(self): 
        pass 

    def log_conditional_probability_lambda(self): 
        return 0 

    def log_conditional_probability_gradient_lambda(self): 
        projection_data = self.PET_scan.get_projection(self._lambda, self._roi, self._alpha) 
        gradient        = self.PET_scan.get_backprojection(self._measurement_data/projection_data+EPS, roi, self._alpha) 
        return gradient 

    def log_conditional_probability_alpha(self): 
        return 0 

    def log_conditional_probability_gradient_alpha(self): 
        return numpy.zeros([100,1])
        
    def sample_conditional_probability_z(self): 
        return 0    



class PET_Dynamic_Poisson(Model): 
    variables = {'lambda':'continuous','alpha':'continuous','roi_1':'continuous','roi_1':'continuous','z_1':'discrete','z_2':'discrete'} 
    dependencies = [['lambda','z_1','directed'],['lambda','z_2','directed'],['alpha','z_1','directed'],['alpha','z_2','directed'],['roi_1','z_1','directed'],['roi_2','z_2','directed']] 

    def __init__(self, PET_scan, name=None): 
        if name == None:  
            name = self.__class__.__name__
        Model.__init__(self, name) 
        # PET scan object: 
        self.PET_scan = PET_scan 
        self.N_time_bins = self.PET_scan.N_time_bins 
        # Variables and dependencies: 
        self._lambda = None 
        self._alpha  = None 
        self._rois   = None 
        self.variables = {'lambda':'continuous', 'alpha':'continuous'} 
        self.dependencies = []
        for t in range(self.N_time_bins): 
            var_name_counts = 'z_'+str(t)
            var_name_roi    = 'roi_'+str(t)
            self.variables[var_name_counts]='discrete'
            self.variables[var_name_roi]='continuous'
            self.dependencies.append(['lambda',var_name_counts,'directed'])
            self.dependencies.append(['alpha',var_name_counts,'directed']) 
            self.dependencies.append([var_name_roi,var_name_counts,'directed']) 

    def set_PET_scan(self, PET_scan): 
        self.PET_scan = PET_scan 

    def init(self): 
        pass 

    def log_conditional_probability_lambda(self): 
        return 0 

    def log_conditional_probability_gradient_lambda(self): 
        return numpy.zeros([100,1])

    def log_conditional_probability_alpha(self): 
        return 0 

    def log_conditional_probability_gradient_alpha(self): 
        return numpy.zeros([100,1])
        
    def sample_conditional_probability_z(self): 
        return 0    





        