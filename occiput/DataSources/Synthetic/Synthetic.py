
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston


__all__ = ['uniform_sphere','uniform_cylinder','uniform_spheres_ring']



from occiput.Core import RigidVolume 
from NiftyCore.NiftyRec import ET_spherical_phantom, ET_cylindrical_phantom, ET_spheres_ring_phantom



def uniform_sphere(voxels,size,center,radius,inner_value,outer_value): 
    """Create volume (3D numpy array) with uniform value within a spherical region. """
    return RigidVolume(ET_spherical_phantom(voxels,size,center,radius,inner_value,outer_value)) 

def uniform_cylinder(voxels,size,center,radius,length,axis,inner_value,outer_value): 
    """Create volume (3D numpy array) with uniform value within a spherical region. """
    return RigidVolume(ET_cylindrical_phantom(voxels,size,center,radius,length,axis,inner_value,outer_value)) 
    
def uniform_spheres_ring(voxels,size,center,ring_radius,min_sphere_radius,max_sphere_radius,N_spheres=6,inner_value=1.0,outer_value=0.0,taper=0,axis=0): 
    """Create volume (3D numpy array) with uniform value within a spherical region. """
    return RigidVolume(ET_spheres_ring_phantom(voxels,size,center,ring_radius,min_sphere_radius,max_sphere_radius,N_spheres,inner_value,outer_value,taper,axis))

