
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston


__all__ = ['uniform_sphere',]



from NiftyCore.NiftyRec import ET_spherical_phantom



def uniform_sphere(voxels,size,center,radius,inner_value,outer_value): 
    """Create volume (3D numpy array) with uniform value within a spherical region. """
    return ET_spherical_phantom(voxels,size,center,radius,inner_value,outer_value) 


