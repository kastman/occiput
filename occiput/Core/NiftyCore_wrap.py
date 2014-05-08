
# occiput  
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# Dec. 2014, Boston, MA 
# April. 2014, Boston, MA


# FIXME: perhaps get rid of this file and place code in occiput.Transformation (mostly), although as it is now, all imports from 
# NiftyCore are in one place (here). This helps the joint development. 


import occiput as __occiput
import numpy as __np
from NiftyCore.NiftyRec import INTERPOLATION_LINEAR, INTERPOLATION_POINT
try:
    from NiftyCore.NiftyRec import TR_resample_grid as             __TR_resample_grid 
    from NiftyCore.NiftyRec import TR_grid_from_box_and_affine as  __TR_grid_from_box_and_affine 
    from NiftyCore.NiftyRec import TR_transform_grid as            __TR_transform_grid 
    from NiftyCore.NiftyRec import INTERPOLATION_LINEAR, INTERPOLATION_POINT
    from NiftyCore.NiftyRec import PET_project_compressed, PET_backproject_compressed 
    from NiftyCore.NiftyRec import SPECT_project_parallelholes, SPECT_backproject_parallelholes 
    has_NiftyCore = True
except: 
    print "NiftyCore could not be loaded: it will not be possible to reconstruct the PET data. "
    has_NiftyCore = False
    PET_project_compressed = None
    PET_backproject_compressed = None
    SPECT_project_parallelholes = None
    SPECT_backproject_parallelholes = None 


## Transformation 

def __make_grid(data,space,ndim):
    if ndim==3: 
        return __occiput.Core.Grid3D(data,space)
    else: 
        return __occiput.Core.GridND(data,space)



def transform_grid(grid, affine_from_grid): 
    # 1) verify if the spaces of the affine map and of the grid are compatible: 
    if not affine_from_grid.can_left_multiply(grid): 
        print "Affine transformation not compatible with grid. " 
        # FIXME: raise error, or warning, depending on a global setting 
    # 2) transform 
    transformed = __TR_transform_grid( grid.data, affine_from_grid.data ) 
    # 3) instantiate a new grid 
    grid = __make_grid(transformed, affine_from_grid.map_to, transformed.ndim-1)
    return grid



def grid_from_box_and_affine(min_coords, max_coords, n_points, affine=None, space="world"): 
    #FIXME: optionally use the affine to transform min_coords and max_coords 
    data = __TR_grid_from_box_and_affine(min_coords,max_coords,n_points)
    ndim = data.ndim-1
    grid = __make_grid(data, space, ndim)
    return grid



def resample_image_on_grid(image, grid, affine_grid_to_world=None, verify_mapping=True, background=0.0, use_gpu=1): 
    if verify_mapping:
        # check if the image, the grid and the affine mapping are compatible: 
        # 1) if affine_grid_to_world is not defined, verify if image and grid are compatible
        if affine_grid_to_world == None: 
            if not image.affine.can_inverse_left_multiply(grid): 
                print "grid and image are not compatible. "
                #FIXME: raise exception
                return       
        # 2) if affine_grid_to_world is defined, verify if image, grid and affine_grid_to_world are compatible 
        else: 
            # a) check mapping from grid to world 
            if not affine_grid_to_world.can_left_multiply(grid): 
                print "grid and affine_grid_to_world are not compatible. "
                #FIXME: raise exception
                return    
            # b) check mapping from image to world 
            if not image.affine.can_inverse_left_multiply(affine_grid_to_world):
                print "image and affine_grid_to_world are not compatible. "
                #FIXME: raise exception
                return
    # compute affine: 
    if affine_grid_to_world == None:
        affine = image.affine
    else: 
         affine = affine_grid_to_world.left_multiply(image.affine)
    # decide sampling mode
    if image.is_mask(): 
        interpolation_mode = INTERPOLATION_POINT
    else:
        interpolation_mode = INTERPOLATION_LINEAR 
    # resample: 
    resampled_data = __TR_resample_grid( image.data, grid.data, affine.data, background, use_gpu, interpolation_mode )
    return resampled_data 
    
    

    
    
