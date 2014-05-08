
# occiput - Computational engine for volumetric imaging 
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# Dec 2013, Boston 


from setuptools import setup, Extension, Library
from glob import glob 



setup(
    name='occiput',
    version='0.1.0',
    author='Stefano Pedemonte',
    author_email='stefano.pedemonte@gmail.com',
    packages=['occiput', 
              'occiput.test', 
              'occiput.notebooks',
              'occiput.Core',
              'occiput.Reconstruction', 
              'occiput.Reconstruction.PET', 
              'occiput.Reconstruction.SPECT', 
              'occiput.Reconstruction.CT', 
              'occiput.Transformation', 
              'occiput.Registration',
              'occiput.Registration.Affine',
              'occiput.Registration.TranslationRotation',
              'occiput.Classification', 
              'occiput.DataSources', 
              'occiput.DataSources.Synthetic', 
              'occiput.DataSources.FileSources', 
              'occiput.Visualization',
              ],
    scripts=[],
    url='http://tomographylab.scienceontheweb.net',
    license='LICENSE.txt',
    description='Tomographic Inference.',
    long_description=open('README.txt').read(),
    keywords = ["PET", "SPECT", "emission tomography", "transmission tomography", "tomographic reconstruction"],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
                     ],
    install_requires=[
        "numpy >= 1.6.0", 
        "simplewrap >= 0.1.0", 
        "interfile >= 0.1.0", 
        "petlink >= 0.1.0", 
        "DisplayNode >= 0.1.0", 
        "ilang >= 0.1.0", 
        "ipy_table >= 1.11.0",
        "pydicom >= 0.9.0",
        "nibabel >= 1.0.0",  
    ], 
) 


