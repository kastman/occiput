
# tomi - Tomographic Inference 
# Stefano Pedemonte
# Harvard University
# Dec 2013, Boston 


from setuptools import setup, Extension, Library
from glob import glob 



setup(
    name='tomi',
    version='0.1.0',
    author='Stefano Pedemonte',
    author_email='stefano.pedemonte@gmail.com',
    packages=['tomi', 
              'tomi.test', 
              'tomi.Reconstruction', 
              'tomi.Reconstruction.PET', 
              'tomi.Reconstruction.SPECT', 
              'tomi.Reconstruction.CT', 
              'tomi.Transformation', 
              'tomi.Classification',                  
              ],
    scripts=[],
    url='http://tomi.com/',
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
        "numpy >= 1.7.1", 
        "simplewrap >= 0.1.0", 
        "interfile >= 0.1.0", 
        "petlink >= 0.1.0", 
        "DisplayNode >= 0.1.0", 
        "ilang >= 0.1.0", 
        "ipy_table >= 1.11.0",  
    ], 
) 


