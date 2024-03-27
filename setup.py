from setuptools import setup

setup(
    name='pyfiremodels',
    version='0.0.0',    
    description='This software is part of a python library to assist in calculation of empirical and analytical relationships in fire safety engineering.',
    url='https://github.com/johodges/pyfiremodels',
    author='Jonathan Hodges',
    author_email='johodges@vt.edu',
    license='MIT',
    packages=['pyfiremodels'],
    
    include_package_data=True,
    install_requires=[
                        'matplotlib>=3.0',
                        'numpy>=1.17',
                        'pandas>=0.25',
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: OS Independent',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
