from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='PyCav',

    version='1.0.0a4',

    description='PyCav module',
    long_description=long_description,

    url='https://github.com/PyCav/PyCav-Module',

    author='PyCav team 2016',
    author_email=' ',

    license='BSD',

    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Education',
        'Topic :: Education',
        'Topic :: Scientific/Engineering :: Physics',

        'License :: OSI Approved :: BSD License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    keywords='physics simulation education',

    packages=['pycav'],

    install_requires=['numpy>=1.1',
    'matplotlib',
    'scipy'],

)
