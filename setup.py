import versioneer
from setuptools import setup, find_packages, Extension
from codecs import open
from os import path
have_cython = False
try:
    from Cython.Distutils import build_ext
    have_cython = True
except ImportError:
    from distutils.command.build_ext import build_ext

mechanics = None
if have_cython:
    mechanics = Extension('pycav.mechanics', ['pycav/mechanics.pyx'])
else:
    mechanics = Extension('pycav.mechanics', ['pycav/mechanics.c'])

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

cmds = versioneer.get_cmdclass()
cmds["build_ext"] = build_ext

setup(
    version=versioneer.get_version(),
    cmdclass=cmds,
    name='PyCav',
   
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

    ext_modules=[mechanics],
    packages=['pycav'],

    install_requires=['numpy>=1.1',
    'matplotlib',
    'scipy'],

)
