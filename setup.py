# Created on Wed Dec 24 15:26:26 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University


import os, sys, lib
from distutils.core import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if (sys.version_info.major != 2) or (sys.version_info.minor < 6):
    print 'PYTHON VERSION MUST BE 2.6 or 2.7. YOU ARE CURRENTLY USING PYTHON ' + sys.version
    sys.exit(2)

# Guarantee Unix Format
text = open('scripts/runHiC', 'rb').read().replace('\r\n', '\n')
open('scripts/runHiC', 'wb').write(text)

setup(
    name = 'runHiC',
    version = lib.__version__,
    author = lib.__author__,
    author_email = 'wangxiaotao868@163.com',
    url = 'https://github.com/XiaoTaoWang/HiC_pipeline',
    description = 'A easy-to-use Hi-C processing software based on hiclib',
    keywords = 'Hi-C HiC ICE Contact',
    package_dir = {'runHiC':'lib'},
    packages = ['runHiC'],
    scripts = ['scripts/runHiC'],
    long_description = read('README.rst'),
    classifiers = [
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Operating System :: POSIX',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ]
    )
