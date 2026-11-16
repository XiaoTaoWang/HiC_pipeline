# Created on Wed Dec 24 15:26:26 2014

# Author: XiaoTao Wang

"""
Setup script for hicpeaks.

This is a free software under GPLv3. Therefore, you can modify, redistribute
or even mix it with other GPL-compatible codes. See the file LICENSE
included with the distribution for more details.

"""
import os, sys, runHiC, glob
import setuptools

def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname), encoding='utf-8') as source:
        return source.read()

if (sys.version_info.major!=3) or (sys.version_info.minor<10):
    print('PYTHON 3.10+ IS REQUIRED. YOU ARE CURRENTLY USING PYTHON {}'.format(sys.version.split()[0]))
    sys.exit(2)

# Guarantee Unix Format
for src in glob.glob('scripts/*'):
    with open(src, 'r', encoding='utf-8') as source:
        text = source.read()
    with open(src, 'w', encoding='utf-8', newline='\n') as destination:
        destination.write(text)

setuptools.setup(
    name = 'runhic',
    version = runHiC.__version__,
    author = runHiC.__author__,
    author_email = 'wangxiaotao686@gmail.com',
    url = 'https://github.com/XiaoTaoWang/HiC_pipeline',
    description = 'A easy-to-use Hi-C processing software supporting distributed computation',
    keywords = 'Hi-C Arima ICE cooler pairs bioinformatics pipeline',
    long_description = read('README.rst'),
    long_description_content_type='text/x-rst',
    scripts = glob.glob('scripts/*'),
    packages = setuptools.find_packages(),
    python_requires = '>=3.10',
    install_requires = ['cooler', 'pypairix'],
    classifiers = [
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ]
    )
