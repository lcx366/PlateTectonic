import setuptools
from setuptools import setup 

setup(
    name='platemotion',
    version='0.1.8',
    description='A package to handle the tectonic plate motion',
    author='Chunxiao Li',
    author_email='lcx366@126.com',
    url='https://github.com/lcx366/PlateTectonic',
    license='MIT',
    long_description_content_type='text/markdown',
    long_description=open('README.md', 'rb').read().decode('utf-8'),
    keywords = ['tectonic plate','plate motion','nnr-morvel56','gsrmv2.1','itrf2014'],
    python_requires = '>=3.8',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: MIT License',
        ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'astropy',
        'pandas',
        'xarray',
        'sphericalpolygon',
        'tqdm',
        'colorama',
        'ftplib',
        'requests',
        'pkg_resources',
        'pathlib',
        'gzip',
        'zipfile'
        ],
    extras_require={
        'files': ['os'],
        'other': ['datetime','time']
        }      
    )
