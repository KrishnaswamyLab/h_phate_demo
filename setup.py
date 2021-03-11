import os
import sys
from setuptools import setup, find_packages

install_requires = [
    'numpy>=1.10.0',
    'scipy>=0.18.0',
    'pandas>=0.25',
    'scprep',
    'tasklogger',
    'graphtools',
    'phate',
    'matplotlib>=3.0',
    'ipywidgets',
    'plotly==3.10.0',
    'graph_coarsening @ git+https://github.com/loukasa/graph-coarsening',
]

package_name = "h_phate"

version_py = os.path.join(os.path.dirname(
    __file__), package_name, 'version.py')
version = open(version_py).read().strip().split(
    '=')[-1].replace('"', '').strip()

readme = open('README.md').read()

setup(name=package_name,
      version=version,
      description=package_name,
      author='Scott Gigante, Krishnaswamy Lab, Yale University',
      author_email='krishnaswamylab@gmail.com',
      packages=find_packages(),
      license='GNU General Public License Version 2',
      install_requires=install_requires,
      long_description=readme,
      url='https://github.com/KrishnaswamyLab/h_phate_demo',
      download_url="https://github.com/KrishnaswamyLab/h_phate_demo/archive/v{}.tar.gz".format(
          version),
      keywords=['big-data',
                'dimensionality-reduction',
                ],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Framework :: Jupyter',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering :: Visualization',
      ]
      )
