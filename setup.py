"""
PyPI pip package following the material by Stephen Hudson found here:
https://betterscientificsoftware.github.io/python-for-hpc/tutorials/python-pypi-packaging/
with release download_url explanation from Joel Barmettler:
https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56
"""

from setuptools import setup, find_packages

setup(
    name='opqua',
    version='v1.2.1',
    description='An epidemiological modeling framework for population ' \
        + 'genetics and evolution.',
    long_description='Opqua is an epidemiological modeling framework for ' \
        + 'population genetics and evolution. Opqua stochastically simulates ' \
        + ' pathogens with specific, evolving genotypes spread through ' \
        + ' populations of hosts that can have specific immune profiles. \n\n' \
        + 'Opqua is a useful tool to test out scenarios, explore hypotheses, ' \
        + 'and make predictions about the relationship between pathogen ' \
        + 'evolution and epidemiology. \n\n Visit ' \
        + 'github.com/pablocarderam/opqua for more information.',
    description_file='README.md',
    url='https://github.com/pablocarderam/opqua',
    download_url='https://github.com/pablocarderam/opqua/archive/v1.2.1.tar.gz',
    author='Pablo Cardenas',
    author_email='pablocarderam@gmail.com',
    keywords=['epidemiology','evolution','biology'],
    license='MIT',
    packages=find_packages(),
    install_requires=['joblib',
                      'textdistance',
                      'numpy',
                      'pandas',
                      'scipy',
                      'matplotlib',
                      'seaborn',
                      'pyvis',
                      ],

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
    ],
    python_requires='>=3.9',
    include_package_data=True, # Needed to include non-Python files, e.g. CSVs
)
