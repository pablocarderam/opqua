"""
PyPI pip package following the material by Stephen Hudson found here:
https://betterscientificsoftware.github.io/python-for-hpc/tutorials/python-pypi-packaging/
with release explanation from Joel Barmettler:
https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56
"""

from setuptools import setup

setup(
    name='opqua',
    version='v0.1.1',
    description='An epidemiological modeling framework for population ' \
        + 'genetics and evolution.',
    long_description='Opqua is an epidemiological modeling framework for ' \
        + 'population genetics and evolution. Opqua stochastically simulates ' \
        + ' pathogens with specific, evolving genotypes spread through ' \
        + ' populations of hosts that can have specific immune profiles. ' \
        + 'Opqua is a useful tool to test out scenarios, explore hypotheses, ' \
        + 'and make predictions about the relationship between pathogen ' \
        + 'evolution and epidemiology.',
    url='https://github.com/pablocarderam/opqua',
    download_url = 'https://github.com/pablocarderam/opqua/archive/v0.1.1.tar.gz',
    author='Pablo Cardenas',
    author_email='pablocarderam@gmail.com',
    keywords=['epidemiology','evolution','biology'],
    license='MIT',
    packages=['opqua'],
    install_requires=['joblib',
                      'numpy',
                      'pandas',
                      'matplotlib',
                      'seaborn',
                      ],

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
