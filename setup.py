"""
Preprocessing ATAC and scATAC Data
"""
from setuptools import find_packages, setup

dependencies = ['click', 'Numpy', 'editdistance', 'pytest', 'python-levenshtein', 'snakemake', 'biopython', 'optparse-pretty', 'regex', 'PyYAML']

setup(
    name='proatac',
    version='0.2.0',
    url='https://github.com/buenrostrolab/parkour',
    license='BSD',
    author='Caleb Lareau and Jason Buenrostro',
    author_email='caleblareau@g.harvard.edu',
    description='Processing and quality control of ATAC data.',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'proatac = proatac.cli:main',
        ],
    },
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'Development Status :: 1 - Planning',
         'Development Status :: 2 - Pre-Alpha',
        # 'Development Status :: 3 - Alpha',
        # 'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        # 'Development Status :: 6 - Mature',
        # 'Development Status :: 7 - Inactive',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
