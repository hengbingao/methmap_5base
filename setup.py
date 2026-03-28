from setuptools import setup, find_packages

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='methmap',
    version='1.3.5',
    description='Methylation alignment tool for 5-base sequencing chemistry (5mC -> T)',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='MethMap',
    url='https://github.com/YOUR_USERNAME/methmap',
    packages=find_packages(exclude=['tests*']),
    python_requires='>=3.7',
    install_requires=[],
    extras_require={
        'test': ['pytest>=6.0'],
    },
    entry_points={
        'console_scripts': [
            'methmap=methmap.__main__:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Intended Audience :: Science/Research',
    ],
    keywords='methylation sequencing epigenetics 5mC bisulfite-free',
)
