import setuptools

setuptools.setup(
    name='batchslopes',
    version='0.0.1',
    author='Ulf Liebal',
    author_email='ulf.liebal@rwth-aachen.de',
    description='batchslopes is a tool to identify the growth rate in growth profiler data. The data is partitioned into decreasing set sizes until the highest correlation coefficient of exponential growth is detected.',
    keywords='optical density, growth curve, growth rate, exponential growth, growth profiler',
    url='https://git.rwth-aachen.de/ulf.liebal/batchslopes.git',
    package_dir={'': 'src'},
    packages=setuptools.find_packages(where='src'),
    classifiers=[
        'Development Status :: 1 - Development/Unstable',
        'Intended Audience :: Biotechnology Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy>=1.18.0',
        'pandas>=1.1.3'
    ],
)
