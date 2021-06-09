import setuptools

setuptools.setup(
    name='batchslopes',
    version='0.0.7',
    author='Ulf Liebal',
    author_email='ulf.liebal@rwth-aachen.de',
    description='batchslopes is a tool to identify the growth rate in growth profiler data. The data is partitioned into decreasing set sizes until the highest correlation coefficient of exponential growth is detected.',
    keywords='optical density, growth curve, growth rate, exponential growth, growth profiler',
    url='https://git.rwth-aachen.de/ulf.liebal/batchslopes.git',
    package_dir={'': 'src'},
    packages=setuptools.find_packages(where='src'),
    python_requires='>=3.6',
    install_requires=[
        'numpy>=1.18.0'
    ],
)
