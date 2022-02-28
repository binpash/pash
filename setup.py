import setuptools

setuptools.setup(
    name='pash',
    version='0.0.1',    
    description='Parallelize POSIX shell scripts.',
    url='https://github.com/binpash/pash',
    license='MIT',
    package_dir={'pash': '.'},
    packages=setuptools.find_packages('.'),
    install_requires=[],
    classifiers=[],
)
