from setuptools import setup, find_packages

setup(
    name='QEhandler',
    version='0.2',
    packages=find_packages(exclude=['*test*']),
    url='https://github.com/woosunjang/QEhandler',
    license='GPL-3.0',
    author='Woosun Jang',
    author_email='jin890@yonsei.ac.kr',
    description='',
    package_data={"pwscf": ["*.yaml"]},
    include_package_data=True
)
