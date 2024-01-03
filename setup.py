# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='metaorf',
    version='0.0.1',
    description="""MetaORF is an sORF focused Ribo-Seq meta-caller""",
    author='An Zheng',
    author_email='an@veliatx.com',
    classifiers=[
        'Programming Language :: Python :: 3.7',
    ],
    keywords='microproteins',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.17.2',
        'pandas',
        'scipy>=1.3.1',
        'pytest>=4.6.6',
        'six>=1.12.0',
        'configparser>=4.0.2',
        'click',
    ],
    entry_points = {
        'console_scripts': ['metaorf=metaorf.pipelines:main'],
    }
)