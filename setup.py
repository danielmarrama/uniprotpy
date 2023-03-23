from setuptools import setup
import re
import os


with open('uniprotpy/version.py', 'r') as f:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        f.read(),
        re.MULTILINE).group(1)

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='uniprotpy',
    version=version,
    description='A Python library that interfaces with UniProt data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/danielmarrama/uniprotpy',
    author='Daniel Marrama',
    author_email='danielmarrama@gmail.com',
    packages=['uniprotpy'],
    install_requires=['pandas>=1.1',
                      'biopython>=1.5'],
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'uniprotpy = uniprotpy.shell:run'
        ],
    },
)