#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2023 Stephen Wasilewski, EPFL
# =======================================================================
# This program is free software: you can redistribute it and/or
# modify it under the terms of theGNU Lesser General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
# =======================================================================

"""The setup script."""
from setuptools import find_packages, setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['raytraverse', 'tifffile']

setup(
    author="Stephen Wasilewski",
    author_email='stephanwaz@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: GNU Lesser General Public License 3.0 (LGPL 3.0)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="merging and managing HDR sequence images",
    python_requires=">=3.7",
    entry_points={
        'console_scripts': ['pylinearhdr=pylinearhdr.cli:main'],
    },
    install_requires=requirements,
    license="GNU Lesser General Public License 3.0 (LGPL 3.0)",
    long_description=readme,
    keywords='pylinearhdr',
    name='pylinearhdr',
    packages=find_packages(),
    url='https://github.com/stephanwaz/linearhdr',
    version='0.2.0',
    zip_safe=False,
)