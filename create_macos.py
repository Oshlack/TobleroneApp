#!/bin/env python3 

"""
Usage:
	python3 create-macOS-app.py py2app
"""

from setuptools import setup

APP = ['toblerone.py']
DATA_FILES = [('', ['assets','bin'])]
OPTIONS = {
	'strip': True,
	'iconfile': 'icon.icns',
	'packages': 'pysam'
}

setup(
	name='Toblerone',
	app=APP,
	data_files=DATA_FILES,
	options={'py2app': OPTIONS},
	setup_requires=['py2app'],
)

