#!/bin/bash 


apt install upx-ucl 
apt install python3-pysam 
python3 -m pip install pyinstaller && pip install -r requirements.txt &&  python3 -m PyInstaller toblerone.py --name Toblerone --onefile --hidden-import pysam.libctabixproxies --add-data "bin:bin" --add-data "assets:assets" --icon icon.icns 



