echo off
IF "%1"=="-install" (
	call pip install -r list.txt > install.txt
	del install.txt
)
python Synapse3D.py
del /S *.pyc