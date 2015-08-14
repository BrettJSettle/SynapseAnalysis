echo off
IF "%1"=="-install" (
	call pip install -r list.txt > install.txt
	del install.txt
)
python Synapse.py
del /S *.pyc
del /S __pycache__
