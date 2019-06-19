# VAFchecker
## Installation
VAFchecker needs a working installation of python 3.
### Clone repository
Clone the VAFchecker repository
```
> git clone git@github.com:ToolsVanBox/VAFchecker.git
```
### Virtual environment
Create a virtual environment
```
> virtualenv venv_3.6 -p /hpc/local/CentOS7/common/lang/python/3.6.1/bin/python
```
Activate the virtual environment
```
> . venv_3.6/bin/activate
```
### Requirements
VAFchecker requires the following modules:
* pyvcf
* pysam
* argparse
* multiprocessing
* queue

Install the required modules:
```
> pip install -r requirements.txt
```
