pycoordtrans
============

A python2 binding to coordtrans, which is written in C++ and provides coordinate conversion between bd09ll, gcj02 and wgs84.

Build
--------
Clone the repository and just run the following scripts. You will get a python module.

```Shell
git clone https://github.com/zxteloiv/pycoordtrans.git
cd pycoordtrans
python2 setup.py build_ext --inplace
```
Usage
---------
After you build this extension, you will have a python module named _coordtrans.
The module is the same as other python modules except that it is written in C++.
The compiled module is a linux shared library. All you get is a _coordtrans.so file residing at the current directory.

Run it as the following:
```Python
from _coordtrans import coordtrans
y, x = 39.981839,116.306411 # x, y are longitude and latitude respectively by convention
(x, y) = coordtrans('gcj02', 'bd09ll', x, y) # yielding a tuple (116.31290580529193, 39.98791010079803)
```
Handful tool
---------
To check the coordinates in commandline, you can also compile the main.cpp as a traditional C++ program.

```Shell
g++ main.cpp coordtrans.cpp -o coordtrans
./coordtrans gcj02 bd09ll 116.306411 39.981839
```

