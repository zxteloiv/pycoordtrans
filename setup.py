from distutils.core import setup, Extension

setup(
    ext_modules=[Extension("_coordtrans", ["_coordtrans.cpp", "coordtrans.cpp"])]
)
