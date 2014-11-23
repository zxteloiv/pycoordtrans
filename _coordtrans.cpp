#include <python2.7/Python.h>
#include "coordtrans.h"

static char module_docstring[] = 
    "This module provides an interface for coordtrans in C.";

static char coordtrans_docstring[] = 
    "Convert coordinates among bd09ll, wgs84, gcj02";

static PyObject *_coordtrans_coordtrans(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"coordtrans", _coordtrans_coordtrans, METH_VARARGS, coordtrans_docstring},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_coordtrans(void) {
    PyObject *m = Py_InitModule3("_coordtrans", module_methods, module_docstring);
    if (m == NULL)
        return;
};

static PyObject *_coordtrans_coordtrans(PyObject *self, PyObject *args) {
    char *from = NULL;
    char *to = NULL;
    double lng = 0;
    double lat = 0;
    double newlng = 0;
    double newlat = 0;
    int ret = -1; // return value
    if (!PyArg_ParseTuple(args, "ssdd", &from, &to, &lng, &lat)) {
        return NULL;
    }

    ret = coordtrans(from, to, lng, lat, newlng, newlat);

    if (0 != ret) {
        PyErr_SetString(PyExc_RuntimeError, "impossible string");
        return NULL;
    }

    PyObject *val = Py_BuildValue("dd", newlng, newlat);
    return val;
}
