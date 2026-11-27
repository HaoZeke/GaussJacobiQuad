// BEGIN_HEADER
// -----------------------------------------------------------------------------
// Gauss-Jacobi Quadrature Implementation
// Authors: Rohit Goswami <rgoswami[at]ieee.org>
// Source: GaussJacobiQuad Library
// License: MIT
// GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
// -----------------------------------------------------------------------------
// END_HEADER
/*
 * CPython extension for GaussJacobiQuad.
 *
 * Calls ISO_C_BINDING C ABI: gauss_jacobi_rule_c / gjp_status_string
 * (same objects as libgjp_cinterp, linked into this module).
 *
 * Dual packaging (both are intentional):
 *
 * 1. Stable ABI (default): Py_LIMITED_API 0x03090000 → one cp39-abi3 wheel
 *    for CPython 3.9+ GIL builds.
 * 2. Free-threaded (GJP_FULL_API / freethreaded interpreter): full C API +
 *    Py_mod_gil = Py_MOD_GIL_NOT_USED → cp313t (and later) wheels. Free-threaded
 *    CPython is a different ABI and cannot load abi3.
 *
 * Common: multi-phase PEP 489 init, module state, heap GaussJacobiError.
 */
#define PY_SSIZE_T_CLEAN
/* Limited API only for non–free-threaded / non–full-API builds (before Python.h). */
#if !defined(GJP_FULL_API) && !defined(Py_LIMITED_API)
#  define Py_LIMITED_API 0x03090000
#endif
#include <Python.h>

#include "GaussJacobiQuadCInterp.h"

#include <limits.h>
#include <string.h>

/* ---- module state ---- */
typedef struct {
    PyObject *error; /* heap exception type */
} gjp_module_state;

static gjp_module_state *
gjp_state(PyObject *module)
{
    return (gjp_module_state *)PyModule_GetState(module);
}

static int
gjp_traverse(PyObject *module, visitproc visit, void *arg)
{
    gjp_module_state *st = gjp_state(module);
    Py_VISIT(st->error);
    return 0;
}

static int
gjp_clear(PyObject *module)
{
    gjp_module_state *st = gjp_state(module);
    Py_CLEAR(st->error);
    return 0;
}

static void
gjp_free(void *module)
{
    gjp_clear((PyObject *)module);
}

/* ---- helpers ---- */
static PyObject *
gjp_set_error(PyObject *module, int status, const char *method)
{
    gjp_module_state *st = gjp_state(module);
    const char *msg = gjp_status_string(status);
    if (!msg)
        msg = "unknown status";
    if (method && method[0])
        return PyErr_Format(st->error, "%s (status=%d, method=%s)", msg, status, method);
    return PyErr_Format(st->error, "%s (status=%d)", msg, status);
}

/*
 * rule(npts, alpha, beta, method=None) -> (nodes: list[float], weights: list[float])
 */
static PyObject *
gjp_rule(PyObject *module, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"npts", "alpha", "beta", "method", NULL};
    Py_ssize_t npts_ss = 0;
    double alpha = 0.0, beta = 0.0;
    const char *method = NULL;
    PyObject *method_obj = Py_None;
    PyObject *method_bytes = NULL; /* keep alive for Limited API UTF-8 */

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ndd|O:rule", kwlist,
                                     &npts_ss, &alpha, &beta, &method_obj))
        return NULL;

    if (npts_ss <= 0 || npts_ss > (Py_ssize_t)INT_MAX) {
        return gjp_set_error(module, GJP_ERR_NPTS, method);
    }

    if (method_obj == Py_None) {
        method = "";
    }
    else if (PyUnicode_Check(method_obj)) {
#if defined(Py_LIMITED_API) && Py_LIMITED_API < 0x030A0000
        /* PyUnicode_AsUTF8 is not in the 3.9 Limited API surface */
        method_bytes = PyUnicode_AsUTF8String(method_obj);
        if (!method_bytes)
            return NULL;
        method = PyBytes_AsString(method_bytes);
        if (!method) {
            Py_DECREF(method_bytes);
            return NULL;
        }
#else
        method = PyUnicode_AsUTF8(method_obj);
        if (!method)
            return NULL;
#endif
        if (strcmp(method, "auto") == 0)
            method = "";
    }
    else {
        PyErr_SetString(PyExc_TypeError, "method must be str or None");
        return NULL;
    }

    int npts = (int)npts_ss;
    double *x = (double *)PyMem_Malloc((size_t)npts * sizeof(double));
    double *w = (double *)PyMem_Malloc((size_t)npts * sizeof(double));
    if (!x || !w) {
        PyMem_Free(x);
        PyMem_Free(w);
        Py_XDECREF(method_bytes);
        return PyErr_NoMemory();
    }

    int status;
    Py_BEGIN_ALLOW_THREADS
    status = gauss_jacobi_rule_c(npts, alpha, beta, x, w, method ? method : "");
    Py_END_ALLOW_THREADS

    if (status != GJP_OK) {
        PyMem_Free(x);
        PyMem_Free(w);
        PyObject *err = gjp_set_error(module, status, method && method[0] ? method : "auto");
        Py_XDECREF(method_bytes);
        return err;
    }

    PyObject *x_list = PyList_New(npts);
    PyObject *w_list = PyList_New(npts);
    if (!x_list || !w_list) {
        Py_XDECREF(x_list);
        Py_XDECREF(w_list);
        PyMem_Free(x);
        PyMem_Free(w);
        Py_XDECREF(method_bytes);
        return NULL;
    }

    for (int i = 0; i < npts; i++) {
        PyObject *xv = PyFloat_FromDouble(x[i]);
        PyObject *wv = PyFloat_FromDouble(w[i]);
        if (!xv || !wv) {
            Py_XDECREF(xv);
            Py_XDECREF(wv);
            Py_DECREF(x_list);
            Py_DECREF(w_list);
            PyMem_Free(x);
            PyMem_Free(w);
            Py_XDECREF(method_bytes);
            return NULL;
        }
        if (PyList_SetItem(x_list, i, xv) < 0 || PyList_SetItem(w_list, i, wv) < 0) {
            Py_DECREF(x_list);
            Py_DECREF(w_list);
            PyMem_Free(x);
            PyMem_Free(w);
            Py_XDECREF(method_bytes);
            return NULL;
        }
    }

    PyMem_Free(x);
    PyMem_Free(w);
    Py_XDECREF(method_bytes);

    PyObject *pair = PyTuple_Pack(2, x_list, w_list);
    Py_DECREF(x_list);
    Py_DECREF(w_list);
    return pair;
}

/*
 * scale(nodes, weights, alpha, beta, a, b) -> (nodes, weights) on [a, b].
 * Weights gain ((b-a)/2)**(alpha+beta+1) for (b-t)**alpha * (t-a)**beta.
 */
static PyObject *
gjp_scale(PyObject *module, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"nodes", "weights", "alpha", "beta", "a", "b", NULL};
    PyObject *nodes_obj = NULL;
    PyObject *weights_obj = NULL;
    double alpha = 0.0, beta = 0.0, a = 0.0, b = 0.0;
    Py_ssize_t n, nw, i;
    double *x = NULL;
    double *w = NULL;
    int status;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOdddd:scale", kwlist, &nodes_obj, &weights_obj,
                                     &alpha, &beta, &a, &b))
        return NULL;

    n = PySequence_Size(nodes_obj);
    nw = PySequence_Size(weights_obj);
    if (n < 0 || nw < 0) {
        PyErr_SetString(PyExc_TypeError, "nodes and weights must be sequences");
        return NULL;
    }
    if (n != nw) {
        PyErr_SetString(PyExc_ValueError, "nodes and weights must have the same length");
        return NULL;
    }
    if (n <= 0 || n > (Py_ssize_t)INT_MAX)
        return gjp_set_error(module, GJP_ERR_NPTS, NULL);

    x = (double *)PyMem_Malloc((size_t)n * sizeof(double));
    w = (double *)PyMem_Malloc((size_t)n * sizeof(double));
    if (!x || !w) {
        PyMem_Free(x);
        PyMem_Free(w);
        return PyErr_NoMemory();
    }

    for (i = 0; i < n; i++) {
        PyObject *xi = PySequence_GetItem(nodes_obj, i);
        PyObject *wi = PySequence_GetItem(weights_obj, i);
        if (!xi || !wi) {
            Py_XDECREF(xi);
            Py_XDECREF(wi);
            PyMem_Free(x);
            PyMem_Free(w);
            return NULL;
        }
        x[i] = PyFloat_AsDouble(xi);
        w[i] = PyFloat_AsDouble(wi);
        Py_DECREF(xi);
        Py_DECREF(wi);
        if (PyErr_Occurred()) {
            PyMem_Free(x);
            PyMem_Free(w);
            return NULL;
        }
    }

    Py_BEGIN_ALLOW_THREADS
    status = gauss_jacobi_scale_c((int)n, alpha, beta, a, b, x, w);
    Py_END_ALLOW_THREADS

    if (status != GJP_OK) {
        PyMem_Free(x);
        PyMem_Free(w);
        return gjp_set_error(module, status, NULL);
    }

    PyObject *x_list = PyList_New(n);
    PyObject *w_list = PyList_New(n);
    if (!x_list || !w_list) {
        Py_XDECREF(x_list);
        Py_XDECREF(w_list);
        PyMem_Free(x);
        PyMem_Free(w);
        return NULL;
    }
    for (i = 0; i < n; i++) {
        PyObject *xv = PyFloat_FromDouble(x[i]);
        PyObject *wv = PyFloat_FromDouble(w[i]);
        if (!xv || !wv) {
            Py_XDECREF(xv);
            Py_XDECREF(wv);
            Py_DECREF(x_list);
            Py_DECREF(w_list);
            PyMem_Free(x);
            PyMem_Free(w);
            return NULL;
        }
        if (PyList_SetItem(x_list, i, xv) < 0 || PyList_SetItem(w_list, i, wv) < 0) {
            Py_DECREF(x_list);
            Py_DECREF(w_list);
            PyMem_Free(x);
            PyMem_Free(w);
            return NULL;
        }
    }
    PyMem_Free(x);
    PyMem_Free(w);
    PyObject *pair = PyTuple_Pack(2, x_list, w_list);
    Py_DECREF(x_list);
    Py_DECREF(w_list);
    return pair;
}

static PyObject *
gjp_status_string_py(PyObject *module, PyObject *args)
{
    int status;
    (void)module;
    if (!PyArg_ParseTuple(args, "i:status_string", &status))
        return NULL;
    const char *msg = gjp_status_string(status);
    if (!msg)
        msg = "unknown status";
    return PyUnicode_FromString(msg);
}

static PyMethodDef gjp_methods[] = {
    {"rule", (PyCFunction)gjp_rule, METH_VARARGS | METH_KEYWORDS,
     "rule(npts, alpha, beta, method=None) -> (nodes, weights)\n"
     "Call gauss_jacobi_rule_c. method None/\"auto\" selects policy."},
    {"scale", (PyCFunction)gjp_scale, METH_VARARGS | METH_KEYWORDS,
     "scale(nodes, weights, alpha, beta, a, b) -> (nodes, weights)\n"
     "Map [-1, 1] nodes/weights onto [a, b] via gauss_jacobi_scale_c."},
    {"status_string", gjp_status_string_py, METH_VARARGS,
     "status_string(status: int) -> str"},
    {NULL, NULL, 0, NULL}
};

/* ---- multi-phase exec ---- */
static int
gjp_mod_exec(PyObject *module)
{
    gjp_module_state *st = gjp_state(module);
    if (!st)
        return -1;

    st->error = PyErr_NewException("gauss_jacobi_quad.GaussJacobiError", PyExc_RuntimeError, NULL);
    if (!st->error)
        return -1;
    Py_INCREF(st->error);
    if (PyModule_AddObject(module, "GaussJacobiError", st->error) < 0) {
        Py_DECREF(st->error);
        st->error = NULL;
        return -1;
    }

    if (PyModule_AddIntConstant(module, "GJP_OK", GJP_OK) < 0)
        return -1;
    if (PyModule_AddIntConstant(module, "GJP_ERR_NPTS", GJP_ERR_NPTS) < 0)
        return -1;
    if (PyModule_AddIntConstant(module, "GJP_ERR_ALPHA", GJP_ERR_ALPHA) < 0)
        return -1;
    if (PyModule_AddIntConstant(module, "GJP_ERR_BETA", GJP_ERR_BETA) < 0)
        return -1;
    if (PyModule_AddIntConstant(module, "GJP_ERR_METHOD", GJP_ERR_METHOD) < 0)
        return -1;
    if (PyModule_AddIntConstant(module, "GJP_ERR_BOGAERT_AB", GJP_ERR_BOGAERT_AB) < 0)
        return -1;
    if (PyModule_AddIntConstant(module, "GJP_ERR_BOGAERT_N", GJP_ERR_BOGAERT_N) < 0)
        return -1;
    if (PyModule_AddIntConstant(module, "GJP_ERR_BOUNDS", GJP_ERR_BOUNDS) < 0)
        return -1;

    if (PyModule_AddStringConstant(module, "__version__", "0.2.6") < 0)
        return -1;

    return 0;
}

static PyModuleDef_Slot gjp_slots[] = {
    {Py_mod_exec, (void *)gjp_mod_exec},
#if defined(Py_mod_multiple_interpreters) && defined(Py_MOD_PER_INTERPRETER_GIL_SUPPORTED)
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#endif
#if defined(Py_mod_gil) && defined(Py_MOD_GIL_NOT_USED)
    /* Free-threaded full-API builds only (not available under Limited API). */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static struct PyModuleDef gjp_module = {
    PyModuleDef_HEAD_INIT,
    "gauss_jacobi_quad._core",
    "CPython extension: GaussJacobiQuad via ISO_C_BINDING (Stable ABI + free-threaded builds)",
    sizeof(gjp_module_state),
    gjp_methods,
    gjp_slots,
    gjp_traverse,
    gjp_clear,
    gjp_free,
};

PyMODINIT_FUNC
PyInit__core(void)
{
    return PyModuleDef_Init(&gjp_module);
}
