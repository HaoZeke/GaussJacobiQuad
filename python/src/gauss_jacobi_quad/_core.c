/*
 * CPython extension for GaussJacobiQuad.
 *
 * Calls ISO_C_BINDING C ABI: gauss_jacobi_rule_c / gjp_status_string
 * (same objects as libgjp_cinterp, linked into this module).
 *
 * - Multi-phase module init (PEP 489): PyModuleDef_Init + Py_mod_exec
 * - Module state (exception type is not a process-global static)
 * - Heap / dynamic exception type via PyErr_NewException
 * - Free-threaded CPython 3.13+: Py_mod_gil = Py_MOD_GIL_NOT_USED
 * - Multi-interpreter: Py_MOD_PER_INTERPRETER_GIL_SUPPORTED (3.12+)
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "GaussJacobiQuadCInterp.h"

#include <limits.h>
#include <string.h>

/* ---- module state ---- */
typedef struct {
    PyObject *error; /* heap exception type */
} gjp_module_state;

static inline gjp_module_state *
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
 * or with numpy if available we still return lists for minimal deps in C;
 * __init__.py can wrap to ndarray.
 */
static PyObject *
gjp_rule(PyObject *module, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"npts", "alpha", "beta", "method", NULL};
    Py_ssize_t npts_ss = 0;
    double alpha = 0.0, beta = 0.0;
    const char *method = NULL;
    PyObject *method_obj = Py_None;

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
        method = PyUnicode_AsUTF8(method_obj);
        if (!method)
            return NULL;
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
        return PyErr_NoMemory();
    }

    int status;
    Py_BEGIN_ALLOW_THREADS
    status = gauss_jacobi_rule_c(npts, alpha, beta, x, w, method ? method : "");
    Py_END_ALLOW_THREADS

    if (status != GJP_OK) {
        PyMem_Free(x);
        PyMem_Free(w);
        return gjp_set_error(module, status, method && method[0] ? method : "auto");
    }

    PyObject *x_list = PyList_New(npts);
    PyObject *w_list = PyList_New(npts);
    if (!x_list || !w_list) {
        Py_XDECREF(x_list);
        Py_XDECREF(w_list);
        PyMem_Free(x);
        PyMem_Free(w);
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
            return NULL;
        }
        PyList_SET_ITEM(x_list, i, xv); /* steals */
        PyList_SET_ITEM(w_list, i, wv);
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

    /* Heap exception type (dynamic), not a static PyTypeObject */
    st->error = PyErr_NewException("gauss_jacobi_quad.GaussJacobiError", PyExc_RuntimeError, NULL);
    if (!st->error)
        return -1;
    Py_INCREF(st->error);
    if (PyModule_AddObject(module, "GaussJacobiError", st->error) < 0) {
        Py_DECREF(st->error);
        st->error = NULL;
        return -1;
    }

    /* Status constants */
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

    if (PyModule_AddStringConstant(module, "__version__", "0.2.0") < 0)
        return -1;

    return 0;
}

static PyModuleDef_Slot gjp_slots[] = {
    {Py_mod_exec, (void *)gjp_mod_exec},
#if defined(Py_mod_multiple_interpreters) && defined(Py_MOD_PER_INTERPRETER_GIL_SUPPORTED)
    /* Isolated subinterpreters: all mutable state is in module state. */
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#endif
#if defined(Py_mod_gil) && defined(Py_MOD_GIL_NOT_USED)
    /* Free-threaded CPython: declare no GIL requirement. Fortran work runs under
     * Py_BEGIN_ALLOW_THREADS; exception type lives in per-module state. */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static struct PyModuleDef gjp_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "gauss_jacobi_quad._core",
    .m_doc = "CPython extension: GaussJacobiQuad via ISO_C_BINDING C ABI",
    .m_size = sizeof(gjp_module_state),
    .m_methods = gjp_methods,
    .m_slots = gjp_slots,
    .m_traverse = gjp_traverse,
    .m_clear = gjp_clear,
    .m_free = gjp_free,
};

PyMODINIT_FUNC
PyInit__core(void)
{
    return PyModuleDef_Init(&gjp_module);
}
