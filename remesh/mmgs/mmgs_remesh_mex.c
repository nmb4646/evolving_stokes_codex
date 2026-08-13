#include "mex.h"
#include "matrix.h"
#include "mmg/mmgs/libmmgs.h"

#include <math.h>
#include <string.h>

static double option_double(const mxArray *options, const char *name, double fallback) {
    const mxArray *field;
    if (options == NULL || !mxIsStruct(options)) {
        return fallback;
    }
    field = mxGetField(options, 0, name);
    if (field == NULL || mxIsEmpty(field)) {
        return fallback;
    }
    if (!mxIsNumeric(field) || mxGetNumberOfElements(field) != 1) {
        mexErrMsgIdAndTxt("mmgs:InvalidOption", "Option '%s' must be a numeric scalar.", name);
    }
    return mxGetScalar(field);
}

static int option_int(const mxArray *options, const char *name, int fallback) {
    return (int)llround(option_double(options, name, (double)fallback));
}

static void set_status(mxArray **output, int return_code, mwSize np_before,
                       mwSize nt_before, MMG5_int np_after, MMG5_int nt_after) {
    const char *fields[] = {"return_code", "success", "low_failure", "strong_failure",
                            "n_vertices_before", "n_faces_before",
                            "n_vertices_after", "n_faces_after"};
    mxArray *status = mxCreateStructMatrix(1, 1, 8, fields);
    mxSetField(status, 0, "return_code", mxCreateDoubleScalar((double)return_code));
    mxSetField(status, 0, "success", mxCreateLogicalScalar(return_code == MMG5_SUCCESS));
    mxSetField(status, 0, "low_failure", mxCreateLogicalScalar(return_code == MMG5_LOWFAILURE));
    mxSetField(status, 0, "strong_failure", mxCreateLogicalScalar(return_code == MMG5_STRONGFAILURE));
    mxSetField(status, 0, "n_vertices_before", mxCreateDoubleScalar((double)np_before));
    mxSetField(status, 0, "n_faces_before", mxCreateDoubleScalar((double)nt_before));
    mxSetField(status, 0, "n_vertices_after", mxCreateDoubleScalar((double)np_after));
    mxSetField(status, 0, "n_faces_after", mxCreateDoubleScalar((double)nt_after));
    *output = status;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *P_in, *M_in, *h_in, *options = NULL;
    const double *P, *M, *h;
    mwSize np, nt, i;
    MMG5_pMesh mesh = NULL;
    MMG5_pSol metric = NULL;
    MMG5_int np_out = 0, nt_out = 0, na_out = 0;
    int return_code = MMG5_STRONGFAILURE;
    double hmin, hmax, hausd, hgrad;
    int verbose, angle_detection;

    if (nrhs < 3 || nrhs > 4) {
        mexErrMsgIdAndTxt("mmgs:InputCount", "Use mmgs_remesh_mex(P,M,h[,options]).");
    }
    if (nlhs < 2 || nlhs > 3) {
        mexErrMsgIdAndTxt("mmgs:OutputCount", "Request [P2,M2] or [P2,M2,status].");
    }
    P_in = prhs[0];
    M_in = prhs[1];
    h_in = prhs[2];
    if (nrhs == 4) {
        options = prhs[3];
        if (!mxIsStruct(options) || mxGetNumberOfElements(options) != 1) {
            mexErrMsgIdAndTxt("mmgs:InvalidOptions", "options must be a scalar struct.");
        }
    }
    if (!mxIsDouble(P_in) || mxIsComplex(P_in) || mxGetN(P_in) != 3 || mxIsEmpty(P_in)) {
        mexErrMsgIdAndTxt("mmgs:InvalidVertices", "P must be a nonempty real double N-by-3 array.");
    }
    if (!mxIsDouble(M_in) || mxIsComplex(M_in) || mxGetN(M_in) != 3 || mxIsEmpty(M_in)) {
        mexErrMsgIdAndTxt("mmgs:InvalidFaces", "M must be a nonempty real double F-by-3 array.");
    }
    if (!mxIsDouble(h_in) || mxIsComplex(h_in) || mxGetN(h_in) != 1) {
        mexErrMsgIdAndTxt("mmgs:InvalidMetric", "h must be a real double N-by-1 array.");
    }

    np = mxGetM(P_in);
    nt = mxGetM(M_in);
    if (mxGetM(h_in) != np) {
        mexErrMsgIdAndTxt("mmgs:InvalidMetric", "h must have one value per vertex.");
    }
    if (np > (mwSize)MMG5_INTMAX || nt > (mwSize)MMG5_INTMAX) {
        mexErrMsgIdAndTxt("mmgs:MeshTooLarge", "Mesh exceeds MMG integer limits.");
    }
    P = mxGetDoubles(P_in);
    M = mxGetDoubles(M_in);
    h = mxGetDoubles(h_in);
    for (i = 0; i < 3 * np; ++i) {
        if (!isfinite(P[i])) {
            mexErrMsgIdAndTxt("mmgs:InvalidVertices", "P contains NaN or Inf.");
        }
    }
    hmin = h[0];
    hmax = h[0];
    for (i = 0; i < np; ++i) {
        if (!isfinite(h[i]) || h[i] <= 0.0) {
            mexErrMsgIdAndTxt("mmgs:InvalidMetric", "All target sizes must be finite and positive.");
        }
        if (h[i] < hmin) hmin = h[i];
        if (h[i] > hmax) hmax = h[i];
    }
    for (i = 0; i < 3 * nt; ++i) {
        double value = M[i];
        if (!isfinite(value) || value != floor(value) || value < 1.0 || value > (double)np) {
            mexErrMsgIdAndTxt("mmgs:InvalidFaces", "M contains an invalid vertex index.");
        }
    }

    hmin = option_double(options, "hmin", hmin);
    hmax = option_double(options, "hmax", hmax);
    hausd = option_double(options, "hausd", 0.05 * hmin);
    hgrad = option_double(options, "hgrad", 1.3);
    verbose = option_int(options, "verbose", -1);
    angle_detection = option_int(options, "angle_detection", 0);
    if (!(hmin > 0.0 && hmax >= hmin && hausd > 0.0 && hgrad >= 1.0)) {
        mexErrMsgIdAndTxt("mmgs:InvalidOption", "Require hmin>0, hmax>=hmin, hausd>0, and hgrad>=1.");
    }

    if (MMGS_Init_mesh(MMG5_ARG_start,
                       MMG5_ARG_ppMesh, &mesh,
                       MMG5_ARG_ppMet, &metric,
                       MMG5_ARG_end) != 1) {
        mexErrMsgIdAndTxt("mmgs:InitializationFailed", "MMGS_Init_mesh failed.");
    }
    if (MMGS_Set_meshSize(mesh, (MMG5_int)np, (MMG5_int)nt, 0) != 1) {
        MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh,
                      MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
        mexErrMsgIdAndTxt("mmgs:MeshAllocationFailed", "MMGS_Set_meshSize failed.");
    }
    for (i = 0; i < np; ++i) {
        if (MMGS_Set_vertex(mesh, P[i], P[i + np], P[i + 2 * np], 0, (MMG5_int)(i + 1)) != 1) {
            MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh,
                          MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
            mexErrMsgIdAndTxt("mmgs:SetVertexFailed", "Failed to insert vertex %llu.", (unsigned long long)(i + 1));
        }
    }
    for (i = 0; i < nt; ++i) {
        if (MMGS_Set_triangle(mesh, (MMG5_int)M[i], (MMG5_int)M[i + nt],
                              (MMG5_int)M[i + 2 * nt], 0, (MMG5_int)(i + 1)) != 1) {
            MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh,
                          MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
            mexErrMsgIdAndTxt("mmgs:SetTriangleFailed", "Failed to insert triangle %llu.", (unsigned long long)(i + 1));
        }
    }
    if (MMGS_Set_solSize(mesh, metric, MMG5_Vertex, (MMG5_int)np, MMG5_Scalar) != 1 ||
        MMGS_Set_scalarSols(metric, (double *)h) != 1) {
        MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh,
                      MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
        mexErrMsgIdAndTxt("mmgs:SetMetricFailed", "Failed to set the scalar size field.");
    }
    if (MMGS_Chk_meshData(mesh, metric) != 1 ||
        MMGS_Set_iparameter(mesh, metric, MMGS_IPARAM_verbose, verbose) != 1 ||
        MMGS_Set_iparameter(mesh, metric, MMGS_IPARAM_angle, angle_detection) != 1 ||
        MMGS_Set_dparameter(mesh, metric, MMGS_DPARAM_hmin, hmin) != 1 ||
        MMGS_Set_dparameter(mesh, metric, MMGS_DPARAM_hmax, hmax) != 1 ||
        MMGS_Set_dparameter(mesh, metric, MMGS_DPARAM_hausd, hausd) != 1 ||
        MMGS_Set_dparameter(mesh, metric, MMGS_DPARAM_hgrad, hgrad) != 1) {
        MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh,
                      MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
        mexErrMsgIdAndTxt("mmgs:ConfigurationFailed", "MMGS rejected mesh data or options.");
    }

    return_code = MMGS_mmgslib(mesh, metric);
    if (return_code == MMG5_STRONGFAILURE ||
        MMGS_Get_meshSize(mesh, &np_out, &nt_out, &na_out) != 1) {
        MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh,
                      MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
        mexErrMsgIdAndTxt("mmgs:RemeshingFailed", "MMGS reported a strong failure.");
    }

    plhs[0] = mxCreateDoubleMatrix((mwSize)np_out, 3, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)nt_out, 3, mxREAL);
    {
        double *P_out = mxGetDoubles(plhs[0]);
        double *M_out = mxGetDoubles(plhs[1]);
        MMG5_int ref, v0, v1, v2;
        int corner, required;
        for (i = 0; i < (mwSize)np_out; ++i) {
            if (MMGS_Get_vertex(mesh, &P_out[i], &P_out[i + np_out], &P_out[i + 2 * np_out],
                                &ref, &corner, &required) != 1) {
                MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh,
                              MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
                mexErrMsgIdAndTxt("mmgs:GetVertexFailed", "Failed to extract an output vertex.");
            }
        }
        for (i = 0; i < (mwSize)nt_out; ++i) {
            if (MMGS_Get_triangle(mesh, &v0, &v1, &v2, &ref, &required) != 1) {
                MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh,
                              MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
                mexErrMsgIdAndTxt("mmgs:GetTriangleFailed", "Failed to extract an output triangle.");
            }
            M_out[i] = (double)v0;
            M_out[i + nt_out] = (double)v1;
            M_out[i + 2 * nt_out] = (double)v2;
        }
    }
    if (nlhs == 3) {
        set_status(&plhs[2], return_code, np, nt, np_out, nt_out);
    }
    MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh,
                  MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
}
