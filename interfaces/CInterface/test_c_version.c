// BEGIN_HEADER
// -----------------------------------------------------------------------------
// Gauss-Jacobi Quadrature Implementation
// Authors: Rohit Goswami <rgoswami[at]ieee.org>
// Source: GaussJacobiQuad Library
// License: MIT
// GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
// -----------------------------------------------------------------------------
// END_HEADER
/* Generated version macros must match the meson project version. */
#include <stdio.h>
#include <string.h>

#include "gjp_version.h"

#ifndef EXPECT_MAJOR
#error "EXPECT_MAJOR must be set from the project version"
#endif
#ifndef EXPECT_MINOR
#error "EXPECT_MINOR must be set from the project version"
#endif
#ifndef EXPECT_PATCH
#error "EXPECT_PATCH must be set from the project version"
#endif
#ifndef EXPECT_STRING
#error "EXPECT_STRING must be set from the project version"
#endif

int main(void) {
    int compact = 10000 * EXPECT_MAJOR + 100 * EXPECT_MINOR + EXPECT_PATCH;

    printf("gjp %d.%d.%d string '%s' compact %d expect %d.%d.%d '%s'\n",
           GJP_VERSION_MAJOR, GJP_VERSION_MINOR, GJP_VERSION_PATCH, GJP_VERSION_STRING, GJP_VERSION,
           EXPECT_MAJOR, EXPECT_MINOR, EXPECT_PATCH, EXPECT_STRING);
    if (GJP_VERSION_MAJOR != EXPECT_MAJOR || GJP_VERSION_MINOR != EXPECT_MINOR ||
        GJP_VERSION_PATCH != EXPECT_PATCH || strcmp(GJP_VERSION_STRING, EXPECT_STRING) != 0 ||
        GJP_VERSION != compact) {
        fprintf(stderr, "version mismatch\n");
        return 1;
    }
    return 0;
}
