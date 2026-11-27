/* Shared pkg-config consumer. gjp_get_version lives in libGaussJacobiQuad. */
#include "gjp_version.h"

extern void __gjp_version_MOD_gjp_get_version(void);

int main(void) {
    void (*fn)(void) = &__gjp_version_MOD_gjp_get_version;
    return fn == 0 || GJP_VERSION_MAJOR < 0;
}
