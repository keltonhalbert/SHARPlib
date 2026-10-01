// Consumer check for the NO_QC build option. Header templates and inline
// functions such as layer_minmax compile in the consumer's translation
// unit, so a consumer must see the library's NO_QC setting. This target only
// links SHARPlib and defines no NO_QC of its own. CMake passes the library's
// setting as SHARPLIB_EXPECT_NO_QC (0 or 1), and this file fails to compile
// if the two differ.
#include <SHARPlib/layer.h>

#if defined(NO_QC) != SHARPLIB_EXPECT_NO_QC
#error "NO_QC did not propagate from the SHARPlib target to this consumer"
#endif

int main() { return 0; }
