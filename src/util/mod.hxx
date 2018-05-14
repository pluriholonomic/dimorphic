#ifndef __MOD_HXX__
#define __MOD_HXX__

// Place to store our hacky mod functions
namespace Dimers {

inline int mod_close(int r, int modulus) {
    return (r < 0) ? r + modulus : r % modulus;
}

}  // namespace Dimers

#endif  // __MOD_HXX__
