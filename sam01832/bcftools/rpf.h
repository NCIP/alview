
#if 1 // VISUAL C BULLSTUFF ...
// #define isnan(x) ((x) != (x))
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#define fpu_error(x) (isinf(x) || isnan(x))
#endif

