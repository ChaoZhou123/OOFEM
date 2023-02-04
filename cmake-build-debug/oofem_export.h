
#ifndef OOFEM_EXPORT_H
#define OOFEM_EXPORT_H

#ifdef OOFEM_STATIC_DEFINE
#  define OOFEM_EXPORT
#  define OOFEM_NO_EXPORT
#else
#  ifndef OOFEM_EXPORT
#    ifdef liboofem_EXPORTS
        /* We are building this library */
#      define OOFEM_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define OOFEM_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef OOFEM_NO_EXPORT
#    define OOFEM_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef OOFEM_DEPRECATED
#  define OOFEM_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef OOFEM_DEPRECATED_EXPORT
#  define OOFEM_DEPRECATED_EXPORT OOFEM_EXPORT OOFEM_DEPRECATED
#endif

#ifndef OOFEM_DEPRECATED_NO_EXPORT
#  define OOFEM_DEPRECATED_NO_EXPORT OOFEM_NO_EXPORT OOFEM_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef OOFEM_NO_DEPRECATED
#    define OOFEM_NO_DEPRECATED
#  endif
#endif

#endif /* OOFEM_EXPORT_H */
