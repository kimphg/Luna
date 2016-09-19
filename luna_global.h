#ifndef LUNA_GLOBAL_H
#define LUNA_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(LUNA_LIBRARY)
#  define LUNASHARED_EXPORT Q_DECL_EXPORT
#else
#  define LUNASHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // LUNA_GLOBAL_H
