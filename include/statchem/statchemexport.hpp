#ifndef STATCHEM_EXPORT_H_
#define STATCHEM_EXPORT_H_

#ifdef _MSC_VER
// We don't want to hear about how sprintf is "unsafe".
#pragma warning(disable : 4996)
// Keep MS VC++ quiet about lack of dll export of private members.
#pragma warning(disable : 4251)
#if defined(STATCHEM_SHARED_LIBRARY)
#define STATCHEM_EXPORT __declspec(dllexport)
#elif defined(STATCHEM_STATIC_LIBRARY) || defined(STATCHEM_USE_STATIC_LIBRARIES)
#define STATCHEM_EXPORT
#else
#define STATCHEM_EXPORT __declspec(dllimport)
#endif
#else
#define STATCHEM_EXPORT
#endif

#endif
