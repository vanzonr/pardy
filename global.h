///
/// @file global.h
///
/// @brief A few global setting for all of pardy.
///
/// @author Ramses van Zon
/// @date 2026
///
#ifndef GLOBALH
#define GLOBALH

#include <array>

// common compile-type constants
constexpr const int DIM = 3;
constexpr const int global_root = 0;

// global_rank and global_size should be set, once, at the start of main.
extern int global_size; 
extern int global_rank;

// common types
template <class T>
using vec = std::array<T,DIM>;
using tiny_int = signed char;

#endif
