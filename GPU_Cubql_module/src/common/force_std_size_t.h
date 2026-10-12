#pragma once
// Force-included (cl /FI) into every CUDA translation unit on MSVC.
//
// nvcc's front end re-emits the system headers into <file>.cudafe1.cpp in its own order. With
// MSVC 14.44 / CUDA 13 the emitted "enum class align_val_t : std::size_t" (vcruntime_new.h) can
// land before "namespace std { using ::size_t; }", giving error C2039 'size_t' is not a member
// of 'std'. Declaring std::size_t first removes the ordering problem.
//
// IMPORTANT: this header must NOT #include anything. The cudafe1.cpp file contains the system
// headers' contents inline (without include guards), so a real header included here would be
// defined twice (errors C2011/C2953 redefinition).
namespace std { typedef decltype(sizeof(0)) size_t; }
