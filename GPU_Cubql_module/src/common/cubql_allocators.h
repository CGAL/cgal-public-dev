// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#pragma once

// Ensure CUBQL_CUDA_CALL or any underlying CUDA types are known
#include "cuBQL/builder/cuda/builder_common.h" 

// --------------------------------------------------------------------
// CENTRAL STREAM-SAFE REPLACEMENTS FOR INTERNAL CUBQL ALLOCATORS
// --------------------------------------------------------------------
#ifndef _ALLOC
#define _ALLOC(ptr, count, stream, memResource) \
    CUBQL_CUDA_CALL(MallocAsync((void**)&(ptr), (size_t)(count) * sizeof(*(ptr)), stream))
#endif

#ifndef _FREE
#define _FREE(ptr, stream, memResource) \
    CUBQL_CUDA_CALL(FreeAsync((ptr), stream))
#endif