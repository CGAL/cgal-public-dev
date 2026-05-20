# --- PATHS ---
CUBQL_ROOT ?= /home/yury/Projects/CPPProj/Cubql
CUDA_PATH  ?= /usr/local/cuda
LIBS_DIR    = $(CUBQL_ROOT)/build

NVCC        = $(CUDA_PATH)/bin/nvcc

# --- FLAGS ---
CC_FLAGS    = -O3 -std=c++17 --expt-relaxed-constexpr -rdc=true --extended-lambda -DCUBQL_HAVE_CUDA=1
GENCODE     = -gencode arch=compute_120,code=sm_120

# --- INCLUDES ---
INCLUDES    = -I$(CUBQL_ROOT) \
              -I$(CUDA_PATH)/include \
              -I$(CUBQL_ROOT)/samples/common \
              -I.

# --- LIBRARIES ---
LIBS        = -L$(LIBS_DIR) -Xlinker -rpath,$(LIBS_DIR) \
              -lcuBQL_cuda_float2 -lcuBQL_cuda_float3 -lcuBQL_samples_common \
              -lcudart -lcudadevrt

# --- TARGETS ---
all: meshIntersectionRegistrationOld

# Reads your actual file, outputs the old binary target name
meshIntersectionRegistrationOld: meshIntersectionRegistration.cu
	$(NVCC) $(CC_FLAGS) $(GENCODE) $(INCLUDES) meshIntersectionRegistration.cu $(LIBS) -o meshIntersectionRegistrationOld

clean:
	rm -f meshIntersectionRegistrationOld