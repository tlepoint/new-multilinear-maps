# makefile version 20150213
# Tancrede Lepoint
# Public domain.

# To adapt to your settings
INCLUDE_DIRS 		= -I/opt/local/include
LIB_DIRS 			= -L/opt/local/lib
CXX 				= g++
PARALLEL 			= no

################################################################################################
# Targets
GENERATE_TARGET 	= generate_pp
KEY_EXCHANGE_TARGET = key_exchange

# flags
CXXFLAGS			= -Wall -Ofast -g -std=c++11 -funroll-loops $(INCLUDE_DIRS)
LDLIBS 				= $(LIB_DIRS) -lgmp -lgmpxx -lmpfr 
LDLIBS_FPLLL 		= $(LIB_DIRS) -lfplll

ifeq ($(PARALLEL),yes)
  CXXFLAGS+= -fopenmp
endif

# Files 
GENERATE_SRCS 		= generate_pp.cpp $(wildcard mmap/prng/*.cpp)
GENERATE_OBJS 		= $(GENERATE_SRCS:%.cpp=%.o)
KEY_EXCHANGE_SRCS 	= key_exchange.cpp $(wildcard mmap/prng/*.cpp)
KEY_EXCHANGE_OBJS 	= $(KEY_EXCHANGE_SRCS:%.cpp=%.o)
ASMS 				= $(wildcard mmap/prng/*.s)
HDRS 				= $(wildcard mmap/*.hpp) $(wildcard mmap/prng/*.h)

# All & targets
all: $(GENERATE_TARGET) $(KEY_EXCHANGE_TARGET)

$(GENERATE_TARGET): $(GENERATE_OBJS) $(ASMS) $(HDRS)
		$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -o $@ $(GENERATE_OBJS) $(ASMS) $(LDLIBS) $(LDLIBS_FPLLL)

$(KEY_EXCHANGE_TARGET): $(KEY_EXCHANGE_OBJS) $(ASMS) $(HDRS)
		$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -o $@ $(KEY_EXCHANGE_OBJS) $(ASMS) $(LDLIBS)

# Dependencies of objs
$(GENERATE_OBJS): $(HDRS) makefile
$(KEY_EXCHANGE_OBJS): $(HDRS) makefile

# cleaning
clean:
	$(RM) -r $(GENERATE_TARGET) $(GENERATE_TARGET).dSYM $(GENERATE_OBJS) $(KEY_EXCHANGE_TARGET) $(KEY_EXCHANGE_TARGET).dSYM $(KEY_EXCHANGE_OBJS)

