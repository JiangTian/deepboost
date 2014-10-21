# SYNOPSIS:
#
#   make - make everything
#   make test - make and run all tests
#   make clean - remove all files generated by make
#   make driver - make the main executable

# LIB_DIR should satisfy the following:
#   LIB_DIR/include/gflags contains Google Commandline Flags include files
#   LIB_DIR/include/glog contains Google Logging include files
#   LIB_DIR/include/gtest contains Google Test include files
#   LIB_DIR/src contains Google Test source files
#   LIB_DIR/lib contains libgflags* and libglog* library files.
LIB_DIR = /home/usyed/googleopensource

# Where to find user code.
USER_DIR = .

# Flags passed to the preprocessor.
CPPFLAGS += -isystem $(LIB_DIR)/include

# Flags passed to the C++ compiler. Add -O3 for the highest optimization level.
# Add -ggdb for GDB debugging info.
CXXFLAGS += -Wall -Wextra -pthread -std=c++0x

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = tree_test boost_test io_test

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(LIB_DIR)/include/gtest/*.h \
                $(LIB_DIR)/include/gtest/internal/*.h

# House-keeping build targets.

test: $(TESTS)
	./tree_test
	./io_test
	./boost_test
clean :
	rm -f $(TESTS) gtest_main.a driver *.o

# Builds gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(LIB_DIR)/src/*.cc $(LIB_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(LIB_DIR) $(CXXFLAGS) -c \
            $(LIB_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(LIB_DIR) $(CXXFLAGS) -c \
            $(LIB_DIR)/src/gtest_main.cc

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds tests.  A test should link with gtest_main.a.

tree.o : $(USER_DIR)/tree.cc $(USER_DIR)/tree.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/tree.cc

tree_test.o : $(USER_DIR)/tree_test.cc \
                     $(USER_DIR)/tree.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/tree_test.cc

tree_test : tree.o tree_test.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -static -lpthread $^ -o $@ -L$(LIB_DIR)/lib -lgflags -lglog

boost.o : $(USER_DIR)/boost.cc $(USER_DIR)/boost.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/boost.cc

boost_test.o : $(USER_DIR)/boost_test.cc \
                     $(USER_DIR)/boost.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/boost_test.cc

boost_test : tree.o boost.o boost_test.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -static -lpthread $^ -o $@ -L$(LIB_DIR)/lib -lgflags -lglog

io.o : $(USER_DIR)/io.cc $(USER_DIR)/io.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/io.cc

io_test.o : $(USER_DIR)/io_test.cc \
                     $(USER_DIR)/io.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/io_test.cc

io_test : tree.o io.o io_test.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -static -lpthread $^ -o $@ -L$(LIB_DIR)/lib -lgflags -lglog

# Build the main executable

driver.o : $(USER_DIR)/driver.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/driver.cc

driver : tree.o boost.o io.o driver.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -static -lpthread $^ -o $@ -L$(LIB_DIR)/lib -lgflags -lglog