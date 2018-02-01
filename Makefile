# Makefile

CXX = g++
# CXXFLAGS = --std=c++11 -fPIC -g -Wall -I include -DF77BLAS -DF77LAPACK -DDEBUG
CXXFLAGS = --std=c++11 -fPIC -g -Wall -I include -DF77BLAS -DF77LAPACK
LDFLAGS = -lblas -llapack

TEST_DBLAS = test/test_dblas.o
TEST_CPPBLAS = test/test_dblascpp.o
# TEST_SPBLAS = test/test_spblas.o
# TEST_POISSON = test/test_poisson.o
# TEST_MEXP = test/test_mexp.o
# TEST_PH = test/test_ph.o
# TEST_CTMC = test/test_ctmc.o
# TEST_GAMMA = test/test_gamma.o

# default target

test: test_dblas test_cppblas

# test dblas

# lib: $(OBJS)
# 	$(CXX) $(LDFLAGS) -shared -o libmarlib.so $^

test_dblas: $(TEST_DBLAS) $(OBJS)
	$(CXX) $(LDFLAGS) -o test/$@.out $^
	test/$@.out | tee test/$@.result

test_cppblas: $(TEST_CPPBLAS) $(OBJS)
	$(CXX) $(LDFLAGS) -o test/$@.out $^
	test/$@.out | tee test/$@.result

# test_spblas: $(TEST_SPBLAS) $(OBJS)
# 	$(CXX) $(LDFLAGS) -o test/$@.out $^
# 	test/$@.out | tee test/$@.result

# test_poisson: $(TEST_POISSON) $(OBJS)
# 	$(CXX) $(LDFLAGS) -o test/$@.out $^
# 	test/$@.out | tee test/$@.result
#
# test_mexp: $(TEST_MEXP) $(OBJS)
# 	$(CXX) $(LDFLAGS) -o test/$@.out $^
# 	test/$@.out | tee test/$@.result
#
# test_ph: $(TEST_PH) $(OBJS)
# 	$(CXX) $(LDFLAGS) -o test/$@.out $^
# 	test/$@.out | tee test/$@.result
#
# test_ctmc: $(TEST_CTMC) $(OBJS)
# 	$(CXX) $(LDFLAGS) -o test/$@.out $^
# 	test/$@.out | tee test/$@.result
#
# test_gamma: $(TEST_GAMMA) $(OBJS)
# 	$(CXX) $(LDFLAGS) -o test/$@.out $^
# 	test/$@.out | tee test/$@.result

# suffix

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

# clean

clean:
	rm -f $(OBJS) $(TEST_DBLAS) \
	$(TEST_CPPBLAS) \
	$(TEST_SPBLAS) \
	$(TEST_POISSON) \
	$(TEST_MEXP) \
	$(TEST_PH) \
	$(TEST_CTMC) \
	$(TEST_GAMMA)
