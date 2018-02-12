#pragma once

#include <cmath>
#include <cassert>
#include <iostream>
#include <memory>
#include <type_traits>

#include "dblas.h"
#include "dsparse_csr.h"
#include "dsparse_csc.h"
#include "dsparse_coo.h"

#include "marlib_traits.hpp"
#include "marlib_dcopy.hpp"
#include "marlib_dnnz.hpp"
#include "marlib_mapapply.hpp"

// #include "cppblas.hpp"
// #include "cppapply.hpp"
//
// #include "poisson.hpp"
//
// #include "mexp_pade.hpp"
// #include "mexp_unif.hpp"
// #include "mexpint_unif.hpp"
// #include "mexpconv_unif.hpp"
//
// #include "gaussinte.hpp"
