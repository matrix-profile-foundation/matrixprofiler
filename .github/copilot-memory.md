# Copilot Memory for matrixprofiler Project

**Last Updated:** 2024-12-16

## Recent Work Completed

### MPX Code Standardization and Enhancement (2024-12-16)

#### src/mpx.cpp - Complete Optimization
**All 7 major functions optimized:**

1. **mpxis_rcpp** (lines 17-189) - Streaming right-side with threshold
   - Removed const from: data, ddf_t, ddg_t, diag_end, compute_order, num_progress, stop, ww, mp_new_size, diff
   - Changed to ::Rf_error
   - IntegerVector loop uses const & pattern

2. **mpxi_rcpp** (lines 194-376) - Streaming incremental
   - Removed const from: profile_len, diag_end, num_progress, stop, mp_new_size, diff
   - Changed to ::Rf_error
   - Maintained: NumericVector const ww for fixed reference

3. **mpxiright_rcpp** (lines 379-520) - Right-side recompute with restart
   - Removed const from: exclusion_zone, profile_len, diag_end, num_progress
   - Changed to ::Rf_error
   - All IntegerVector loops use const & pattern

4. **mpxileft_rcpp** (lines 522-667) - Left-side recompute
   - Removed const from: exclusion_zone, profile_len, num_progress
   - Changed to ::Rf_error
   - All IntegerVector loops use const & pattern

5. **mpx_rcpp_new** (lines 675-907) - Experimental 3-profile version
   - Removed const from: exclusion_zone, profile_len, diag_start, num_progress
   - Changed to ::Rf_error
   - **Features:** Calculates mmp + rmmp + lmmp, wild sigma check (sig > 60), mp_time_constraint parameter
   - **Status:** NOT exported (no [[Rcpp::export]]) - under active development
   - **FRANZ TEST block:** Experimental averaged profiles code (commented out, not validated)
   - **Documentation added:** Comprehensive comparison comments explaining differences vs mpx_rcpp

6. **mpx_rcpp** (lines 915-1032) - Stable exported version
   - Removed const from: exclusion_zone, profile_len, num_progress
   - Changed to ::Rf_error
   - **NEW FEATURE:** wild_sigma_threshold parameter (default R_PosInf, disabled)
   - Skip pairs with abnormally high sigma: `if (sig > threshold) continue;`
   - Backward compatible (opt-in feature)

7. **mpxab_rcpp** (lines 1034-1210) - AB-join variant
   - Removed const from: exclusion_zone, a_len, b_len, profile_len_a, profile_len_b, num_progress
   - Changed to ::Rf_error
   - Both IntegerVector loops (AB and BA) use const & pattern

#### All 8 IntegerVector Iteration Loops Fixed
**Pattern applied:** `for (int32_t const &diag : compute_order)`
**Reason:** Avoids Rcpp proxy→int32_t conversion overhead
**Locations:** Lines 87, 309, 457, 604, 749, 974, 1104, 1152

#### src/mpx.h - Header Updated
- Added `wild_sigma_threshold = R_PosInf` parameter to mpx_rcpp declaration (line 20)
- Headers synchronized with implementation

### Code Standardization in windowfunc files (2024-12-16)

#### R/windowfunc.R
**Completed changes:**
1. ✅ Standardized all `window_size < 2` comparisons to `window_size < 2L` (integer literal)
2. ✅ Added `# nolint` comments to Rcpp function calls in:
   - `mov_std()` line 292
   - `movmean_std()` line 327
   - `muinvn()` lines 387 and 389
3. ✅ Standardized array indices to use integer literals:
   - `data[1:window_size]` → `data[1L:window_size]`
   - Applied in `mov_std()` and `movmean_std()`
4. ✅ Explicit numeric literals for clarity:
   - `data^2` → `data^2.0`
   - `data_mean^2` → `data_mean^2.0`
   - `< 0` → `< 0.0`
   - `sqrt(1 / ...)` → `sqrt(1.0 / ...)`
   - Applied in `mov_std()` else branch and `movmean_std()` else branch

**Pattern:** All validation checks now use explicit integer literals (`2L`) for consistency with parameter type conversion via `as.integer()`.

#### src/windowfunc.cpp
**Completed changes:**
1. ✅ Removed `const` qualifiers from local `NumericVector` variables in:
   - `movstd_rcpp()` - removed from mu, data2_sum, data_var, data_sd
   - `movmean_std_rcpp()` - removed from data_sum, data_mean, data2, data2_sum, data_var, data_sd, data_sig
   - `movvar_rcpp()` - removed from mu, data2_sum, data_var
   - `muinvn_rcpp()` - removed from mu, data2_sum, sig
2. ✅ Kept `const` on:
   - Function parameters (e.g., `const NumericVector data`)
   - Primitive types where meaningful (e.g., `uint32_t const data_size`)
3. ✅ Fixed `const` correctness in `muinvn_rcpp_parallel()`:
   - Removed `const` from `sig` variable (line 501) - it's modified by Worker
   - Kept `const` on `mu` and `data2_sum` (read-only in Worker)
4. ✅ Changed loop variable types in `MuinWorker::operator()`:
   - `uint32_t i` → matches std::size_t parameter semantics (line 494)

**Pattern:** Local `NumericVector` temporaries no longer use `const` for consistency. Only output variables modified by parallel Workers are non-const; input variables remain const.

#### src/mpx_parallel.cpp
**Completed changes:**
1. ✅ Changed loop variable types in Worker::operator() methods to `std::size_t`:
   - `MatrixProfileP::operator()` - line 70: `for (std::size_t dd = begin; dd < end; dd++)`
   - `MatrixProfilePAB::operator()` - line 241: `for (std::size_t dd = begin; dd < end; dd++)`
   - `MatrixProfilePAB::operator()` - line 265: `for (std::size_t dd = begin; dd < end; dd++)`
2. ✅ Reasoning: `std::size_t` matches the type of `begin` and `end` parameters in Worker operators
3. ✅ Type safety: Avoids implicit conversions between uint32_t and std::size_t

**Pattern:** All Worker::operator() methods use `std::size_t` for loop variables to match RcppParallel signature.

### Key Technical Patterns Discovered

#### Rcpp Proxy Efficiency
**IntegerVector iteration:** `for (int32_t const &diag : vector)` is MORE efficient than `for (int32_t diag : vector)`
- Reason: Rcpp containers use proxy objects; `const &` avoids proxy→value conversion
- Applied to all 8 IntegerVector loops in mpx.cpp
- Performance benefit without code complexity increase

#### Const Correctness Rules
1. **Remove const from:** Variables that are computed or reassigned (profile_len, diag_end, num_progress, stop, etc.)
2. **Keep const for:** Fixed reference vectors that are never modified (NumericVector const ww)
3. **Documents intent:** const on fixed references shows immutability clearly

#### Wild Sigma Check Feature
**Purpose:** Prevents misleading correlations with abnormally high standard deviation
**Implementation:** `if ((sig[offset] > threshold) || (sig[off_diag] > threshold)) continue;`
**Threshold:** 60.0 used in mpx_rcpp_new (experimental), configurable in mpx_rcpp (default disabled)
**Use cases:** Sensor data, financial data, IoT data prone to outliers

## Code Style Guidelines Established

### R Code
- Use `2L` for integer comparisons with `window_size`
- Use `1L` for array start indices
- Use `.0` suffix for numeric literals in mathematical operations for clarity
- Add `# nolint` to suppress false positive linter warnings on Rcpp calls

### C++ Code
- Remove `const` from local `NumericVector` temporary variables that are reassigned
- Keep `const` on fixed reference vectors (e.g., `NumericVector const ww` for window reference)
- Keep `const` on function parameters
- Keep `const` on primitive types (uint32_t, double) when they shouldn't change
- Variables modified by RcppParallel Workers must not be const
- Use `std::size_t` for loop variables in Worker::operator() methods
- Use `int32_t const &` for IntegerVector range-based for loops (avoids proxy conversion)
- Prefer `::Rf_error()` over `Rcpp::stop()` for C API efficiency

## Known Issues/Next Steps
- mpx.cpp optimization complete - ready for unit testing
- wild_sigma_threshold feature needs validation with test datasets
- mpx_rcpp_new FRANZ TEST block needs validation before export
- Consider documenting wild sigma threshold scientific basis (UCR papers)

## Project Context
- Working on matrixprofiler R package
- Focus on code consistency and C++ const correctness
- Comparing optimizations from external R project for potential integration
