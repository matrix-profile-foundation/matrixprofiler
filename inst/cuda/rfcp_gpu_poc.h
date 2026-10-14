#ifndef MATRIXPROFILER_RFCP_GPU_POC_H
#define MATRIXPROFILER_RFCP_GPU_POC_H

#include <cstdint>

namespace matrixprofiler
{
  namespace gpu
  {

    enum class RfcpGpuStatus : int
    {
      ok = 0,
      invalid_argument = 1,
      cuda_unavailable = 2,
      cuda_error = 3,
    };

    // Computes the RFCP rms_profile for query indices [query_begin, query_end).
    //
    // All input arrays are host-resident; the function manages device memory
    // internally.  Statistics (mean, sig, valid_window) must be the stable,
    // recurrence-safe values produced by the package's muinvn_na helper — the
    // same pre-computation used by the CPU NA-aware path.
    //
    // Outputs:
    //   rms_profile[i]       - RMSC for query (query_begin + i), or NaN if
    //                          the query window is invalid.
    //   paired_frequency[i]  - number of paired (AA, AB) neighbor ranks
    //                          accepted for query (query_begin + i).
    //   plato_index          - 0-based global index of the query with the
    //                          highest rms_profile (set only if plato found).
    //   plato_rms            - rms_profile value at plato_index (NaN if none).
    //   plato_dist_aa        - host buffer of positive_profiles doubles
    //                          containing the raw AA Euclidean distances at the
    //                          plato query (optional; pass nullptr to skip).
    //   plato_dist_ab        - host buffer of negative_profiles doubles
    //                          containing the raw AB Euclidean distances at the
    //                          plato query (optional; pass nullptr to skip).
    //
    // Optional profile matrices (each column = one query, each row = one rank).
    // All five must be non-null together, or all nullptr to skip.
    //   profile_aa[max_freq × query_count]   - raw AA distances per rank/query
    //   profile_ab[max_freq × query_count]   - raw AB distances per rank/query
    //   profile_rfcp[max_freq × query_count] - clipped contrast per rank/query
    //   profile_index_aa[max_freq × query_count] - 1-based AA neighbor indices
    //   profile_index_ab[max_freq × query_count] - 1-based AB neighbor indices
    RfcpGpuStatus rfcp_gpu_compute(
        const double *positive, const double *mean_positive,
        const double *sig_positive, const int *valid_positive,
        std::uint32_t positive_size,
        const double *negative, const double *mean_negative,
        const double *sig_negative, const int *valid_negative,
        std::uint32_t negative_size,
        std::uint32_t window_size, std::uint32_t max_freq,
        std::uint32_t exclusion_radius,
        std::uint32_t query_begin, std::uint32_t query_end,
        double *rms_profile, int *paired_frequency,
        std::uint32_t *plato_index, double *plato_rms,
        double *plato_dist_aa, double *plato_dist_ab,
        double *profile_aa       = nullptr,
        double *profile_ab       = nullptr,
        double *profile_rfcp     = nullptr,
        int    *profile_index_aa = nullptr,
        int    *profile_index_ab = nullptr);

    bool rfcp_gpu_available();

  } // namespace gpu
} // namespace matrixprofiler

#endif // MATRIXPROFILER_RFCP_GPU_POC_H
