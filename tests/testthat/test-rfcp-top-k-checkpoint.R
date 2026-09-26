test_that("rfcp_top_k produces identical results with chunked execution and resumes correctly", {
  set.seed(123)
  pos <- rnorm(400, sd = 0.05)
  pulse <- sin(seq(0, pi, length.out = 15)) * 2
  for (idx in c(30, 80, 150, 220, 310)) {
    pos[idx:(idx + 14)] <- pos[idx:(idx + 14)] + pulse
  }
  neg <- rnorm(400, sd = 0.05)

  # 1. Execucao de referencia (monolitica)
  ref_res <- rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = 15L,
    k_platos = 2L,
    max_freq = 4L,
    query_block_size = NULL,
    backend = "cpu"
  )

  # 2. Execucao fatiada em blocos pequenos (ex: blocos de 50 queries)
  chunked_res <- rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = 15L,
    k_platos = 2L,
    max_freq = 4L,
    query_block_size = 50L,
    backend = "cpu"
  )

  expect_equal(ref_res$k_found, chunked_res$k_found)
  expect_equal(ref_res$platos[[1]]$plato_index, chunked_res$platos[[1]]$plato_index)
  expect_equal(ref_res$platos[[1]]$plato_rms, chunked_res$platos[[1]]$plato_rms, tolerance = 1e-7)
  if (ref_res$k_found > 1) {
    expect_equal(ref_res$platos[[2]]$plato_index, chunked_res$platos[[2]]$plato_index)
    expect_equal(ref_res$platos[[2]]$plato_rms, chunked_res$platos[[2]]$plato_rms, tolerance = 1e-7)
  }

  # 3. Teste de Retoma Simulada com checkpoint_callback e resume_state
  captured_states <- list()
  checkpoint_cb <- function(st) {
    captured_states[[length(captured_states) + 1L]] <<- st
  }

  # Corre ate gerar pelo menos 2 checkpoints e simula interrupcao
  rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = 15L,
    k_platos = 2L,
    max_freq = 4L,
    query_block_size = 40L,
    checkpoint_callback = checkpoint_cb,
    backend = "cpu"
  )

  expect_true(length(captured_states) >= 2L)
  intermediate_state <- captured_states[[2L]]

  # Retoma a partir do 2o checkpoint
  resumed_res <- rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = 15L,
    k_platos = 2L,
    max_freq = 4L,
    query_block_size = 40L,
    resume_state = intermediate_state,
    backend = "cpu"
  )

  expect_equal(resumed_res$k_found, ref_res$k_found)
  expect_equal(resumed_res$platos[[1]]$plato_index, ref_res$platos[[1]]$plato_index)
  expect_equal(resumed_res$platos[[1]]$plato_rms, ref_res$platos[[1]]$plato_rms, tolerance = 1e-7)
  if (ref_res$k_found > 1) {
    expect_equal(resumed_res$platos[[2]]$plato_index, ref_res$platos[[2]]$plato_index)
    expect_equal(resumed_res$platos[[2]]$plato_rms, ref_res$platos[[2]]$plato_rms, tolerance = 1e-7)
  }
})

test_that("rfcp_top_k emits structured progress events correctly", {
  set.seed(42)
  pos <- rnorm(200, sd = 0.05)
  pos[20:34] <- pos[20:34] + sin(seq(0, pi, length.out = 15)) * 3
  neg <- rnorm(200, sd = 0.05)

  events <- list()
  prog_cb <- function(event) {
    events[[length(events) + 1L]] <<- event
  }

  res <- rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = 15L,
    k_platos = 1L,
    max_freq = 2L,
    query_block_size = 50L,
    progress_callback = prog_cb,
    backend = "cpu"
  )

  event_types <- vapply(events, `[[`, character(1), "event")
  expect_true("round_start" %in% event_types)
  expect_true("block_done" %in% event_types)
  expect_true("round_done" %in% event_types)

  r_start <- events[[which(event_types == "round_start")[[1L]]]]
  expect_equal(r_start$k, 1L)
  expect_equal(r_start$max_freq, 2L)
  expect_equal(r_start$total_queries, 200L - 15L + 1L)
  expect_equal(r_start$n_blocks, ceiling((200L - 15L + 1L) / 50L))

  b_done <- events[[which(event_types == "block_done")[[1L]]]]
  expect_equal(b_done$k, 1L)
  expect_equal(b_done$block_index, 1L)
  expect_equal(b_done$query_begin, 1L)
  expect_equal(b_done$query_end, 50L)
  expect_true(is.numeric(b_done$elapsed_seconds))
  expect_true(is.numeric(b_done$block_plato_rms))

  r_done <- events[[which(event_types == "round_done")[[1L]]]]
  expect_equal(r_done$k, 1L)
  expect_equal(r_done$plato_index, res$platos[[1L]]$plato_index)
  expect_equal(r_done$plato_rms, res$platos[[1L]]$plato_rms)
})

test_that("rfcp_top_k can resume from level 1 checkpoint (between rounds)", {
  set.seed(42)
  pos <- rnorm(300, sd = 0.05)
  sine_pulse <- sin(seq(0, pi, length.out = 15)) * 2
  for (idx in c(30, 80, 150, 220)) {
    pos[idx:(idx + 14)] <- pos[idx:(idx + 14)] + sine_pulse
  }
  neg <- rnorm(300, sd = 0.05)

  ref_res <- rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = 15L,
    k_platos = 2L,
    max_freq = 3L,
    query_block_size = 60L,
    backend = "cpu"
  )

  round_checkpoints <- list()
  ckpt_cb <- function(st) {
    if (st$completed_query_cursor == 0L) {
      round_checkpoints[[length(round_checkpoints) + 1L]] <<- st
    }
  }

  rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = 15L,
    k_platos = 2L,
    max_freq = 3L,
    query_block_size = 60L,
    checkpoint_callback = ckpt_cb,
    backend = "cpu"
  )

  # Checkpoint after round 1 has k_iter = 2L and completed_query_cursor = 0L
  expect_true(length(round_checkpoints) >= 1L)
  r1_state <- round_checkpoints[[1L]]
  expect_equal(r1_state$k_iter, 2L)
  expect_equal(r1_state$count_found, 1L)

  resumed_res <- rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = 15L,
    k_platos = 2L,
    max_freq = 3L,
    query_block_size = 60L,
    resume_state = r1_state,
    backend = "cpu"
  )

  expect_equal(resumed_res$k_found, ref_res$k_found)
  expect_equal(resumed_res$platos[[1]]$plato_index, ref_res$platos[[1]]$plato_index)
  expect_equal(resumed_res$platos[[2]]$plato_index, ref_res$platos[[2]]$plato_index)
})

test_that("rfcp_top_k validates new parameters and resume_state integrity", {
  pos <- rnorm(100)
  neg <- rnorm(100)

  expect_error(
    rfcp_top_k(pos, neg, window_size = 10L, query_block_size = 0L),
    "query_block_size"
  )
  expect_error(
    rfcp_top_k(pos, neg, window_size = 10L, query_block_size = -5L),
    "query_block_size"
  )
  expect_error(
    rfcp_top_k(pos, neg, window_size = 10L, progress_callback = "not_a_function"),
    "progress_callback"
  )
  expect_error(
    rfcp_top_k(pos, neg, window_size = 10L, checkpoint_callback = 42L),
    "checkpoint_callback"
  )
  expect_error(
    rfcp_top_k(pos, neg, window_size = 10L, resume_state = "invalid"),
    "resume_state"
  )

  # Incompatible resume_state (wrong window_size)
  bad_state <- list(
    signature = list(
      n_positive = length(pos),
      n_negative_initial = length(neg),
      window_size = 20L, # Mismatch (actual is 10)
      k_platos = 2L,
      exclusion_radius = 10L,
      adaptive_max_freq = FALSE,
      min_freq = 2L,
      backend = "cpu"
    )
  )
  expect_error(
    rfcp_top_k(pos, neg, window_size = 10L, resume_state = bad_state),
    "resume_state is incompatible with current input parameters."
  )
})
