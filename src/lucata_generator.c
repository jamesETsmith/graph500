/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */
/*           Anton Korzh                                                   */

/* These need to be before any possible inclusions of stdint.h or inttypes.h.
 * */

/*
 __          __     _____  _   _ _____ _   _  _____
 \ \        / /\   |  __ \| \ | |_   _| \ | |/ ____|
  \ \  /\  / /  \  | |__) |  \| | | | |  \| | |  __
   \ \/  \/ / /\ \ |  _  /| . ` | | | | . ` | | |_ |
    \  /\  / ____ \| | \ \| |\  |_| |_| |\  | |__| |
     \/  \/_/    \_\_|  \_\_| \_|_____|_| \_|\_____|

This was ripped from the graph500 generator to dump out graphs with a header
that beedrill and graph_analyzer can process.

There are A LOT of code paths that don't make sense anymore. Since we're in a
rush right now, I haven't cleaned them up yet.

You should be able to run the file with input like this:

  SEED0=19 TMPFILE=matrix.bin ./lucata_generator 10


James Smith (2023-10-20)
*/

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../generator/make_graph.h"
#include "../generator/utils.h"
#include "aml.h"
#include "common.h"

struct lucata_settings_t {
  size_t filesize;
  size_t header_size;
  char header[1024];

} lucata_settings;

int main(int argc, char** argv) {
  aml_init(&argc, &argv);  // includes MPI_Init inside
  setup_globals();

  /* Parse arguments. */
  int SCALE = 16;
  int edgefactor = 16; /* nedges / nvertices, i.e., 2*avg. degree */
  if (argc >= 2) SCALE = atoi(argv[1]);
  if (argc >= 3) edgefactor = atoi(argv[2]);
  if (argc <= 1 || argc >= 4 || SCALE == 0 || edgefactor == 0) {
    if (rank == 0) {
      fprintf(stderr,
              "Usage: %s SCALE edgefactor\n  SCALE = log_2(# vertices) "
              "[integer, required]\n  edgefactor = (# edges) / (# vertices) = "
              ".5 * (average vertex degree) [integer, defaults to 16]\n(Random "
              "number seed and Kronecker initiator are in main.c)\n",
              argv[0]);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  uint64_t seed1 = 2, seed2 = 3;
  char* tmp;
  if (getenv("SEED0")) {
    seed1 = strtoul(getenv("SEED0"), &tmp, 10);
  }
  if (getenv("SEED1")) {
    seed2 = strtoul(getenv("SEED1"), &tmp, 10);
  }

  const char* filename = getenv("TMPFILE");
  const int reuse_file = getenv("REUSEFILE") ? 1 : 0;
  /* If filename is NULL, store data in memory */

  tuple_graph tg;
  tg.nglobaledges = (int64_t)(edgefactor) << SCALE;
  int64_t nglobalverts = (int64_t)(1) << SCALE;

  tg.data_in_file = (filename != NULL);
  tg.write_file = 1;

  if (tg.data_in_file) {
    int is_opened = 0;
    int mode = MPI_MODE_RDWR | MPI_MODE_EXCL | MPI_MODE_UNIQUE_OPEN;

    if (!reuse_file) {
      mode |= MPI_MODE_CREATE;
      //   mode |= MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE;
      if (rank == 0) {
        if (access(filename, F_OK) == 0) {
          fprintf(stdout, "deleting %s\n", filename);

          if (remove(filename)) {
            fprintf(stdout, "couldn't delete file\n");
          }
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    } else {
      MPI_File_set_errhandler(MPI_FILE_NULL, MPI_ERRORS_RETURN);

      if (MPI_File_open(MPI_COMM_WORLD, (char*)filename, mode, MPI_INFO_NULL,
                        &tg.edgefile)) {
        if (0 == rank && getenv("VERBOSE"))
          fprintf(stderr, "%d: failed to open %s, creating\n", rank, filename);

        mode |= MPI_MODE_RDWR | MPI_MODE_CREATE;
      } else {
        MPI_Offset size;
        MPI_File_get_size(tg.edgefile, &size);
        if (size == tg.nglobaledges * sizeof(packed_edge)) {
          is_opened = 1;
          tg.write_file = 0;
        } else /* Size doesn't match, assume different parameters. */
          MPI_File_close(&tg.edgefile);
      }
    }
    MPI_File_set_errhandler(MPI_FILE_NULL, MPI_ERRORS_ARE_FATAL);
    if (!is_opened) {
      MPI_File_open(MPI_COMM_WORLD, (char*)filename, mode, MPI_INFO_NULL,
                    &tg.edgefile);
      //   MPI_File_set_size(tg.edgefile, tg.nglobaledges *
      //   sizeof(packed_edge));
      MPI_File_set_size(tg.edgefile, lucata_settings.filesize);
    }
    // MPI_File_set_view(tg.edgefile, 0, packed_edge_mpi_type,
    //                   packed_edge_mpi_type, "native", MPI_INFO_NULL);
    if (rank == 0) {
      size_t sz =
          snprintf(lucata_settings.header, 1023,
                   "--format el64 --num_edges %lu --num_vertices "
                   "%lu --is_undirected --seed0 %lu --seed1 %lu\n",
                   tg.nglobaledges, tg.nglobaledges / edgefactor, seed1, seed2);
      lucata_settings.filesize = tg.nglobaledges * sizeof(packed_edge) + sz;
      lucata_settings.header_size = sz;
      printf("HERE'S the HEADER:\n%s", lucata_settings.header);
      MPI_File_write(tg.edgefile, lucata_settings.header,
                     lucata_settings.header_size, MPI_CHAR, MPI_STATUS_IGNORE);
      //   MPI_File_write_at(tg.edgefile, start_edge_index, actual_buf,
      //   edge_count,
      // packed_edge_mpi_type, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_set_view(tg.edgefile, 0, packed_edge_mpi_type,
                      packed_edge_mpi_type, "native", MPI_INFO_NULL);
    MPI_File_set_atomicity(tg.edgefile, 0);
  }  // end if (tg.data_in_file)

  /* Make the raw graph edges. */
  /* Get roots for BFS runs, plus maximum vertex with non-zero degree (used by
   * validator). */
  int num_bfs_roots = 64;
  int64_t* bfs_roots = (int64_t*)xmalloc(num_bfs_roots * sizeof(int64_t));

  double make_graph_start = MPI_Wtime();
  if (!tg.data_in_file || tg.write_file) {
    /* Spread the two 64-bit numbers into five nonzero values in the correct
     * range. */
    uint_fast32_t seed[5];
    make_mrg_seed(seed1, seed2, seed);

    /* As the graph is being generated, also keep a bitmap of vertices with
     * incident edges.  We keep a grid of processes, each row of which has a
     * separate copy of the bitmap (distributed among the processes in the
     * row), and then do an allreduce at the end.  This scheme is used to avoid
     * non-local communication and reading the file separately just to find BFS
     * roots. */
    MPI_Offset nchunks_in_file =
        (tg.nglobaledges + FILE_CHUNKSIZE - 1) / FILE_CHUNKSIZE;
    int64_t bitmap_size_in_bytes =
        int64_min(BITMAPSIZE, (nglobalverts + CHAR_BIT - 1) / CHAR_BIT);
    if (bitmap_size_in_bytes * size * CHAR_BIT < nglobalverts) {
      bitmap_size_in_bytes =
          (nglobalverts + size * CHAR_BIT - 1) / (size * CHAR_BIT);
    }
    int ranks_per_row = tg.data_in_file
                            ? ((nglobalverts + CHAR_BIT - 1) / CHAR_BIT +
                               bitmap_size_in_bytes - 1) /
                                  bitmap_size_in_bytes
                            : 1;
    int nrows = size / ranks_per_row;
    int my_row = -1, my_col = -1;
    MPI_Comm cart_comm;
    {
      int dims[2] = {size / ranks_per_row, ranks_per_row};
      int periods[2] = {0, 0};
      MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
    }
    int in_generating_rectangle = 0;
    if (cart_comm != MPI_COMM_NULL) {
      in_generating_rectangle = 1;
      {
        int dims[2], periods[2], coords[2];
        MPI_Cart_get(cart_comm, 2, dims, periods, coords);
        my_row = coords[0];
        my_col = coords[1];
      }
      MPI_Comm this_col;
      MPI_Comm_split(cart_comm, my_col, my_row, &this_col);
      MPI_Comm_free(&cart_comm);
      /* Every rank in a given row creates the same vertices (for updating the
       * bitmap); only one writes them to the file (or final memory buffer). */
      packed_edge* buf =
          (packed_edge*)xmalloc(FILE_CHUNKSIZE * sizeof(packed_edge));
      MPI_Offset block_limit = (nchunks_in_file + nrows - 1) / nrows;
      /* fprintf(stderr, "%d: nchunks_in_file = %" PRId64 ", block_limit = %"
       * PRId64 " in grid of %d rows, %d cols\n", rank,
       * (int64_t)nchunks_in_file, (int64_t)block_limit, nrows, ranks_per_row);
       */
      if (tg.data_in_file) {
        tg.edgememory_size = 0;
        tg.edgememory = NULL;
      } else {
        int my_pos = my_row + my_col * nrows;
        int last_pos =
            (tg.nglobaledges %
                 ((int64_t)FILE_CHUNKSIZE * nrows * ranks_per_row) !=
             0)
                ? (tg.nglobaledges / FILE_CHUNKSIZE) % (nrows * ranks_per_row)
                : -1;
        int64_t edges_left = tg.nglobaledges % FILE_CHUNKSIZE;
        int64_t nedges =
            FILE_CHUNKSIZE * (tg.nglobaledges / ((int64_t)FILE_CHUNKSIZE *
                                                 nrows * ranks_per_row)) +
            FILE_CHUNKSIZE * (my_pos < (tg.nglobaledges / FILE_CHUNKSIZE) %
                                           (nrows * ranks_per_row)) +
            (my_pos == last_pos ? edges_left : 0);
        /* fprintf(stderr, "%d: nedges = %" PRId64 " of %" PRId64 "\n", rank,
         * (int64_t)nedges, (int64_t)tg.nglobaledges); */
        tg.edgememory_size = nedges;
        tg.edgememory = (packed_edge*)xmalloc(nedges * sizeof(packed_edge));
      }
      MPI_Offset block_idx;
      for (block_idx = 0; block_idx < block_limit; ++block_idx) {
        /* fprintf(stderr, "%d: On block %d of %d\n", rank, (int)block_idx,
         * (int)block_limit); */
        MPI_Offset start_edge_index = int64_min(
            FILE_CHUNKSIZE * (block_idx * nrows + my_row), tg.nglobaledges);
        MPI_Offset edge_count =
            int64_min(tg.nglobaledges - start_edge_index, FILE_CHUNKSIZE);
        packed_edge* actual_buf =
            (!tg.data_in_file && block_idx % ranks_per_row == my_col)
                ? tg.edgememory + FILE_CHUNKSIZE * (block_idx / ranks_per_row)
                : buf;
        /* fprintf(stderr, "%d: My range is [%" PRId64 ", %" PRId64 ") %swriting
         * into index %" PRId64 "\n", rank, (int64_t)start_edge_index,
         * (int64_t)(start_edge_index + edge_count), (my_col == (block_idx %
         * ranks_per_row)) ? "" : "not ", (int64_t)(FILE_CHUNKSIZE * (block_idx
         * / ranks_per_row))); */
        if (!tg.data_in_file && block_idx % ranks_per_row == my_col) {
          assert(FILE_CHUNKSIZE * (block_idx / ranks_per_row) + edge_count <=
                 tg.edgememory_size);
        }
        if (tg.write_file) {
          generate_kronecker_range(seed, SCALE, start_edge_index,
                                   start_edge_index + edge_count, actual_buf);
          if (tg.data_in_file &&
              my_col ==
                  (block_idx %
                   ranks_per_row)) { /* Try to spread writes among ranks */
            MPI_File_write_at(tg.edgefile, start_edge_index, actual_buf,
                              edge_count, packed_edge_mpi_type,
                              MPI_STATUS_IGNORE);
          }
        } else {
          /* All read rather than syncing up for a row broadcast. */
          MPI_File_read_at(tg.edgefile, start_edge_index, actual_buf,
                           edge_count, packed_edge_mpi_type, MPI_STATUS_IGNORE);
        }
      }
      free(buf);
      MPI_Comm_free(&this_col);
    } else {
      tg.edgememory = NULL;
      tg.edgememory_size = 0;
    }
    MPI_Allreduce(&tg.edgememory_size, &tg.max_edgememory_size, 1, MPI_INT64_T,
                  MPI_MAX, MPI_COMM_WORLD);
    if (tg.data_in_file && tg.write_file) {
      MPI_File_sync(tg.edgefile);
    }
  }

  free(bfs_roots);

  if (tg.data_in_file) {
    MPI_File_close(&tg.edgefile);
  } else {
    free(tg.edgememory);
    tg.edgememory = NULL;
  }

  cleanup_globals();
  aml_finalize();  // includes MPI_Finalize()
  return 0;
}
