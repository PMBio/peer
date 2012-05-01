/*
 *  io.h
 *  peer
 *
 *  Created by Leopold Parts on 13/03/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PEERUTIL
#define PEERUTIL
#include "vbfa.h"
#include "sparsefa.h"

using namespace PEER;

struct sPeerArgs{
	bool no_residuals, no_X, no_W, no_Alpha,no_Z, add_mean, has_header, transpose, has_rownames;
	string out_dir, expr_file, expr_file_std, cov_file, prior_file;
	int n_factors, n_iter;
	float alpha_pa, alpha_pb, eps_pa, eps_pb, bound_tolerance, var_tolerance, sigma_off;
};


sPeerArgs parseCmdlineArgs(int argc, char * const argv[]);

PMatrix parseCsv(const char *filename, bool header);

cPEER getInstance(sPeerArgs args);

void write_output(cPEER& vb, sPeerArgs args);


#endif
