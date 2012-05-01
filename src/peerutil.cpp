/*
 *  io.cpp
 *  peer
 *
 *  Created by Leopold Parts on 13/03/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include "peerutil.h"
#include "sparsefa.h"
#include "tclap/CmdLine.h"
#include "csv_parser/include/csv_parser/csv_parser.hpp"

using namespace TCLAP;
using namespace PEER;

sPeerArgs parseCmdlineArgs(int argc, char * const argv[]){
	sPeerArgs args;
	CmdLine cmd("Probabilistic estimation of expression residuals (PEER)", ' ', "1.0");
	
	try {  
		SwitchArg no_residuals("","no_res_out","No output of residual values", cmd, false);
		SwitchArg no_X("","no_x_out","No output of estimated factors", cmd, false);
		SwitchArg no_W("","no_w_out","No output of estimated factor weights", cmd, false);
		SwitchArg no_Z("","no_z_out","No output of posterior sparsity prior",cmd,false);
		SwitchArg no_Alpha("","no_a_out","No output of weight precision", cmd, false);
		SwitchArg add_mean("","add_mean","Add a covariate to model mean effect", cmd, false);
		SwitchArg has_header("","has_header","Expression and covariates files have a header", cmd, false);
		SwitchArg has_rownames("","has_rownames","Expression files have rownames", cmd, false);
		SwitchArg transpose("","transpose","data is probes*samples (default is rows of samples)", cmd, false);
		
		ValueArg<std::string> out_dir("o","out_dir","Output directory",false,"peer_out","string", cmd);
		ValueArg<std::string> expr_file("f","file","Expression data file",true,"","string", cmd);
		ValueArg<std::string> expr_file_std("","var_file","Expression uncertainty (variance) data file",false,"","string", cmd);
		
		ValueArg<std::string> cov_file("c","cov_file","Covariate data file",false,"","string", cmd);
		ValueArg<std::string> prior_file("","prior","Factor prior file",false,"","string", cmd);
		ValueArg<int> n_factors("n","n_factors","Number of hidden factors",false,5,"int", cmd);
		ValueArg<int> n_iter("i","n_iter","Number of iterations",false,50,"int", cmd);
		ValueArg<float> alpha_pa("","a_pa","Alpha node prior parameter a",false,1.,"float", cmd);
		ValueArg<float> alpha_pb("","a_pb","Alpha node prior parameter b",false,1.,"float", cmd);
		ValueArg<float> eps_pa("","e_pa","Eps node prior parameter a",false,1.,"float", cmd);
		ValueArg<float> eps_pb("","e_pb","Eps node prior parameter b",false,1.,"float", cmd);
		ValueArg<float> bound_tolerance("","bound_tol","Bound tolerance",false,0.001,"float", cmd);
		ValueArg<float> var_tolerance("","var_tol","Variation tolerance",false,0.000001,"float", cmd);
		ValueArg<float> sigma_off("","sigma_off","Variance inactive component",false,0.001,"float", cmd);

		
		cmd.parse( argc, argv );

		args.no_residuals = no_residuals.getValue();
		args.no_X = no_X.getValue();
		args.no_W = no_W.getValue();
		args.no_Alpha = no_Alpha.getValue();
		args.no_Z = no_Z.getValue();
		args.transpose = transpose.getValue();
		args.add_mean = add_mean.getValue();
		args.has_header = has_header.getValue();
		args.has_rownames = has_rownames.getValue();
		args.out_dir = out_dir.getValue();
		args.expr_file = expr_file.getValue();
		args.expr_file_std = expr_file_std.getValue();
		args.cov_file = cov_file.getValue();
		args.prior_file = prior_file.getValue();
		args.n_factors = n_factors.getValue();
		args.n_iter = n_iter.getValue();
		args.alpha_pa = alpha_pa.getValue();
		args.alpha_pb = alpha_pb.getValue();
		args.eps_pa = eps_pa.getValue();
		args.eps_pb = eps_pb.getValue();
		args.bound_tolerance = bound_tolerance.getValue();
		args.var_tolerance = var_tolerance.getValue();
		args.sigma_off = sigma_off.getValue();
	} catch (ArgException &e)  // catch any exceptions
	{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
	
	return args;
}



PMatrix parseCsv(string filename, bool header, bool rownames){
	csv_parser file_parser;	
	vector< vector<float> > data;
	file_parser.set_skip_lines(header);
	file_parser.set_line_term_char('\n');	
	if (filename.rfind(".csv") != string::npos) {
		file_parser.set_field_term_char(',');	
	}
	else if (filename.rfind(".tab") != string::npos) {
		file_parser.set_field_term_char('\t');
	}
	else {
		file_parser.set_field_term_char(' ');		
	}
	file_parser.init(filename.c_str());
	
	int n_rows = 0, n_cols = 0;
	
	while(file_parser.has_more_rows()){
		vector<float> row_data;
		csv_row row = file_parser.get_row();
		
		n_rows++;
		n_cols = row.size();
		for (int i = (int)rownames; i < n_cols; i++) {row_data.push_back(atof(row[i].c_str()));}
		data.push_back(row_data);
	}
	
	// cast result into PMatrix
    n_cols -= (int)rownames;
	PMatrix result = PMatrix::Zero(n_rows, n_cols);
	for (int i = 0; i < n_rows; ++i){
		for (int j = 0; j < n_cols; ++j){
			result(i,j) = data[i][j];
		}
	}
	
	return result;
}


void writeCsv(string filename, PMatrix m){
	ofstream ofh(filename.c_str());
	
	for(int i = 0; i < m.rows(); ++i){
		for(int j = 0; j < m.cols() - 1; ++j){
			ofh << m(i,j) << ",";
		}
		ofh << m(i, m.cols() -1 ) << "\n";
	}
	ofh.close();
}


/**
 Get an instance of PEER FA object from command line arguments
 */
cPEER getInstance(sPeerArgs args){
	
	//main expr matrix
	PMatrix expr = parseCsv(args.expr_file, args.has_header, args.has_rownames);
	if(args.transpose){
		expr.transposeInPlace();
	}

	//optional matrices
	PMatrix expr_std = PMatrix();
	PMatrix covs = PMatrix();
	PMatrix prior = PMatrix();
	if (args.cov_file.length() > 0) covs = parseCsv(args.cov_file, args.has_header, false);
	if (args.expr_file_std.length() > 0) {
        expr_std = parseCsv(args.expr_file_std, args.has_header, false);
        if(args.transpose){
            expr_std.transposeInPlace();
        }
    }
	if (args.prior_file.length() > 0) prior = parseCsv(args.prior_file, args.has_header, false);
	

	
	cSPARSEFA vb;
	
	vb.setPhenoMean(expr);
	vb.setCovariates(covs);
	vb.setSparsityPrior(prior);
	vb.setPhenoVar(expr_std);
	vb.setNk(args.n_factors);
	
	//(expr, PMatrix(), covs, args.n_factors);

	vb.setAdd_mean(args.add_mean);
	vb.setNmax_iterations(args.n_iter);
	vb.setTolerance(args.bound_tolerance);
	vb.setVarTolerance(args.var_tolerance);

	vb.setPriorAlphaA(args.alpha_pa);
	vb.setPriorAlphaB(args.alpha_pb);
	vb.setPriorEpsA(args.eps_pa);
	vb.setPriorEpsB(args.eps_pb);

	return vb;
}


string concat(string str1, string str2){
	char result[1001];
	sprintf(result, "%s%s",str1.c_str(), str2.c_str());
	return string(result);
}


void write_output(cPEER& vb, sPeerArgs args){
    // 1. create out dir
	system(concat("mkdir ", args.out_dir).c_str());

	// 2. output all matrices that are not specified to be skipped
    string out_dir = args.out_dir;
	if (out_dir[strlen(out_dir.c_str()) - 1] != '/') {out_dir = concat(out_dir, "/");}

    if(args.transpose)
    {
        if (!args.no_residuals) writeCsv(concat(out_dir, "residuals.csv"), vb.getResiduals());
        if (!args.no_W) writeCsv(concat(out_dir, "W.csv"), vb.W->E1);
        if (!args.no_X) writeCsv(concat(out_dir, "X.csv"), vb.X->E1);
    }
    else {
        if (!args.no_residuals) writeCsv(concat(out_dir, "residuals.csv"), vb.getResiduals().transpose());
        if (!args.no_W) writeCsv(concat(out_dir, "W.csv"), vb.W->E1.transpose());
        if (!args.no_X) writeCsv(concat(out_dir, "X.csv"), vb.X->E1.transpose());
    }
	
	//if sparse mode on, export indicator Z also
	if(!isnull(vb.getSparsityPrior()))
	{
		if (!args.no_Z) writeCsv(concat(out_dir, "Z.csv"), vb.getZ());
		
	}
	else if (!args.no_Alpha) // if no sparse mode, export posterior weight precision
	{
		writeCsv(concat(out_dir, "Alpha.csv"), vb.Alpha->E1);
	}
}
		
