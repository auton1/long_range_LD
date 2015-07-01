//
//  main.cpp
//  long_range_LD
//
//  Created by Adam Auton on 1/07/15.
//  Copyright (c) 2015 Adam Auton. All rights reserved.
//

#include "main.h"

output_log LOG;


void read_VCF_file(parameters &params, vector<string> &out_chr, vector<int> &out_pos, vector<string> &ref, vector<string> &alt, vector< double > &out_freq, vector<string> &indv)
{
    variant_file *vf;
    if (!params.bcf_format)
        vf = new vcf_file(params);
    else
        vf = new bcf_file(params);
    
    vf->apply_filters(params);
    LOG.printLOG("After filtering, kept " + output_log::int2str(vf->N_kept_individuals()) + " out of " + output_log::int2str(vf->meta_data.N_indv) + " Individuals\n");
    //vf->read_PL_data(params, 1, out_pos, ref, alt, out_PLs);
    LOG.printLOG("After filtering, kept " + header::int2str(vf->N_kept_sites()) + " out of a possible " + header::int2str(vf->N_total_sites()) + " Sites\n");
    if (vf->N_total_sites() <= 0)
        LOG.warning("File does not contain any sites");
    else if (vf->N_kept_sites() <= 0)
        LOG.warning("No data left for analysis!");
    
    indv = (vf->meta_data).indv;
    delete vf;
}


int main(int argc, char *argv[])
{

    time_t start, end;
    time(&start);
    
    // The following turns off sync between C and C++ streams.
    // Apparently it's faster to turn sync off, and as I don't use C streams, it's okay to turn off.
    ios_base::sync_with_stdio(false);
    
    parameters params(argc, argv);
    params.print_help();
    params.read_parameters();
    
    LOG.open(params.stream_out, params.stream_err, params.output_prefix);
    
    LOG.printLOG("\nlong_range_LD - " + VCFTOOLS_VERSION + "\n");
    LOG.printLOG("(C) Adam Auton 2015\n\n");
    
    params.print_params();
    
    vector<string> chr;
    vector<int> pos;
    vector<string> ref, alt;
    vector< double > freq;
    vector<string> indv;
    read_VCF_file(params, chr, pos, ref, alt, freq, indv);
    int N_pos = pos.size();
    
    
    time(&end);
    double running_time = difftime(end,start);
    LOG.printLOG("Run Time = " + output_log::dbl2str_fixed(running_time, 2) + " seconds\n");
    LOG.close();
    return 0;
}
