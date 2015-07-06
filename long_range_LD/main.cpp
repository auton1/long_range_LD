//
//  main.cpp
//  long_range_LD
//
//  Created by Adam Auton on 1/07/15.
//  Copyright (c) 2015 Adam Auton. All rights reserved.
//

#include "main.h"

output_log LOG;

void read_VCF_file(parameters &params, vector<tuple<int, int, double, vector<bool> > > &out_chridx_pos_freq_data, map<int,string> &chridx_to_chr)
{
    variant_file *vf;
    if (!params.bcf_format)
        vf = new vcf_file(params);
    else
        vf = new bcf_file(params);
    
    vf->apply_filters(params);
    LOG.printLOG("After filtering, kept " + output_log::int2str(vf->N_kept_individuals()) + " out of " + output_log::int2str(vf->meta_data.N_indv) + " Individuals\n");
    vf->read_data(params, out_chridx_pos_freq_data, chridx_to_chr);
    LOG.printLOG("After filtering, kept " + header::int2str(vf->N_kept_sites()) + " out of a possible " + header::int2str(vf->N_total_sites()) + " Sites\n");
    if (vf->N_total_sites() <= 0)
        LOG.warning("File does not contain any sites");
    else if (vf->N_kept_sites() <= 0)
        LOG.warning("No data left for analysis!");
    
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
    
    map<int,string> chridx_to_chr;
    vector<tuple<int, int, double, vector<bool> > > chridx_pos_freq_data;
    
    read_VCF_file(params, chridx_pos_freq_data, chridx_to_chr);
    unsigned int N_pos = (unsigned int)chridx_pos_freq_data.size();
    
    int tested_count = 0;
    int passed_count = 0;
    int bound_test = 0;
    
    sort(chridx_pos_freq_data.begin(), chridx_pos_freq_data.end(), [](const tuple<int,int,double, vector<bool>>& a,
                                                      const tuple<int,int,double, vector<bool>>& b) -> bool
                                                    {
                                                        return std::get<2>(a) > std::get<2>(b);
                                                    }); // Sort tuple by third column
    
    
    LOG.printLOG("Taper 2D loop for positive correlations.\n");
    
    cout << "CHR1\tPOS1\tFREQ1\tCHR2\tPOS2\tFREQ2\tPHI" << endl;
    
    int startposi = 1;
    for (unsigned int i=0; i<N_pos-1; i++)
    {
        double A = get<2>(chridx_pos_freq_data[i]);
        vector<bool> data1 = get<3>(chridx_pos_freq_data[i]);
        
        for (unsigned int j=i+1; j<N_pos; j++)
        {
            double B = get<2>(chridx_pos_freq_data[j]);
            if (j >= startposi)
            {
                double upper = sqrt(B*(1-A)/(A*(1-B)));
                bound_test++;
            
                if (upper < params.threshold)
                {   // Pruning by the monotone property
                    if ((j>i+1) || (startposi == (N_pos-1)))
                    {
                        startposi = j;
                    }
                    else
                    {
                        startposi = j+1;
                    }
                    break;
                }
            }
            
            int CHR1idx = get<0>(chridx_pos_freq_data[i]);
            int CHR2idx = get<0>(chridx_pos_freq_data[j]);
            int POS1 = get<1>(chridx_pos_freq_data[i]);
            int POS2 = get<1>(chridx_pos_freq_data[j]);
            
            if ((CHR1idx == CHR2idx) && (abs(POS1 - POS2) < params.min_dist))
            {   // SNPs too close
                continue;
            }
            
            // Refine here...
            vector<bool> data2 = get<3>(chridx_pos_freq_data[j]);

            double C = 0;
            for (unsigned int ui=0; ui<data1.size(); ui++)
                C += (data1[ui] && data2[ui]);
            C = C / data1.size();
            
            double phi = (C - A*B) / sqrt(A*B*(1-A)*(1-B));
            tested_count++;
        
            if (phi >= params.threshold)
            {
                string CHR1 = chridx_to_chr[CHR1idx];
                string CHR2 = chridx_to_chr[CHR2idx];
                
                cout << CHR1 << "\t" << POS1 << "\t" << A << "\t" << CHR2 << "\t" << POS2 << "\t" << B << "\t" << phi << endl;
                passed_count++;
            }
        }
        if (startposi == (i+1))
            startposi++;
    }
    
    LOG.printLOG("Upper Bound Tests " + LOG.int2str(bound_test) + "\n");
    LOG.printLOG("Tested " + LOG.int2str(tested_count) + " candidates out of a possible " + LOG.int2str(N_pos*(N_pos-1)/2) + "\n");
    LOG.printLOG("Found " + LOG.int2str(passed_count) + " sites with Pearson >= " + LOG.dbl2str(params.threshold, 3) + "\n");
    
    /*
    LOG.printLOG("Taper 1D loop.\n");
    tested_count = 0;
    passed_count = 0;
    bound_test = 0;
    for (unsigned int i=0; i<N_pos-1; i++)
    {
        double A = get<2>(chr_pos_freq_data[i]);
        vector<bool> data1 = get<3>(chr_pos_freq_data[i]);
        
        for (unsigned int j=i+1; j<N_pos; j++)
        {
            double B = get<2>(chr_pos_freq_data[j]);
            double upper = sqrt(B*(1-A)/(A*(1-B)));
            bound_test++;
            if (upper < params.threshold)
            {   // Pruning by the monotone property
                break;
            }
            
            // Refine here...
            vector<bool> data2 = get<3>(chr_pos_freq_data[j]);
            
            double C = 0;
            for (unsigned int ui=0; ui<data1.size(); ui++)
                C += (data1[ui] && data2[ui]);
            C = C / data1.size();
            
            double phi = (C - A*B) / sqrt(A*B*(1-A)*(1-B));

            if (phi >= params.threshold)
            {
                //      string CHR1 = get<0>(chr_pos_freq_data[i]);
                //     string CHR2 = get<0>(chr_pos_freq_data[j]);
                //    int POS1 = get<1>(chr_pos_freq_data[i]);
                //   int POS2 = get<1>(chr_pos_freq_data[j]);
                
                // cout << CHR1 << ":" << POS1 << "\t" << CHR2 << ":" << POS2 << "\t" << phi << endl;
                passed_count++;
            }
            
            tested_count++;
        }
    }

    LOG.printLOG("Upper Bound Tests " + LOG.int2str(bound_test) + "\n");
    LOG.printLOG("Tested " + LOG.int2str(tested_count) + " candidates out of a possible " + LOG.int2str(N_pos*(N_pos-1)/2) + "\n");
    LOG.printLOG("Found " + LOG.int2str(passed_count) + " sites with Pearson > " + LOG.dbl2str(params.threshold, 3) + "\n");
    */

    if (params.skip_neg == false)
    {
        LOG.printLOG("TAPER 2D loop for negative correlations.\n");
        tested_count = 0;
        passed_count = 0;
        bound_test=0;
        
        for (unsigned int i=0; i<N_pos-1; i++)
        {
            double A = get<2>(chridx_pos_freq_data[i]);
            vector<bool> data1 = get<3>(chridx_pos_freq_data[i]);
            
            for (unsigned int j=i+1; j<N_pos; j++)
            {
                double B = get<2>(chridx_pos_freq_data[j]);
                double lower;
                
                if ((A+B) > 1.0)
                { // These will be tested first due to the (descending) ordering of the list
                    lower = -sqrt(((1-A)*(1-B))/(A*B));
                    bound_test++;
                    if (lower > -params.threshold)
                    {   // Pruning by the monotone property
                        continue;   // Skip this one.
                    }
                }
                else
                {   // A+B <= 1.0, which will necessarily be tested after the A+B > 1.0 sites.
                    bound_test++;
                    lower = -sqrt((A*B)/((1-A)*(1-B)));
                    if (lower > -params.threshold)
                    {   // Pruning by the monotone property
                        break;  // Finally, no need to continue.
                    }
                }
                
                int CHR1idx = get<0>(chridx_pos_freq_data[i]);
                int CHR2idx = get<0>(chridx_pos_freq_data[j]);
                int POS1 = get<1>(chridx_pos_freq_data[i]);
                int POS2 = get<1>(chridx_pos_freq_data[j]);
                
                if ((CHR1idx == CHR2idx) && (abs(POS1 - POS2) < params.min_dist))
                {   // SNPs too close
                    continue;
                }
                
                // Refine here...
                vector<bool> data2 = get<3>(chridx_pos_freq_data[j]);
                
                double C = 0;
                for (unsigned int ui=0; ui<data1.size(); ui++)
                    C += (data1[ui] && data2[ui]);
                C = C / data1.size();
                
                double phi = (C - A*B) / sqrt(A*B*(1-A)*(1-B));
                tested_count++;

                if (phi <= -params.threshold)
                {
                    int CHR1idx = get<0>(chridx_pos_freq_data[i]);
                    int CHR2idx = get<0>(chridx_pos_freq_data[j]);
                    int POS1 = get<1>(chridx_pos_freq_data[i]);
                    int POS2 = get<1>(chridx_pos_freq_data[j]);
                    
                    string CHR1 = chridx_to_chr[CHR1idx];
                    string CHR2 = chridx_to_chr[CHR2idx];
                    
                    cout << CHR1 << "\t" << POS1 << "\t" << A << "\t" << CHR2 << "\t" << POS2 << "\t" << B << "\t" << phi << endl;
                    passed_count++;
                }
                
            }
        }

        LOG.printLOG("Lower Bound Tests " + LOG.int2str(bound_test) + "\n");
        LOG.printLOG("Tested " + LOG.int2str(tested_count) + " candidates out of a possible " + LOG.int2str(N_pos*(N_pos-1)/2) + "\n");
        LOG.printLOG("Found " + LOG.int2str(passed_count) + " sites with Pearson <= " + LOG.dbl2str(-params.threshold, 3) + "\n");
    }

    time(&end);
    double running_time = difftime(end,start);
    LOG.printLOG("Run Time = " + output_log::dbl2str_fixed(running_time, 2) + " seconds\n");
    LOG.close();
    return 0;
}
