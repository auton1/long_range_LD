//
//  main.cpp
//  long_range_LD
//
//  Created by Adam Auton on 1/07/15.
//  Copyright (c) 2015 Adam Auton. All rights reserved.
//

#include "main.h"

output_log LOG;

void read_VCF_file(parameters &params, vector<tuple<string, int, double, vector<bool> > > &out_chr_pos_freq_data)
{
    variant_file *vf;
    if (!params.bcf_format)
        vf = new vcf_file(params);
    else
        vf = new bcf_file(params);
    
    vf->apply_filters(params);
    LOG.printLOG("After filtering, kept " + output_log::int2str(vf->N_kept_individuals()) + " out of " + output_log::int2str(vf->meta_data.N_indv) + " Individuals\n");
    vf->read_data(params, out_chr_pos_freq_data);
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
    
    vector<tuple<string, int, double, vector<bool> > > chr_pos_freq_data;
    
    read_VCF_file(params, chr_pos_freq_data);
    unsigned int N_pos = (unsigned int)chr_pos_freq_data.size();
    
    int tested_count = 0;
    int passed_count = 0;
    
    sort(chr_pos_freq_data.begin(), chr_pos_freq_data.end(), [](const tuple<string,int,double, vector<bool>>& a,
                                                      const tuple<string,int,double, vector<bool>>& b) -> bool
                                                    {
                                                        return std::get<2>(a) > std::get<2>(b);
                                                    }); // Sort tuple by third column
    
    
    LOG.printLOG("First loop.\n");
    int startposi = 1;
    double threshold = params.threshold*params.threshold;
    bool flag=true;
    
    for (unsigned int i=0; i<N_pos-1; i++)
    {
        double A = get<2>(chr_pos_freq_data[i]);
        vector<bool> data1 = get<3>(chr_pos_freq_data[i]);
        
        for (unsigned int j=i+1; j<N_pos; j++)
        {
            flag = false;
            double B = get<2>(chr_pos_freq_data[j]);
            if (j >= startposi)
            {
                double upper = sqrt(B*(1-A)/(A*(1-B)));
                //double upper2 = (B*(1-A)/(A*(1-B)));
            
                //if (upper2 < threshold)
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
            
            // Refine here...
            vector<bool> data2 = get<3>(chr_pos_freq_data[j]);

            double C = 0;
            for (unsigned int ui=0; ui<data1.size(); ui++)
                C += (data1[ui] && data2[ui]);
            C = C / data1.size();
            
            double phi = (C - A*B) / sqrt(A*B*(1-A)*(1-B));
            //double phi2 = (C - A*B)*(C - A*B) / (A*B*(1-A)*(1-B));
            
            //if ((phi2 >= threshold))
            if ((phi >= params.threshold))
            {
          //      string CHR1 = get<0>(chr_pos_freq_data[i]);
           //     string CHR2 = get<0>(chr_pos_freq_data[j]);
            //    int POS1 = get<1>(chr_pos_freq_data[i]);
             //   int POS2 = get<1>(chr_pos_freq_data[j]);
                
               // cout << CHR1 << ":" << POS1 << "\t" << CHR2 << ":" << POS2 << "\t" << phi << endl;
                passed_count++;
            }
            
            tested_count++;
            
            if ((startposi == (i+1)) && (flag == false))
                startposi++;
        }
    }
    
    LOG.printLOG("Tested " + LOG.int2str(tested_count) + " candidates out of a possible " + LOG.int2str(N_pos*(N_pos-1)/2) + "\n");
    LOG.printLOG("Found " + LOG.int2str(passed_count) + " sites with Pearson > " + LOG.dbl2str(params.threshold, 3) + "\n");
    
     
    tested_count = 0;
    passed_count = 0;
    for (unsigned int i=0; i<N_pos-1; i++)
    {
        double A = get<2>(chr_pos_freq_data[i]);
        vector<bool> data1 = get<3>(chr_pos_freq_data[i]);
        
        for (unsigned int j=i+1; j<N_pos; j++)
        {
            double B = get<2>(chr_pos_freq_data[j]);
            double upper = sqrt(B*(1-A)/(A*(1-B)));
            
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
            //double phi2 = (C - A*B)*(C - A*B) / (A*B*(1-A)*(1-B));
            
            //if ((phi2 >= threshold))
            if ((phi >= params.threshold))
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

    
    LOG.printLOG("Tested " + LOG.int2str(tested_count) + " candidates out of a possible " + LOG.int2str(N_pos*(N_pos-1)/2) + "\n");
    LOG.printLOG("Found " + LOG.int2str(passed_count) + " sites with Pearson > " + LOG.dbl2str(params.threshold, 3) + "\n");
    
    //reverse(chr_pos_freq_data.begin(), chr_pos_freq_data.end());
    exit(0);
    // The following is currently broken.
    
    tested_count = 0;
    passed_count = 0;
    LOG.printLOG("Second loop.\n");
    for (unsigned int i=0; i<N_pos-1; i++)
    {
        double A = get<2>(chr_pos_freq_data[i]);
        vector<bool> data1 = get<3>(chr_pos_freq_data[i]);
        
        for (unsigned int j=i+1; j<N_pos; j++)
        {
            double B = get<2>(chr_pos_freq_data[j]);
            
            double lower = 0;
            if (A+B <= 1)
            {
               lower = -sqrt(A*B / ((1-A)*(1-B)));
            }
            else
            {
                lower = -sqrt((1-A)*(1-B)/(A*B));
            }
                
            if ((A>=B) && (lower > -params.threshold))  // Not sure about this.
            {   // Pruning by the monotone property
                continue;
            }
            else if (lower > -params.threshold)
            {
                break;
            }
            
            // Refine here...
            vector<bool> data2 = get<3>(chr_pos_freq_data[j]);
            
            double C = 0;
            for (unsigned int ui=0; ui<data1.size(); ui++)
                C += (data1[ui] && data2[ui]);
            C = C / data1.size();
            
            double phi = (C - A*B) / sqrt(A*B*(1-A)*(1-B));
            
            if ((phi <= -params.threshold))
            {
           //     string CHR1 = get<0>(chr_pos_freq_data[i]);
            //    string CHR2 = get<0>(chr_pos_freq_data[j]);
            //    int POS1 = get<1>(chr_pos_freq_data[i]);
            //    int POS2 = get<1>(chr_pos_freq_data[j]);
                
               // cout << CHR1 << ":" << POS1 << "\t" << CHR2 << ":" << POS2 << "\t" << phi << endl;
                passed_count++;
            }
            tested_count++;
        }
    }
    
    LOG.printLOG("Tested " + LOG.int2str(tested_count) + " candidates out of a possible " + LOG.int2str(N_pos*(N_pos-1)/2) + "\n");
    LOG.printLOG("Found " + LOG.int2str(passed_count) + " sites with Pearson > " + LOG.dbl2str(params.threshold, 3) + "\n");

    time(&end);
    double running_time = difftime(end,start);
    LOG.printLOG("Run Time = " + output_log::dbl2str_fixed(running_time, 2) + " seconds\n");
    LOG.close();
    return 0;
}
