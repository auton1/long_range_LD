/*
 * variant_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "variant_file.h"

variant_file::~variant_file() {}

void variant_file::read_data(const parameters &params, vector<tuple<int, int, double, vector<bool>> > &out_chridx_pos_freq_data, map<int,string> &chridx_to_chr)
{
    vector<char> variant_line;
    entry *e = get_entry_object();
    istringstream ss;
  
    vector<int> allele_counts;
    unsigned int N_non_missing_chr;
    int chr_idx = 0;
    int current_chr_idx = 0;
    string last_CHROM = "";
    
    while(!eof())
    {
        get_entry(variant_line);
        e->reset(variant_line);
        N_entries += e->apply_filters(params);
        
        if(!e->passed_filters)
            continue;
        N_kept_entries++;
        e->parse_basic_entry(true);
        
        if (e->get_N_alleles() != 2)
        {
            LOG.one_off_warning("\t: Only using biallelic loci.");
            continue;
        }
        
        e->parse_genotype_entries(true);
        
        string CHROM = e->get_CHROM();
        int POS = e->get_POS();
        
        if (CHROM != last_CHROM)
        {
            chridx_to_chr[chr_idx] = CHROM;
            current_chr_idx = chr_idx;
            chr_idx++;
            last_CHROM = CHROM;
            LOG.printLOG("\t" + CHROM + "\n");
        }
        
        e->get_allele_counts(allele_counts, N_non_missing_chr);
        
        if (N_non_missing_chr != e->N_indv*2)
        {
            LOG.one_off_warning("\t: Only using sites with no missing data.");
            continue;
        }
        
        double freq = allele_counts[1]/(double)N_non_missing_chr;
        
        pair<int,int> genotype;
        vector<bool> data(N_non_missing_chr);
        for (unsigned int ui=0; ui<e->N_indv; ui++)
        {
            e->get_indv_GENOTYPE_ids(ui, genotype);
            data[2*ui]=(bool)genotype.first;
            data[(2*ui)+1]=(bool)genotype.second;
        }
        
        out_chridx_pos_freq_data.push_back(make_tuple(current_chr_idx, POS, freq, data));
    }
    delete e;
}


void variant_file::read_PL_data(const parameters &params, int GL_or_PL, vector<int> &out_pos, vector<string> &ref, vector<string> &alt, vector< vector<vector<double> > > &out_matrix)
{
	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in VCF file in order to output BEAGLE genotype likelihoods.");

	if (GL_or_PL == 0)
		LOG.printLOG("Getting GLs\n");
	else if (GL_or_PL == 1)
		LOG.printLOG("Getting PLs\n");
	else
		LOG.error("Unknown GL or PL option.");

	string GL_entry, tmp_string;
	vector<char> variant_line;
	entry *e = get_entry_object();
	double lk1, lk2, lk3;
	bool found_GL=false;
	istringstream ss;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\t: Only using biallelic loci.");
			continue;
		}

		e->parse_full_entry(true);

		if (GL_or_PL == 0)
			if (e->FORMAT_id_exists("GL") == false)
				continue;
		if (GL_or_PL == 1)
			if (e->FORMAT_id_exists("PL") == false)
				continue;

		found_GL = true;

		string CHROM = e->get_CHROM();
		int POS = e->get_POS();
		string a1 = e->get_REF();
		string a2 = e->get_ALT();

		vector<vector<double> > all_indv_PL;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			vector<double> indv_PL(3, 1.0);
			if (include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == true)
			{
				if (GL_or_PL == 0)
					e->read_indv_generic_entry(ui, "GL", GL_entry);
				else
					e->read_indv_generic_entry(ui, "PL", GL_entry);
				ss.clear();
				ss.str(GL_entry);
				getline(ss, tmp_string, ',');
				lk1 = atof(tmp_string.c_str());
				getline(ss, tmp_string, ',');
				lk2 = atof(tmp_string.c_str());
				getline(ss, tmp_string);
				lk3 = atof(tmp_string.c_str());
				if (GL_or_PL == 0)
				{
					indv_PL[0] = lk1*log(10.0);
					indv_PL[1] = lk2*log(10.0);
					indv_PL[2] = lk3*log(10.0);
				}
				else
				{
					indv_PL[0] = -0.1*lk1*log(10.0);
					indv_PL[1] = -0.1*lk2*log(10.0);
					indv_PL[2] = -0.1*lk3*log(10.0);
				}
			}
			else
			{	// Mark as unknown
				indv_PL[0] = 1.0;
				indv_PL[1] = 1.0;
				indv_PL[2] = 1.0;
			}
			all_indv_PL.push_back(indv_PL);
		}
		out_matrix.push_back(all_indv_PL);
		out_pos.push_back(POS);
		ref.push_back(a1);
		alt.push_back(a2);
	}
	delete e;
	if (found_GL == false)
		LOG.error("Require GL or PL FORMAT tags in VCF file to output BEAGLE input.");
}

// Return the number of individuals that have not been filtered out
int variant_file::N_kept_individuals() const
{
	int N_kept = 0;
	for (unsigned int ui=0; ui<include_indv.size(); ui++)
		if (include_indv[ui] == true)
			N_kept++;
	return N_kept;
}

// Return the number of sites that have not been filtered out
int variant_file::N_kept_sites() const
{
	return N_kept_entries;
}

// Return the total number of sites in the file
int variant_file::N_total_sites() const
{
	return N_entries;
}

void variant_file::ByteSwap(unsigned char *b, int n) const
{
   int i = 0;
   int j = n-1;
   while (i<j)
   {
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

void variant_file::get_contigs(const string &contigs_file, vector<string> &contig_vector)
{
	if (contigs_file == "")
		LOG.error("Contig declarations in header are necessary for BCF conversion. Use --contigs <filename> to add contigs to the header.");

	ifstream contigs(contigs_file.c_str());
	if (!contigs.is_open())
		LOG.error("Could not open contigs file: " + contigs_file);

	string line;
	int contig_lines = 0;
	contig_vector.resize(0);

	while (getline(contigs, line))
	{
		if (line.find("##contig=")==string::npos)
			LOG.error("Contigs file must contain only contig header lines.");

		contig_vector.push_back(line);
		contig_lines++;
	}

	contigs.close();
	LOG.printLOG("Including "+header::int2str(contig_lines)+" header lines from the contig file.\n");
}
