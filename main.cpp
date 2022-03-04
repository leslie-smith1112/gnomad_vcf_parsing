#include <iostream>
#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>
#include <hts.h>
#include <vcf.h>
#include <vector>


int main(int argc, char **argv)
{

    CLI::App app("Frequencies reader");
    std::string vcf_file_name;
    std::string out_file;
    app.add_option("-i,--input", vcf_file_name, "Input VCF file")->required();
    app.add_option("-o,--output", out_file, "Output file")->required();
    CLI11_PARSE(app, argc, argv);

    htsFile* inf = bcf_open(vcf_file_name.c_str(), "r");
    if (inf == NULL)
    {
        spdlog::error("Can't open vcf file: {}", vcf_file_name);
        std::exit(EXIT_FAILURE);
    }

    spdlog::info("Parsing vcf: {}", vcf_file_name);
 
    bcf_hdr_t* hdr = bcf_hdr_read(inf);

    bcf1_t* rec = bcf_init();
    if (rec == NULL)
    {
        spdlog::error("Error while parsing vcf file: {}", vcf_file_name);
        bcf_close(inf);
        bcf_hdr_destroy(hdr);
        std::exit(EXIT_FAILURE);
    }


    std::ofstream out_file_stream(out_file);
    if (not out_file_stream.is_open()) { spdlog::error("Can't open output file {}", out_file); }
    else { spdlog::info("Writing output to {}", out_file); }
    
    while (bcf_read(inf, hdr, rec) == 0)
    {
    	int chrom = rec->rid;
        int position = rec->pos;
        const char* chr_name = bcf_hdr_id2name(hdr, rec->rid);
	bcf_unpack(rec,BCF_UN_STR);
        bcf_dec_t qhat = rec ->d;
        char* ID = qhat.id;
        std::string someString(ID);
        out_file_stream << chr_name << "\t"<< position << "\t"<< ID << "\t";

	std::string line;
        bcf_unpack(rec, BCF_UN_ALL);
        bcf_info_t* AN  = bcf_get_info(hdr, rec, "AN");
        int ANT = AN ->v1.i;

	bcf_info_t* amr_male  = bcf_get_info(hdr, rec, "AN_amr_male");
        int AN_amr_male = amr_male ->v1.i;
        float prob_amr_male = (float) AN_amr_male/ (float) ANT;
        out_file_stream << "AN_amr_male=" << prob_amr_male << ";";

	bcf_info_t* amr_female  = bcf_get_info(hdr, rec, "AN_amr_female");
        int AN_amr_female = amr_female ->v1.i;
        float prob_amr_female = (float) AN_amr_female/ (float) ANT;
        out_file_stream << "AN_amr_female=" << prob_amr_female << ";";

        bcf_info_t* eas_male  = bcf_get_info(hdr, rec, "AN_eas_male");
        int AN_eas_male = eas_male ->v1.i;
        float prob_eas_male = (float) AN_eas_male/ (float) ANT;
        out_file_stream << "AN_eas_male=" << prob_eas_male << ";";

        bcf_info_t* eas_female  = bcf_get_info(hdr, rec, "AN_eas_female");
        int AN_eas_female = eas_female ->v1.i;
        float prob_eas_female = (float) AN_eas_female/ (float) ANT;
        out_file_stream << "AN_eas_female=" << prob_eas_female << ";";

        bcf_info_t* nfe_male  = bcf_get_info(hdr, rec, "AN_nfe_male");
        int AN_nfe_male = nfe_male ->v1.i;
        float prob_nfe_male = (float) AN_nfe_male/ (float) ANT;
        out_file_stream << "AN_nfe_male=" << prob_nfe_male << ";";

        bcf_info_t* nfe_female  = bcf_get_info(hdr, rec, "AN_nfe_female");
        int AN_nfe_female = nfe_female ->v1.i;
        float prob_nfe_female = (float) AN_nfe_female/ (float) ANT;
        out_file_stream << "AN_nfe_female=" << prob_nfe_female << ";";

        bcf_info_t* asj_male  = bcf_get_info(hdr, rec, "AN_asj_male");
        int AN_asj_male = asj_male ->v1.i;
        float prob_asj_male = (float) AN_asj_male/ (float) ANT;
        out_file_stream << "AN_asj_male=" << prob_asj_male << ";";

        bcf_info_t* asj_female  = bcf_get_info(hdr, rec, "AN_asj_female");
        int AN_asj_female = asj_female ->v1.i;
        float prob_asj_female = (float) AN_asj_female/ (float) ANT;
        out_file_stream << "AN_asj_female=" << prob_asj_female << ";";

        bcf_info_t* oth_male  = bcf_get_info(hdr, rec, "AN_oth_male");
        int AN_oth_male = oth_male ->v1.i;
        float prob_oth_male = (float) AN_oth_male/ (float) ANT;
        out_file_stream << "AN_oth_male=" << prob_oth_male << ";";

        bcf_info_t* oth_female  = bcf_get_info(hdr, rec, "AN_oth_female");
        int AN_oth_female = oth_female ->v1.i;
        float prob_oth_female = (float) AN_oth_female/ (float) ANT;
        out_file_stream << "AN_oth_female=" << prob_oth_female << ";";

        bcf_info_t* fin_male  = bcf_get_info(hdr, rec, "AN_fin_male");
        int AN_fin_male = fin_male ->v1.i;
        float prob_fin_male = (float) AN_fin_male/ (float) ANT;
        out_file_stream << "AN_fin_male=" << prob_fin_male << ";";

        bcf_info_t* fin_female  = bcf_get_info(hdr, rec, "AN_fin_female");
        int AN_fin_female = fin_female ->v1.i;
        float prob_fin_female = (float) AN_fin_female/ (float) ANT;
        out_file_stream << "AN_fin_female=" << prob_fin_female << ";";

        bcf_info_t* afr_male  = bcf_get_info(hdr, rec, "AN_afr_male");
        int AN_afr_male = afr_male ->v1.i;
        float prob_afr_male = (float) AN_afr_male/ (float) ANT;
        out_file_stream << "AN_afr_male=" << prob_afr_male << ";";

        bcf_info_t* afr_female  = bcf_get_info(hdr, rec, "AN_afr_female");
        int AN_afr_female = afr_female ->v1.i;
        float prob_afr_female = (float) AN_afr_female/ (float) ANT;
        out_file_stream << "AN_afr_female=" << prob_afr_female << ";";

        out_file_stream << std::endl;
    }



    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    bcf_destroy(rec);
}













